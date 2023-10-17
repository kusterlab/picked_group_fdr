from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt, QObject, QThread, pyqtSignal # https://realpython.com/python-pyqt-qthread/

import sys

# mock matplotlib imports
sys.modules['matplotlib'] = __import__('matplotlib_mock')

import os
import logging
from logging.handlers import QueueListener
import time
import datetime
import multiprocessing

import mokapot
import numpy as np
from joblib import parallel_backend

import picked_group_fdr.pipeline.andromeda2pin as andromeda2pin
import picked_group_fdr.pipeline.update_evidence_from_pout as update_evidence
import picked_group_fdr.pipeline.merge_pout as merge_pout
import picked_group_fdr.picked_group_fdr as picked_group_fdr
import picked_group_fdr.utils.multiprocessing_pool as pool
import picked_group_fdr.digest as digest
import picked_group_fdr.digestion_params as digestion_params

logger = logging.getLogger()
logger.setLevel(logging.INFO)


def run_picked_group_fdr_all(evidence_files, pout_files, fasta_file, output_dir, digest_params, input_type):
    try:
        if len(output_dir) == 0:
            raise RuntimeError("Please specify an output folder")
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        
        if input_type == "rescoring":
            run_update_evidence_rescoring(evidence_files, pout_files, output_dir)
            run_picked_group_fdr(output_dir, fasta_file, digest_params)
        elif input_type == "percolator":
            run_merge_pout(pout_files, fasta_file, output_dir, digest_params)
            run_picked_group_fdr_percolator_input(output_dir, fasta_file)
        else:
            run_andromeda_to_pin(evidence_files, fasta_file, output_dir, digest_params)
            run_mokapot(output_dir)
            run_update_evidence(evidence_files, output_dir)
            run_picked_group_fdr(output_dir, fasta_file, digest_params)
        
    except SystemExit as e:
        logger.error(f"Error while running Picked Group FDR, exited with error code: {e}.")
    except Exception as e:
        logger.error(f"Error while running Picked Group FDR: {e}")


def run_andromeda_to_pin(evidence_files, fasta_file, output_dir, digest_params):
    andromeda2pin.main(evidence_files + ['--outputTab', f'{output_dir}/pin.tab', '--databases', fasta_file] + digest_params)
    

def run_mokapot(output_dir):
    np.random.seed(0) # TODO: Make seed configurable
    psms = mokapot.read_pin(f'{output_dir}/pin.tab')
    with parallel_backend("threading"):
        results, models = mokapot.brew(psms)
    results.to_txt(dest_dir=output_dir, decoys = True)
    

def run_update_evidence(evidence_files, output_dir):
    update_evidence.main(
        ['--mq_evidence'] + evidence_files + 
        ['--perc_results', f'{output_dir}/mokapot.psms.txt', f'{output_dir}/mokapot.decoys.psms.txt', 
         '--mq_evidence_out', f'{output_dir}/evidence_percolator.txt'])


def run_update_evidence_rescoring(evidence_files, pout_files, output_dir):
    update_evidence.main(
        ['--mq_evidence'] + evidence_files + 
        ['--perc_results'] + pout_files +
        ['--mq_evidence_out', f'{output_dir}/evidence_percolator.txt',
         '--pout_input_type', 'prosit'])


def run_picked_group_fdr(output_dir, fasta_file, digest_params):
    picked_group_fdr.main(
        ['--mq_evidence', f'{output_dir}/evidence_percolator.txt',
         '--methods', 'picked_protein_group_mq_input', 
         '--do_quant',
         '--protein_groups_out', f'{output_dir}/proteinGroups_percolator.txt',
         '--fasta', fasta_file] + digest_params)


def run_merge_pout(pout_files, fasta_file, output_dir, digest_params):
    merge_pout.main(
        ['--perc_results'] + pout_files +
        ['--perc_merged', f'{output_dir}/pout_merged.txt', 
         '--fasta', fasta_file] + digest_params)


def run_picked_group_fdr_percolator_input(output_dir, fasta_file):
    picked_group_fdr.main(
        ['--perc_evidence', f'{output_dir}/pout_merged.txt',
         '--methods', 'picked_protein_group_no_remap', 
         '--protein_groups_out', f'{output_dir}/proteinGroups_percolator.txt',
         '--fasta', fasta_file])


# https://stackoverflow.com/questions/28655198/best-way-to-display-logs-in-pyqt#60528393
class QTextEditLogger(logging.Handler, QObject):
    appendPlainText = pyqtSignal(str)
    
    def __init__(self, parent):
        super().__init__()        
        formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        formatter.converter = time.gmtime
        self.setFormatter(formatter)

        QObject.__init__(self)
        self.widget = QtWidgets.QPlainTextEdit(parent)
        self.widget.setReadOnly(True)
        self.appendPlainText.connect(self.widget.appendPlainText)

    def emit(self, record):
        self.appendPlainText.emit(self.format(record))


# https://stackoverflow.com/questions/53288877/python-multiprocessing-sending-child-process-logging-to-gui-running-in-parent
class LogEmitter(QObject):
    sigLog = pyqtSignal(str)


class LogHandler(logging.Handler):
    def __init__(self):
        super().__init__()
        formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        formatter.converter = time.gmtime
        self.setFormatter(formatter)

        self.emitter = LogEmitter()

    def emit(self, record):
        msg = self.format(record)
        self.emitter.sigLog.emit(msg)


class FileSelect(QtWidgets.QWidget):
    def __init__(self, file_type, file_extensions, file_hint='', folder_select=False, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.file_type = file_type
        self.file_hint = file_hint
        self.file_extensions = file_extensions
        
        self.label_text = f"Select {self.file_type} file"
        if folder_select:
            self.label_text = f"Select {self.file_type} folder"
        
        self.file_hint_text = ""
        if len(file_hint) > 0:
            self.file_hint_text = f'<br><font color="grey">{self.file_hint}</font>'
        self.label = QtWidgets.QLabel(self.label_text + self.file_hint_text)
        
        self.hbox_layout = QtWidgets.QHBoxLayout()
        self.hbox_layout.setContentsMargins(0, 0, 0, 0)
        
        self.line_edit = QtWidgets.QLineEdit()
        self.browse_button = QtWidgets.QPushButton("Browse")
        if folder_select:
            self.browse_button.clicked.connect(self.select_dir)
        else:
            self.browse_button.clicked.connect(self.select_file)
        self.hbox_layout.addWidget(self.line_edit, stretch = 1)
        self.hbox_layout.addWidget(self.browse_button)
        
        self.setLayout(self.hbox_layout)
    
    def select_file(self):
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(self, self.label_text, '', self.file_extensions)
        self.line_edit.setText(filename)
    
    def select_dir(self):
        output_dir = QtWidgets.QFileDialog.getExistingDirectory(self, self.label_text , '', QtWidgets.QFileDialog.ShowDirsOnly)
        self.line_edit.setText(output_dir)
    
    def get_file(self):
        return self.line_edit.text()
    
    def setButtonsEnabled(self, enable):
        self.browse_button.setEnabled(enable)
        

class MultiFileSelect(QtWidgets.QWidget):
    def __init__(self, file_type, file_extensions, file_hint='', *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.file_type = file_type
        self.file_hint = file_hint
        self.file_extensions = file_extensions
        
        self.label_text = f"Select {self.file_type} file(s)"
        
        self.file_hint_text = ""
        if len(file_hint) > 0:
            self.file_hint_text = f'<br><font color="grey">{self.file_hint}</font>'
        self.label = QtWidgets.QLabel(self.label_text + self.file_hint_text)
        
        self.hbox_layout = QtWidgets.QHBoxLayout()
        self.hbox_layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(self.hbox_layout)
        
        self.line_edit = QtWidgets.QListWidget()
        self.line_edit.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        
        self.vbox_layout = QtWidgets.QVBoxLayout()
        self.vbox_layout.setContentsMargins(0, 0, 0, 0)
        self.vbox_layout.setAlignment(Qt.AlignTop)
        
        self.browse_button = QtWidgets.QPushButton("Add")
        self.browse_button.clicked.connect(self.add_files)
        
        self.remove_button = QtWidgets.QPushButton("Remove")
        self.remove_button.clicked.connect(self.remove_files)
        
        self.vbox_layout.addWidget(self.browse_button)
        self.vbox_layout.addWidget(self.remove_button)
        
        self.hbox_layout.addWidget(self.line_edit, stretch = 1)
        self.hbox_layout.addLayout(self.vbox_layout)
    
    def add_files(self):
        filenames, _ = QtWidgets.QFileDialog.getOpenFileNames(self, self.label_text, '', self.file_extensions)
        self.line_edit.addItems(filenames)
    
    def remove_files(self):
        selected_items = self.line_edit.selectedItems()
        if not selected_items:
            return
        
        for item in selected_items:
            self.line_edit.takeItem(self.line_edit.row(item))
    
    def get_files(self):
        return [str(self.line_edit.item(i).text()) for i in range(self.line_edit.count())]
    
    def setButtonsEnabled(self, enable):
        self.browse_button.setEnabled(enable)
        self.remove_button.setEnabled(enable)
        

class DigestionParametersGroup(QtWidgets.QGroupBox):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
        self.digestion_group_layout = QtWidgets.QGridLayout()
        self.setLayout(self.digestion_group_layout)
        
        self.min_length_label = QtWidgets.QLabel("Min peptide length")
        self.min_length_spinbox = QtWidgets.QSpinBox()
        self.min_length_spinbox.setValue(digestion_params.MIN_PEPLEN_DEFAULT)
        self.min_length_spinbox.setRange(1,20)
        
        self.max_length_label = QtWidgets.QLabel("Max peptide length")
        self.max_length_spinbox = QtWidgets.QSpinBox()
        self.max_length_spinbox.setValue(digestion_params.MAX_PEPLEN_DEFAULT)
        self.max_length_spinbox.setRange(1,100)
        
        self.max_cleavages_label = QtWidgets.QLabel("Max miscleavages")
        self.max_cleavages_spinbox = QtWidgets.QSpinBox()
        self.max_cleavages_spinbox.setValue(digestion_params.CLEAVAGES_DEFAULT)
        self.max_cleavages_spinbox.setRange(1,10)
        
        self.enzyme_label = QtWidgets.QLabel("Enzyme")
        self.enzyme_select = QtWidgets.QComboBox()
        self.enzyme_select.addItems(digest.ENZYME_CLEAVAGE_RULES.keys())
        self.enzyme_select.setCurrentText(digestion_params.ENZYME_DEFAULT)
        
        self.digestion_label = QtWidgets.QLabel("Digestion")
        self.digestion_select = QtWidgets.QComboBox()
        self.digestion_select.addItems(["full", "semi", "none"])
        self.digestion_select.setCurrentText(digestion_params.DIGESTION_DEFAULT)
        
        self.special_aas_label = QtWidgets.QLabel("Special AAs")
        self.special_aas_line_edit = QtWidgets.QLineEdit()
        self.special_aas_line_edit.setText(digestion_params.SPECIAL_AAS_DEFAULT)
        
        self.digestion_group_layout.addWidget(self.enzyme_label, 0, 0)
        self.digestion_group_layout.addWidget(self.enzyme_select, 0, 1)
        self.digestion_group_layout.addWidget(self.min_length_label, 0, 2)
        self.digestion_group_layout.addWidget(self.min_length_spinbox, 0, 3)
        self.digestion_group_layout.addWidget(self.digestion_label, 0, 4)
        self.digestion_group_layout.addWidget(self.digestion_select, 0, 5)
        
        self.digestion_group_layout.addWidget(self.max_cleavages_label, 1, 0)
        self.digestion_group_layout.addWidget(self.max_cleavages_spinbox, 1, 1)
        self.digestion_group_layout.addWidget(self.max_length_label, 1, 2)
        self.digestion_group_layout.addWidget(self.max_length_spinbox, 1, 3)
        self.digestion_group_layout.addWidget(self.special_aas_label, 1, 4)
        self.digestion_group_layout.addWidget(self.special_aas_line_edit, 1, 5)
        
        for col in range(6):
            self.digestion_group_layout.setColumnStretch(col, 1)
    
    def get_params(self):
        return ["--min-length", str(self.min_length_spinbox.value()),
                "--max-length", str(self.max_length_spinbox.value()),
                "--cleavages", str(self.max_cleavages_spinbox.value()),
                "--enzyme", self.enzyme_select.currentText(),
                "--digestion", self.digestion_select.currentText(),
                "--special-aas", self.special_aas_line_edit.text()]

class MainWindow(QtWidgets.QWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
        self.setWindowTitle("Picked group FDR")
        
        self.tabs = QtWidgets.QTabWidget()
        self._add_mq_input_tab()
        self._add_percolator_input_tab()
        self._add_rescoring_input_tab()
        
        self.fasta_widget = FileSelect('fasta', 'Fasta file (*.fasta *.fa)')
        self.output_dir_widget = FileSelect('output', '', folder_select=True)        
        self.digestion_group = DigestionParametersGroup("Digestion parameters")
        
        layout = QtWidgets.QFormLayout()
        self.setLayout(layout)
        
        layout.addRow(self.tabs)
        layout.addRow(self.fasta_widget.label, self.fasta_widget)
        layout.addRow(self.output_dir_widget.label, self.output_dir_widget)
        layout.addRow(self.digestion_group)
    
        self._add_run_button(layout)
        self._add_log_textarea(layout)
        
        # sets up handler that will be used by QueueListener
        # which will update the LogDialog
        handler = LogHandler()
        handler.setLevel(logging.INFO) # TODO: this doesn't work, debug messages still make it through
        handler.emitter.sigLog.connect(self.log_text_area.widget.appendPlainText)
        
        self.q = multiprocessing.Queue()
        self.ql = QueueListener(self.q, handler)
        self.ql.start()

        self.pool = pool.JobPool(processes=1, warningFilter="default", queue=self.q)
        
        self.resize(800, self.height())

    def _add_mq_input_tab(self):
        self.mq_tab = QtWidgets.QWidget()
        
        self.mq_layout = QtWidgets.QFormLayout()
        self.evidence_widget = MultiFileSelect('evidence.txt', 'Tab separated file (*.txt)')
        self.mq_layout.addRow(self.evidence_widget.label, self.evidence_widget)
        
        self.mq_tab.setLayout(self.mq_layout)
        self.tabs.addTab(self.mq_tab, "MaxQuant input")
    
    def _add_percolator_input_tab(self):
        self.percolator_tab = QtWidgets.QWidget()
        
        self.percolator_layout = QtWidgets.QFormLayout()
        self.pout_widget = MultiFileSelect('percolator output', '', 'Add both target and decoy results!')
        self.percolator_layout.addRow(self.pout_widget.label, self.pout_widget)
        
        self.percolator_tab.setLayout(self.percolator_layout)
        self.tabs.addTab(self.percolator_tab, "Percolator input")
    
    def _add_rescoring_input_tab(self):
        self.rescoring_tab = QtWidgets.QWidget()
        
        self.rescoring_layout = QtWidgets.QFormLayout()
        self.evidence_widget_rescoring = MultiFileSelect('evidence.txt', 'Tab separated file (*.txt)')
        self.pout_widget_rescoring = MultiFileSelect('percolator output', '', 'Add both target and decoy results!')
        self.rescoring_layout.addRow(self.evidence_widget_rescoring.label, self.evidence_widget_rescoring)
        self.rescoring_layout.addRow(self.pout_widget_rescoring.label, self.pout_widget_rescoring)
        
        self.rescoring_tab.setLayout(self.rescoring_layout)
        self.tabs.addTab(self.rescoring_tab, "Rescoring input (e.g. Prosit)")
        
    def _add_run_button(self, layout):    
        self.run_button = QtWidgets.QPushButton("Run")
        self.run_button.clicked.connect(self.run_picked)
        self.run_button.setContentsMargins(20,100,20,100)
        
        layout.addRow(self.run_button)
    
    def _add_log_textarea(self, layout):    
        self.log_text_area = QTextEditLogger(self)
        self.log_text_area.setLevel(logging.INFO)
        self.log_text_area.widget.setMinimumHeight(self.log_text_area.widget.sizeHint().height())
        logger.addHandler(self.log_text_area)
        
        self.log_text_widget = QtWidgets.QWidget()
        
        self.log_text_hbox_layout = QtWidgets.QHBoxLayout()
        self.log_text_hbox_layout.setContentsMargins(0, 0, 0, 0)
        
        self.log_text_button_vbox_layout = QtWidgets.QVBoxLayout()
        self.log_text_button_vbox_layout.setContentsMargins(0, 0, 0, 0)
        self.log_text_button_vbox_layout.setAlignment(Qt.AlignTop)
        
        self.log_text_clear_button = QtWidgets.QPushButton("Clear")
        self.log_text_clear_button.clicked.connect(lambda: self.log_text_area.widget.clear())
        
        self.log_text_save_button = QtWidgets.QPushButton("Save")
        self.log_text_save_button.clicked.connect(self.save_log_file)
        
        self.log_text_button_vbox_layout.addWidget(self.log_text_clear_button)
        self.log_text_button_vbox_layout.addWidget(self.log_text_save_button)
        
        self.log_text_hbox_layout.addWidget(self.log_text_area.widget, stretch = 1)
        self.log_text_hbox_layout.addLayout(self.log_text_button_vbox_layout)
        
        self.log_text_widget.setLayout(self.log_text_hbox_layout)
             
        layout.addRow(self.log_text_widget)
    
    def get_output_dir(self):
        output_dir = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select output folder' , '', QtWidgets.QFileDialog.ShowDirsOnly)
        self.output_dir_line_edit.setText(output_dir)
    
    def save_log_file(self):
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Save log to file', 'picked_protein_group_fdr.log', 'Log file (*.txt *.log)' )
        if not filename:
            return
        
        logger.info(f"Saved log messages to {filename}")
        with open(filename, 'w') as f:
            f.write(self.log_text_area.widget.toPlainText())
        
    def set_buttons_enabled_state(self, enable):
        self.evidence_widget.setButtonsEnabled(enable)
        self.pout_widget.setButtonsEnabled(enable)
        self.evidence_widget_rescoring.setButtonsEnabled(enable)
        self.pout_widget_rescoring.setButtonsEnabled(enable)
        
        self.fasta_widget.setButtonsEnabled(enable)
        self.output_dir_widget.setButtonsEnabled(enable)

        # Cannot stop a QThread if it doesn't have an own event loop
        self.run_button.clicked.disconnect()
        if enable:
            self.run_button.setText("Run")
            self.run_button.clicked.connect(self.run_picked)
        else:
            self.run_button.setText("Stop")
            self.run_button.clicked.connect(self.stop_picked)
        
    def run_picked(self):
        fasta_file = self.fasta_widget.get_file()
        output_dir = self.output_dir_widget.get_file()
        digest_params = self.digestion_group.get_params()
        
        evidence_files, pout_files = list(), list()
        if self.tabs.currentIndex() == 2:
            input_type = "rescoring"
            evidence_files = self.evidence_widget_rescoring.get_files()
            pout_files = self.pout_widget_rescoring.get_files()
        elif self.tabs.currentIndex() == 1:
            input_type = "percolator"
            pout_files = self.pout_widget.get_files()
        else:
            input_type = "mq"
            evidence_files = self.evidence_widget.get_files()
        
        self.set_buttons_enabled_state(False)
        self.pool.applyAsync(run_picked_group_fdr_all, (evidence_files, pout_files, fasta_file, output_dir, digest_params, input_type), callback=self.on_picked_finished)
    
    def on_picked_finished(self, return_code):
        self.set_buttons_enabled_state(True)
        
    def stop_picked(self):
        self.pool.stopPool()
        self.on_picked_finished(-2)
        
        logger.info("Picked Group FDR stopped by user")
        
        self.pool = pool.JobPool(processes=1, warningFilter="default", queue=self.q)
    
    def closeEvent(self, _):
        self.stop_picked()


if __name__ == '__main__':
    if sys.platform.startswith('win'):
        # On Windows calling this function is necessary when combined with pyinstaller: 
        # https://stackoverflow.com/questions/24944558/pyinstaller-built-windows-exe-fails-with-multiprocessing
        multiprocessing.freeze_support()
    else:
        # On Linux calling this function is necessary when combined with pyqt: 
        # https://stackoverflow.com/questions/29556291/multiprocessing-with-qt-works-in-windows-but-not-linux
        multiprocessing.set_start_method('spawn', force=True)

    app = QtWidgets.QApplication(sys.argv)
    w = MainWindow()
    w.show()
    app.exec()
