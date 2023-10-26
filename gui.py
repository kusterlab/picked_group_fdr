from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt, QObject, pyqtSignal # https://realpython.com/python-pyqt-qthread/

import sys

# mock matplotlib imports
sys.modules['matplotlib'] = __import__('matplotlib_mock')

import logging
from logging.handlers import QueueListener
import time
import multiprocessing

import picked_group_fdr.utils.multiprocessing_pool as pool
import picked_group_fdr.digest as digest
import picked_group_fdr.digestion_params as digestion_params
import picked_group_fdr.pipeline as pipeline

logger = logging.getLogger()
logger.setLevel(logging.INFO)


# https://stackoverflow.com/questions/28655198/best-way-to-display-logs-in-pyqt#60528393
class QTextEditLogger(logging.Handler, QObject):
    """Text area which prints logs from both the pipeline and GUI.

    Messages from the pipeline are passed through with the pyqtSignal.
    """    
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


class LogHandler(logging.Handler):
    def __init__(self):
        super().__init__()
        formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        formatter.converter = time.gmtime
        self.setFormatter(formatter)

        self.emitter = self.LogEmitter()

    def emit(self, record):
        msg = self.format(record)
        self.emitter.sigLog.emit(msg)
    
    # https://stackoverflow.com/questions/53288877/python-multiprocessing-sending-child-process-logging-to-gui-running-in-parent
    class LogEmitter(QObject):
        sigLog = pyqtSignal(str)



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
        

class MultiFileSelectTable(QtWidgets.QWidget):
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
        
        table_headers = {
            "Evidence file": 310, 
            "Min len": 70, 
            "Max len": 70,
            "Miscleav.": 70,
            "Enzyme": 100, 
            "Digestion": 100, 
            "Special AAs": 90
        }
        self.table_edit = QtWidgets.QTableWidget(0, len(table_headers))
        self.table_edit.setMinimumHeight(150);
        self.table_edit.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.table_edit.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.table_edit.setTextElideMode(Qt.ElideMiddle)
        self.table_edit.setWordWrap(False)
        for column_idx, column_width in enumerate(table_headers.values()):
            self.table_edit.setColumnWidth(column_idx, column_width)
        self.table_edit.setHorizontalHeaderLabels(table_headers.keys())
        self.table_edit.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.Stretch)
        self.table_edit.horizontalHeader().setHighlightSections(False)
        self.table_edit.verticalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeToContents)
        self.table_edit.verticalHeader().setMaximumSectionSize(22)
        
        self.add_remove_buttons = QtWidgets.QHBoxLayout()
        self.add_remove_buttons.setContentsMargins(0, 0, 0, 0)
        self.add_remove_buttons.setAlignment(Qt.AlignTop)
        
        self.browse_button = QtWidgets.QPushButton("Add")
        self.browse_button.setSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        self.browse_button.clicked.connect(self.add_files)
        
        self.remove_button = QtWidgets.QPushButton("Remove")
        self.remove_button.setSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        self.remove_button.clicked.connect(self.remove_files)
        
        self.add_remove_buttons.addWidget(self.browse_button)
        self.add_remove_buttons.addWidget(self.remove_button)
        
        self.hbox_layout.addWidget(self.table_edit, stretch = 1)
    
    def add_files(self):
        filenames, _ = QtWidgets.QFileDialog.getOpenFileNames(self, self.label_text, '', self.file_extensions)
        for filename in filenames:
            self.table_edit.insertRow(self.table_edit.rowCount())
            row_idx = self.table_edit.rowCount()-1
            self.table_edit.setItem(row_idx, 0, QtWidgets.QTableWidgetItem(filename))
            
            min_length_spinbox = QtWidgets.QSpinBox()
            min_length_spinbox.setValue(digestion_params.MIN_PEPLEN_DEFAULT)
            min_length_spinbox.setRange(1,20)
            self.table_edit.setCellWidget(row_idx, 1, min_length_spinbox)

            max_length_spinbox = QtWidgets.QSpinBox()
            max_length_spinbox.setValue(digestion_params.MAX_PEPLEN_DEFAULT)
            max_length_spinbox.setRange(1,100)
            self.table_edit.setCellWidget(row_idx, 2, max_length_spinbox)

            max_cleavages_spinbox = QtWidgets.QSpinBox()
            max_cleavages_spinbox.setValue(digestion_params.CLEAVAGES_DEFAULT)
            max_cleavages_spinbox.setRange(1,10)
            self.table_edit.setCellWidget(row_idx, 3, max_cleavages_spinbox)
            
            enzyme_select = QtWidgets.QComboBox()
            enzyme_select.addItems(digest.ENZYME_CLEAVAGE_RULES.keys())
            enzyme_select.setCurrentText(digestion_params.ENZYME_DEFAULT)
            self.table_edit.setCellWidget(row_idx, 4, enzyme_select)

            digestion_select = QtWidgets.QComboBox()
            digestion_select.addItems(["full", "semi", "none"])
            digestion_select.setCurrentText(digestion_params.DIGESTION_DEFAULT)
            self.table_edit.setCellWidget(row_idx, 5, digestion_select)
            
            special_aas_line_edit = QtWidgets.QLineEdit()
            special_aas_line_edit.setText(digestion_params.SPECIAL_AAS_DEFAULT)
            self.table_edit.setCellWidget(row_idx, 6, special_aas_line_edit)
    
    def remove_files(self):
        selected_items = self.table_edit.selectedItems()
        if not selected_items:
            return
        
        for item in selected_items:
            self.table_edit.removeRow(self.table_edit.row(item))
    
    def get_files(self):
        return [str(self.table_edit.item(i, 0).text()) for i in range(self.table_edit.rowCount())]
    
    def setButtonsEnabled(self, enable):
        self.browse_button.setEnabled(enable)
        self.remove_button.setEnabled(enable)
    
    def get_digestion_params(self):
        digestion_params_list = list()
        for row_idx in range(self.table_edit.rowCount()):
            digestion_params_list.append(
                digestion_params.DigestionParams(
                    self.table_edit.cellWidget(row_idx, 4).currentText(), # enzyme
                    self.table_edit.cellWidget(row_idx, 5).currentText(), # digestion
                    self.table_edit.cellWidget(row_idx, 1).value(), # min_length
                    self.table_edit.cellWidget(row_idx, 2).value(), # max_length
                    self.table_edit.cellWidget(row_idx, 3).value(), # max_cleavages
                    self.table_edit.cellWidget(row_idx, 6).text(), # special_aas
                    fasta_contains_decoys=False,
                )
            )
        return digestion_params_list

        # return digestion_params.DigestionParams(
        #     self.enzyme_select.currentText(),
        #     self.digestion_select.currentText(),
        #     self.min_length_spinbox.value(),
        #     self.max_length_spinbox.value(),
        #     self.max_cleavages_spinbox.value(),
        #     self.special_aas_line_edit.text(),
        #     False
        # )


class MainWindow(QtWidgets.QWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
        self.setWindowTitle("Picked group FDR")
        
        self.tabs = QtWidgets.QTabWidget()
        self._add_mq_input_tab()
        self._add_percolator_input_tab()
        
        self.output_dir_widget = FileSelect('output', '', folder_select=True)
        
        layout = QtWidgets.QFormLayout()
        self.setLayout(layout)
        
        layout.addRow(self.tabs)
        layout.addRow(self.output_dir_widget.label, self.output_dir_widget)
    
        self._add_run_button(layout)
        self._add_log_textarea(layout)
        
        # sets up handler that will be used by QueueListener
        # which will update the LogDialog with messages from the pipeline
        handler = LogHandler()
        handler.setLevel(logging.INFO) # TODO: this doesn't work, debug messages still make it through
        handler.emitter.sigLog.connect(self.log_text_area.widget.appendPlainText)
        
        self.q = multiprocessing.Queue()
        self.ql = QueueListener(self.q, handler)
        self.ql.start()

        self.pool = pool.JobPool(processes=1, warningFilter="default", queue=self.q)
        
        self.resize(900, self.height())

    def _add_mq_input_tab(self):
        self.mq_tab = QtWidgets.QWidget()
        
        self.mq_layout = QtWidgets.QFormLayout()
        self.evidence_widget = MultiFileSelectTable('evidence.txt', 'Tab separated file (*.txt)')
        self.mq_layout.addRow(self.evidence_widget.label)
        self.mq_layout.addRow(self.evidence_widget)
        self.mq_layout.addRow(self.evidence_widget.add_remove_buttons)

        self.fasta_widget = MultiFileSelect('fasta', 'Fasta file (*.fasta *.fa)')
        self.mq_layout.addRow(self.fasta_widget.label, self.fasta_widget)

        self.rescoring_checkbox = QtWidgets.QCheckBox("Use rescoring results (e.g. Prosit)")
        self.rescoring_checkbox.stateChanged.connect(self._disable_enable_rescoring_panel)

        self.mq_layout.addRow(self.rescoring_checkbox)

        self.pout_widget_rescoring = MultiFileSelect('percolator output', '', 'Add both target and decoy results!')
        self.pout_widget_rescoring.setDisabled(True)
        self.pout_widget_rescoring.label.setDisabled(True)
        self.mq_layout.addRow(self.pout_widget_rescoring.label, self.pout_widget_rescoring)

        self.do_quant_checkbox = QtWidgets.QCheckBox("Do quantification")
        self.do_quant_checkbox.setChecked(True)
        self.do_quant_checkbox.stateChanged.connect(self._disable_enable_quant_panel)
        self.mq_layout.addRow(self.do_quant_checkbox)

        self.min_lfq_peptides_label = QtWidgets.QLabel("LFQ min peptides")
        self.min_lfq_peptides_spinbox = QtWidgets.QSpinBox()
        self.min_lfq_peptides_spinbox.setValue(2)
        self.min_lfq_peptides_spinbox.setRange(1,5)
        self.mq_layout.addRow(self.min_lfq_peptides_label, self.min_lfq_peptides_spinbox)

        self.mq_tab.setLayout(self.mq_layout)
        self.tabs.addTab(self.mq_tab, "MaxQuant/Prosit input")
    
    def _disable_enable_rescoring_panel(self):
        self.pout_widget_rescoring.setDisabled(not self.rescoring_checkbox.isChecked())
        self.pout_widget_rescoring.label.setDisabled(not self.rescoring_checkbox.isChecked())

    def _disable_enable_quant_panel(self):
        self.min_lfq_peptides_spinbox.setDisabled(not self.do_quant_checkbox.isChecked())
        self.min_lfq_peptides_label.setDisabled(not self.do_quant_checkbox.isChecked())

    def _add_percolator_input_tab(self):
        self.percolator_tab = QtWidgets.QWidget()
        
        self.percolator_layout = QtWidgets.QFormLayout()
        self.pout_widget = MultiFileSelect('percolator output', '', 'Add both target and decoy results!')
        self.percolator_layout.addRow(self.pout_widget.label, self.pout_widget)
        
        self.percolator_tab.setLayout(self.percolator_layout)
        self.tabs.addTab(self.percolator_tab, "Percolator/MSFragger input")
        
    def _add_run_button(self, layout):    
        self.run_button = QtWidgets.QPushButton("Run")
        self.run_button.setStyleSheet("background-color: CornflowerBlue")
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
        fasta_files = self.fasta_widget.get_files()
        output_dir = self.output_dir_widget.get_file()
        
        evidence_files, pout_files, digest_params = list(), list(), list()
        if self.tabs.currentIndex() == 1:
            input_type = "percolator"
            pout_files = self.pout_widget.get_files()
        else:
            input_type = "mq"
            if self.rescoring_checkbox.isChecked():
                input_type = "rescoring"
                pout_files = self.pout_widget_rescoring.get_files()
            
            evidence_files = self.evidence_widget.get_files()
            digest_params = self.evidence_widget.get_digestion_params()
        
        self.set_buttons_enabled_state(False)
        self.pool.applyAsync(pipeline.run_picked_group_fdr_all, (evidence_files, pout_files, fasta_files, output_dir, digest_params, input_type), callback=self.on_picked_finished)
    
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
