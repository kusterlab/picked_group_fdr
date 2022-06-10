from PyQt5 import QtWidgets
from PyQt5.QtCore import QObject, QThread, pyqtSignal # https://realpython.com/python-pyqt-qthread/

import sys

# mock matplotlib imports
sys.modules['matplotlib'] = __import__('matplotlib_mock')

import os
import logging
from logging.handlers import QueueListener
import time
import datetime
import multiprocessing

import picked_group_fdr.pipeline.andromeda2pin as andromeda2pin
import picked_group_fdr.pipeline.update_evidence_from_pout as update_evidence
import picked_group_fdr.picked_group_fdr as picked_group_fdr
import picked_group_fdr.utils.multiprocessing_pool as pool

import mokapot
import numpy as np
from joblib import parallel_backend


logger = logging.getLogger()
logger.setLevel(logging.INFO)


def run_picked_group_fdr_all(evidence_file, fasta_file, output_dir, digest_params):
    try:
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        run_andromeda_to_pin(evidence_file, fasta_file, output_dir, digest_params)
        run_mokapot(output_dir)
        run_update_evidence(evidence_file, output_dir)
        run_picked_group_fdr(fasta_file, output_dir, digest_params)
    except SystemExit as e:
        logger.info(f"Error while running Picked Group FDR, exited with error code {e}.")
    except Exception as e:
        logger.info(f"Error while running Picked Group FDR: {e}")


def run_andromeda_to_pin(evidence_file, fasta_file, output_dir, digest_params):
    andromeda2pin.main([evidence_file, '--outputTab', f'{output_dir}/pin.tab', '--databases', fasta_file] + digest_params.split())
    

def run_mokapot(output_dir):
    np.random.seed(0) # TODO: Make seed configurable
    psms = mokapot.read_pin(f'{output_dir}/pin.tab')
    with parallel_backend("threading"):
        results, models = mokapot.brew(psms)
    results.to_txt(dest_dir=output_dir, decoys = True)
    

def run_update_evidence(evidence_file, output_dir):
    update_evidence.main(['--mq_evidence', evidence_file, '--perc_results', f'{output_dir}/mokapot.psms.txt', f'{output_dir}/mokapot.decoy.psms.txt', '--mq_evidence_out', f'{output_dir}/evidence_percolator.txt'])


def run_picked_group_fdr(fasta_file, output_dir, digest_params):
    picked_group_fdr.main(['--mq_evidence', f'{output_dir}/evidence_percolator.txt', '--protein_groups_out', f'{output_dir}/proteinGroups_percolator.txt', '--do_quant', '--fasta', fasta_file] + digest_params.split())


# https://stackoverflow.com/questions/28655198/best-way-to-display-logs-in-pyqt#60528393
class QTextEditLogger(logging.Handler, QObject):
    appendPlainText = pyqtSignal(str)
    
    def __init__(self, parent):
        super().__init__()
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
        self.emitter = LogEmitter()

    def emit(self, record):
        msg = self.format(record)
        self.emitter.sigLog.emit(msg)


class MainWindow(QtWidgets.QWidget):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
        layout = QtWidgets.QFormLayout()
        
        self.setWindowTitle("Picked group FDR")
        
        self._add_evidence_field(layout)
        self._add_fasta_field(layout)
        self._add_output_dir_field(layout)
        self._add_digestion_params_field(layout)
        self._add_run_button(layout)
        self._add_log_textarea(layout)
        
        self.setLayout(layout)
        
        # sets up handler that will be used by QueueListener
        # which will update the LogDialoag
        handler = LogHandler()
        handler.emitter.sigLog.connect(self.log_text_area.widget.appendPlainText)
        
        self.q = multiprocessing.Queue()
        self.ql = QueueListener(self.q, handler)
        self.ql.start()

        self.pool = pool.JobPool(processes=1, warningFilter="default", queue=self.q)
        
        self.resize(700, self.height())

    def _add_evidence_field(self, layout):
        # evidence.txt input
        self.evidence_label = QtWidgets.QLabel("Select an evidence.txt file")
        #self.evidence_label.setMargin(10)
        
        self.evidence_widget = QtWidgets.QWidget()
        
        self.evidence_hbox_layout = QtWidgets.QHBoxLayout()
        self.evidence_hbox_layout.setContentsMargins(0, 0, 0, 0)
        
        self.evidence_line_edit = QtWidgets.QLineEdit()
        self.evidence_browse_button = QtWidgets.QPushButton("Browse")
        self.evidence_browse_button.clicked.connect(self.get_evidence_file)
        self.evidence_hbox_layout.addWidget(self.evidence_line_edit, stretch = 1)
        self.evidence_hbox_layout.addWidget(self.evidence_browse_button)
        
        self.evidence_widget.setLayout(self.evidence_hbox_layout)
        
        layout.addRow(self.evidence_label, self.evidence_widget)

    def _add_fasta_field(self, layout):
        # fasta file input
        self.fasta_label = QtWidgets.QLabel("Select a fasta file")
        #self.fasta_label.setMargin(10)
        
        self.fasta_widget = QtWidgets.QWidget()
        
        self.fasta_hbox_layout = QtWidgets.QHBoxLayout()
        self.fasta_hbox_layout.setContentsMargins(0, 0, 0, 0)
        
        self.fasta_line_edit = QtWidgets.QLineEdit()
        self.fasta_browse_button = QtWidgets.QPushButton("Browse")
        self.fasta_browse_button.clicked.connect(self.get_fasta_file)
        self.fasta_hbox_layout.addWidget(self.fasta_line_edit, stretch = 1)
        self.fasta_hbox_layout.addWidget(self.fasta_browse_button)
        
        self.fasta_widget.setLayout(self.fasta_hbox_layout)
        
        layout.addRow(self.fasta_label, self.fasta_widget)
    
    def _add_output_dir_field(self, layout):
        # fasta file input
        self.output_dir_label = QtWidgets.QLabel("Select output folder")
        #self.output_dir_label.setMargin(10)
        
        self.output_dir_widget = QtWidgets.QWidget()
        
        self.output_dir_hbox_layout = QtWidgets.QHBoxLayout()
        self.output_dir_hbox_layout.setContentsMargins(0, 0, 0, 0)
        
        self.output_dir_line_edit = QtWidgets.QLineEdit()
        self.output_dir_browse_button = QtWidgets.QPushButton("Browse")
        self.output_dir_browse_button.clicked.connect(self.get_output_dir)
        self.output_dir_hbox_layout.addWidget(self.output_dir_line_edit, stretch = 1)
        self.output_dir_hbox_layout.addWidget(self.output_dir_browse_button)
        
        self.output_dir_widget.setLayout(self.output_dir_hbox_layout)
        
        layout.addRow(self.output_dir_label, self.output_dir_widget)

    def _add_digestion_params_field(self, layout):
        # additional parameters input, TODO: make user friendly options for each digestion parameter
        self.args_label = QtWidgets.QLabel("Digestion parameters")
        #self.args_label.setMargin(10)
        
        self.args_line_edit = QtWidgets.QLineEdit()
        
        layout.addRow(self.args_label, self.args_line_edit)
    
    def _add_run_button(self, layout):    
        self.run_button = QtWidgets.QPushButton("Run")
        self.run_button.clicked.connect(self.run_picked)
        self.run_button.setContentsMargins(20,100,20,100)
        
        layout.addRow(self.run_button)
    
    def _add_log_textarea(self, layout):    
        self.log_text_area = QTextEditLogger(self)
        self.log_text_area.setLevel(logging.INFO)
        logger.addHandler(self.log_text_area)
             
        layout.addRow(self.log_text_area.widget)
        
    def get_evidence_file(self):
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Open evidence file' , '','Tab separated file (*.txt)' )
        self.evidence_line_edit.setText(filename)
    
    def get_fasta_file(self):
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Open fasta file' , '','Fasta file (*.fasta *.fa)' )
        self.fasta_line_edit.setText(filename)
    
    def get_output_dir(self):
        output_dir = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select output folder' , '', QtWidgets.QFileDialog.ShowDirsOnly)
        self.output_dir_line_edit.setText(output_dir)
    
    def set_buttons_enabled_state(self, enable):
        self.evidence_browse_button.setEnabled(enable)
        self.fasta_browse_button.setEnabled(enable)
        self.output_dir_browse_button.setEnabled(enable)
        #self.run_button.setEnabled(enable)
        # Cannot stop a QThread if it doesn't have an own event loop
        self.run_button.clicked.disconnect()
        if enable:
            self.run_button.setText("Run")
            self.run_button.clicked.connect(self.run_picked)
        else:
            self.run_button.setText("Stop")
            self.run_button.clicked.connect(self.stop_picked)
        
    def run_picked(self):
        evidence_file = self.evidence_line_edit.text()
        fasta_file = self.fasta_line_edit.text()
        output_dir = self.output_dir_line_edit.text()
        digest_params = self.args_line_edit.text()
        
        self.set_buttons_enabled_state(False)
        self.pool.applyAsync(run_picked_group_fdr_all, (evidence_file, fasta_file, output_dir, digest_params), callback=self.on_picked_finished)
    
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
        # On Windows calling this function is necessary when combined with pyinstaller: https://stackoverflow.com/questions/24944558/pyinstaller-built-windows-exe-fails-with-multiprocessing
        multiprocessing.freeze_support()
    else:
        # On Linux calling this function is necessary when combined with pyqt: https://stackoverflow.com/questions/29556291/multiprocessing-with-qt-works-in-windows-but-not-linux
        multiprocessing.set_start_method('spawn')
    app = QtWidgets.QApplication(sys.argv)
    w = MainWindow()
    w.show()
    app.exec()
