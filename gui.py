from PyQt5 import QtWidgets

import sys
import logging
import time
import datetime

import picked_group_fdr.pipeline.andromeda2pin as andromeda2pin
import picked_group_fdr.pipeline.update_evidence_from_pout as update_evidence
import picked_group_fdr.picked_group_fdr as picked_group_fdr

import mokapot
import numpy as np


logging.basicConfig(level=logging.INFO)


class GuiLogger(logging.Handler):
    def emit(self, record):
        self.edit.append(self.format(record))  # implementation of append_line omitted
        QtWidgets.QApplication.processEvents()


class MainWindow(QtWidgets.QWidget):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
        layout = QtWidgets.QFormLayout()
        
        self.setWindowTitle("Picked group FDR")
        
        self._add_evidence_field(layout)
        self._add_fasta_field(layout)
        self._add_digestion_params_field(layout)
        self._add_run_button(layout)
        self._add_log_textarea(layout)
        
        self.setLayout(layout)
        
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

    def _add_digestion_params_field(self, layout):
        # additional parameters input, TODO: make user friendly options for each digestion parameter
        self.args_label = QtWidgets.QLabel("Digestion parameters")
        #self.args_label.setMargin(10)
        
        self.args_line_edit = QtWidgets.QLineEdit()
        
        layout.addRow(self.args_label, self.args_line_edit)
    
    def _add_run_button(self, layout):    
        self.run_button = QtWidgets.QPushButton("Run")
        self.run_button.clicked.connect(self.process)
        self.run_button.setContentsMargins(20,100,20,100)
        
        layout.addRow(self.run_button)
    
    def _add_log_textarea(self, layout):    
        self.log_text_edit = QtWidgets.QTextEdit()
        
        self.logger = GuiLogger()
        self.logger.edit = self.log_text_edit
        logging.getLogger().addHandler(self.logger)
        
        layout.addRow(self.log_text_edit)
        
    def get_evidence_file(self):
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Open evidence file' , '','Tab separated file (*.txt)' )
        self.evidence_line_edit.setText(filename)
    
    def get_fasta_file(self):
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Open fasta file' , '','Fasta file (*.fasta *.fa)' )
        self.fasta_line_edit.setText(filename)
    
    def process(self):
        start = time.time()
        
        self.log_text_edit.append('Running...')
        
        self.log_text_edit.append('Converting evidence.txt to percolator input')
        QtWidgets.QApplication.processEvents()
        self.run_andromeda_to_pin()
        
        self.log_text_edit.append('Running mokapot')
        QtWidgets.QApplication.processEvents()
        self.run_mokapot()
        
        self.log_text_edit.append('Updating evidence.txt file')
        QtWidgets.QApplication.processEvents()
        self.run_update_evidence()
        
        self.log_text_edit.append('Calculating picked group FDRs')
        QtWidgets.QApplication.processEvents()
        self.run_picked_group_fdr()
        
        total_time = round(time.time() - start)
        total_time = str(datetime.timedelta(seconds=total_time))
        
        self.log_text_edit.append('=== Done ===')
        QtWidgets.QApplication.processEvents()
        
        logging.info("Picked group FDR analysis completed in %s", total_time)
        
    def run_andromeda_to_pin(self):
        evidence_file = self.evidence_line_edit.text()
        fasta_file = self.fasta_line_edit.text()
        digest_params = self.args_line_edit.text()
        andromeda2pin.main([evidence_file, '--outputTab', 'pin.tab', '--databases', fasta_file] + digest_params.split())
        
    def run_mokapot(self):
        np.random.seed(0) # TODO: Make seed configurable
        psms = mokapot.read_pin('pin.tab')
        results, models = mokapot.brew(psms)
        results.to_txt(decoys = True)
        
    def run_update_evidence(self):
        evidence_file = self.evidence_line_edit.text()
        update_evidence.main(['--mq_evidence', evidence_file, '--perc_results', 'mokapot.psms.txt', 'mokapot.decoys.psms.txt', '--mq_evidence_out', 'evidence_percolator.txt'])
    
    def run_picked_group_fdr(self):
        fasta_file = self.fasta_line_edit.text()
        digest_params = self.args_line_edit.text()
        picked_group_fdr.main(['--mq_evidence', 'evidence_percolator.txt', '--protein_groups_out', 'proteinGroups_percolator.txt', '--do_quant', '--fasta', fasta_file] + digest_params.split())
        

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    w = MainWindow()
    w.show()
    app.exec()
