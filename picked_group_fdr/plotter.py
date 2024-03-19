import logging

import numpy as np

from . import fdr
from . import methods


logger = logging.getLogger(__name__)


class NoPlotter:
    def __init__(self):
        pass
    
    def init_plots(self):
        pass
    
    def decorate_plots(self):
        pass
        
    def save_plots(self):
        pass
    
    def show(self):
        pass
    
    def set_series_label_base(self, figure_label):
        pass
    
    def set_series_label(self, method_config: methods.MethodConfig, rescue_step: bool):
        pass
    
    def plot_qval_calibration_and_performance(self, reported_qvals, observed_qvals, absent_ratio = 1.0):
        pass


class Plotter:
    max_prot: float
    max_plotted_qval: float
    figure_base_fn: str
    plot_figures: bool
    label: str
    label_base: str
    
    import matplotlib.pyplot as plt
    
    def __init__(self, figure_base_fn, plot_figures):
        self.max_plotted_qval = 0.1 # 0.02
        self.max_prot = 0
        self.figure_base_fn = figure_base_fn
        self.plot_figures = plot_figures
        self.label = ""
    
    def _updateMaxProt(self, reported_qvals):
        self.max_prot = max([self.max_prot, fdr.count_below_threshold(reported_qvals, self.max_plotted_qval)])
        
    def init_plots(self):
        self.plt.figure(1, figsize=(8,6))
        self.plt.figure(2, figsize=(8,6))
        self.plt.figure(3, figsize=(8,6))
    
    def decorate_plots(self):
        increments = 200 if self.max_prot < 2000 else 1000
        max_prot = int(self.max_prot / increments + 1) * increments
        
        # calibration
        self.plt.figure(1)
        #self.plt.subplot(2,2,1, label = 'calibration')
        self.plt.plot([0, 1], [0, 1], 'k-', alpha = 0.5, linewidth = 1)
        self.plt.plot([0, 1.5], [0, 1], 'k--', alpha = 0.5, linewidth = 1)
        self.plt.plot([0, 0.67], [0, 1], 'k--', alpha = 0.5, linewidth = 1)
        self.plt.xlabel("Decoy q-value", fontsize = 14)
        self.plt.ylabel("Entrapment q-value", fontsize = 14)
        
        self.plt.xlim([0,self.max_plotted_qval])
        self.plt.ylim([0,self.max_plotted_qval])
        legend = self.plt.legend(fontsize = 14)
        self._right_align_legend(legend)
        
        for tick in self.plt.gca().xaxis.get_major_ticks():
            tick.label.set_fontsize(14) 
        for tick in self.plt.gca().yaxis.get_major_ticks():
            tick.label.set_fontsize(14) 
            
        self.plt.tight_layout()
        
        # performance entrapment
        self.plt.figure(2)
        #self.plt.subplot(2,2,2, label = 'performance entrapment')
        self.plt.plot([0.01, 0.01], [0, max_prot], 'k--', alpha = 0.5, linewidth = 1)
        self.plt.xlabel("Entrapment q-value", fontsize = 14)
        self.plt.ylabel("# proteins / protein groups", fontsize = 14)
        
        self.plt.xlim([0, self.max_plotted_qval])
        self.plt.ylim([0, max_prot])
        legend = self.plt.legend(fontsize = 14)
        self._right_align_legend(legend)
        
        for tick in self.plt.gca().xaxis.get_major_ticks():
            tick.label.set_fontsize(14) 
        for tick in self.plt.gca().yaxis.get_major_ticks():
            tick.label.set_fontsize(14) 
        
        self.plt.tight_layout()
        
        # performance decoy
        self.plt.figure(3)
        #self.plt.subplot(2,2,3, label = 'performance')
        self.plt.plot([0.01, 0.01], [0, max_prot], 'k--', alpha = 0.5, linewidth = 1)
        self.plt.xlabel("Decoy q-value", fontsize = 14)
        self.plt.ylabel("# proteins / protein groups", fontsize = 14)
        
        self.plt.xlim([0, self.max_plotted_qval])
        self.plt.ylim([0, max_prot])
        legend = self.plt.legend(fontsize = 14)
        self._right_align_legend(legend)
        
        for tick in self.plt.gca().xaxis.get_major_ticks():
            tick.label.set_fontsize(14) 
        for tick in self.plt.gca().yaxis.get_major_ticks():
            tick.label.set_fontsize(14) 
        
        self.plt.tight_layout()

    def save_plots(self):
        if self.figure_base_fn:
            logger.info(f"Saving calibration plot: {self.figure_base_fn}-calibration.png")
            self.plt.figure(1)
            self.plt.savefig(self.figure_base_fn + "-calibration.pdf")
            self.plt.savefig(self.figure_base_fn + "-calibration.png")
            
            logger.info(f"Saving entrapment performance plot: {self.figure_base_fn}-entrapment-performance.png")
            self.plt.figure(2)
            self.plt.savefig(self.figure_base_fn + "-entrapment-performance.pdf")
            self.plt.savefig(self.figure_base_fn + "-entrapment-performance.png")
            
            logger.info(f"Saving performance plot: {self.figure_base_fn}-performance.png")
            self.plt.figure(3)
            self.plt.savefig(self.figure_base_fn + "-performance.pdf")
            self.plt.savefig(self.figure_base_fn + "-performance.png")
    
    def show(self):
        if self.plot_figures:
            self.plt.show()
    
    def plot_qval_calibration_and_performance(self, reported_qvals, observed_qvals, absent_ratio = 1.0):
        self.plt.figure(1)
        #self.plt.subplot(2,2,1, label = 'calibration')
        self.plt.plot(absent_ratio*np.array(reported_qvals), np.array(observed_qvals), '-', label = self.label)
        
        self.plt.figure(2)
        #self.plt.subplot(2,2,2, label = 'performance entrapment')
        self.plt.plot(observed_qvals, range(len(reported_qvals)), '-', label = self.label)
        
        self.plt.figure(3)
        #self.plt.subplot(2,2,3, label = 'performance')
        self.plt.plot(reported_qvals, range(len(reported_qvals)), '-', label = self.label)
        
        self._updateMaxProt(reported_qvals)

    def _right_align_legend(self, legend):    
        # get the width of your widest label, since every label will need 
        # to shift by this amount after we align to the right
        renderer = self.plt.gcf().canvas.get_renderer()
        shift = max([t.get_window_extent(renderer).width for t in legend.get_texts()])
        for t in legend.get_texts():
                t.set_ha('right') # ha is alias for horizontalalignment
                t.set_position((shift,0))
    
    def set_series_label_base(self, figure_label):
        self.label_base = figure_label
    
    def set_series_label(self, method_config: methods.MethodConfig, rescue_step: bool):
        label = self.label_base + " " if self.label_base else ""
        short_desc = method_config.short_description(rescue_step, sep=",")
        label += f'({short_desc})'
        self.label = label


class PlotterFactory:
    @staticmethod
    def get_plotter(figure_base_fn, plot_figures):
        plotter = NoPlotter()
        if figure_base_fn or plot_figures:
            plotter = Plotter(figure_base_fn, plot_figures)            
        plotter.init_plots()
        return plotter


