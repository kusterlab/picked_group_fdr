import logging

import numpy as np

from . import fdr
from . import methods


logger = logging.getLogger(__name__)


class NoPlotter:
  def __init__(self):
    pass
  
  def initPlots(self):
    pass
  
  def decoratePlots(self):
    pass
    
  def savePlots(self):
    pass
  
  def show(self):
    pass
  
  def set_series_label_base(self, figureLabel):
    pass
  
  def set_series_label(self, scoreType, groupingStrategy, pickedStrategy, rescue_step):
    pass
  
  def plotQvalCalibrationAndPerformance(self, reportedQvals, observedQvals, label, absentRatio = 1.0):
    pass


class Plotter:
  maxProt: float
  maxPlottedQval: float
  figure_base_fn: str
  plot_figures: bool
  label: str
  label_base: str
  
  import matplotlib.pyplot as plt
  
  def __init__(self, figure_base_fn, plot_figures):
    self.maxPlottedQval = 0.02
    self.maxProt = 0
    self.figure_base_fn = figure_base_fn
    self.plot_figures = plot_figures
    self.label = ""
  
  def _updateMaxProt(self, reportedQvals):
    self.maxProt = max([self.maxProt, fdr.countBelowThreshold(reportedQvals, self.maxPlottedQval)])
    
  def initPlots(self):
    self.plt.figure(1, figsize=(8,6))
    self.plt.figure(2, figsize=(8,6))
    self.plt.figure(3, figsize=(8,6))
  
  def decoratePlots(self):
    increments = 200 if self.maxProt < 2000 else 1000
    maxProt = int(self.maxProt / increments + 1) * increments
    
    # calibration
    self.plt.figure(1)
    #self.plt.subplot(2,2,1, label = 'calibration')
    self.plt.plot([0, 1], [0, 1], 'k-', alpha = 0.5, linewidth = 1)
    self.plt.plot([0, 1.5], [0, 1], 'k--', alpha = 0.5, linewidth = 1)
    self.plt.plot([0, 0.67], [0, 1], 'k--', alpha = 0.5, linewidth = 1)
    self.plt.xlabel("Decoy q-value", fontsize = 14)
    self.plt.ylabel("Entrapment q-value", fontsize = 14)
    
    self.plt.xlim([0,self.maxPlottedQval])
    self.plt.ylim([0,self.maxPlottedQval])
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
    self.plt.plot([0.01, 0.01], [0, maxProt], 'k--', alpha = 0.5, linewidth = 1)
    self.plt.xlabel("Entrapment q-value", fontsize = 14)
    self.plt.ylabel("# proteins / protein groups", fontsize = 14)
    
    self.plt.xlim([0, self.maxPlottedQval])
    self.plt.ylim([0, maxProt])
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
    self.plt.plot([0.01, 0.01], [0, maxProt], 'k--', alpha = 0.5, linewidth = 1)
    self.plt.xlabel("Decoy q-value", fontsize = 14)
    self.plt.ylabel("# proteins / protein groups", fontsize = 14)
    
    self.plt.xlim([0, self.maxPlottedQval])
    self.plt.ylim([0, maxProt])
    legend = self.plt.legend(fontsize = 14)
    self._right_align_legend(legend)
    
    for tick in self.plt.gca().xaxis.get_major_ticks():
      tick.label.set_fontsize(14) 
    for tick in self.plt.gca().yaxis.get_major_ticks():
      tick.label.set_fontsize(14) 
    
    self.plt.tight_layout()

  def savePlots(self):
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
  
  def plotQvalCalibrationAndPerformance(self, reportedQvals, observedQvals, absentRatio = 1.0):
    self.plt.figure(1)
    #self.plt.subplot(2,2,1, label = 'calibration')
    self.plt.plot(absentRatio*np.array(reportedQvals), np.array(observedQvals), '-', label = self.label)
    
    #self.plt.hist([decoyScores, entrapmentScores, poolScores], bins = np.arange(0, 400, 10), label = ['decoy', 'entrapment', 'pool'])
    
    self.plt.figure(2)
    #self.plt.subplot(2,2,2, label = 'performance entrapment')
    self.plt.plot(observedQvals, range(len(reportedQvals)), '-', label = self.label)
    
    self.plt.figure(3)
    #self.plt.subplot(2,2,3, label = 'performance')
    self.plt.plot(reportedQvals, range(len(reportedQvals)), '-', label = self.label)
    
    self._updateMaxProt(reportedQvals)

  def _right_align_legend(self, legend):  
    # get the width of your widest label, since every label will need 
    # to shift by this amount after we align to the right
    renderer = self.plt.gcf().canvas.get_renderer()
    shift = max([t.get_window_extent(renderer).width for t in legend.get_texts()])
    for t in legend.get_texts():
        t.set_ha('right') # ha is alias for horizontalalignment
        t.set_position((shift,0))
  
  def set_series_label_base(self, figureLabel):
    self.label_base = figureLabel
  
  def set_series_label(self, scoreType, groupingStrategy, pickedStrategy, rescue_step):
    label = self.label_base + " " if self.label_base else ""
    short_desc = methods.short_description(scoreType, groupingStrategy, pickedStrategy, rescue_step, sep=",")
    label += f'({short_desc})'
    self.label = label


class PlotterFactory:
  @staticmethod
  def get_plotter(figure_base_fn, plot_figures):
    if figure_base_fn or plot_figures:
      return Plotter(figure_base_fn, plot_figures)
    else:
      return NoPlotter()


