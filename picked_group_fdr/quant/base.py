from abc import ABC, abstractmethod


class ProteinGroupColumns(ABC):
  @abstractmethod
  def append_headers(self, proteinGroupResults, experiments):
    pass
  
  @abstractmethod
  def append_columns(self, proteinGroupResults, experimentToIdxMap, postErrProbCutoff):
    pass

