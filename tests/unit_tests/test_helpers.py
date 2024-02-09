import pytest
import numpy as np

import picked_group_fdr.helpers as helpers


def test_is_contaminant_true():
  assert helpers.is_contaminant(["CON__contaminantA"]) == True

  
def test_is_contaminant_false():
  assert helpers.is_contaminant(["CON__contaminantA", "NotAContaminant"]) == False


def test_is_decoy_true():
  assert helpers.is_decoy(["REV__decoyA"]) == True

  
def test_is_decoy_false():
  assert helpers.is_decoy(["REV__decoyA", "proteinA"]) == False


def test_is_shared_peptide_true():
  assert helpers.is_shared_peptide({1,2,3}) == True


def test_is_shared_peptide_false():
  assert helpers.is_shared_peptide({3}) == False


def test_is_shared_peptide_false():
  assert helpers.is_shared_peptide({-1}) == True
  

def test_remove_decoy_proteins_from_target_peptides():
  assert helpers.remove_decoy_proteins_from_target_peptides(["proteinA", "REV__proteinA"]) == ["proteinA"]
  

class TestCleanPeptide:
  def test_clean_peptide_flanks(self):
    assert helpers.clean_peptide("_APEPTIDE_") == "APEPTIDE"
  
  def test_clean_peptide_mods_square_brackets(self):
    assert helpers.clean_peptide("_APEPT[AMODIFICATION]IDE_") == "APEPTIDE"
  
  def test_clean_peptide_mods_parentheses(self):
    assert helpers.clean_peptide("_APEPT(AMODIFICATION)IDE_") == "APEPTIDE"
  
  def test_clean_peptide_mods_nested_parentheses(self):
    assert helpers.clean_peptide("_APEPT(Oxidation (M))IDE_") == "APEPTIDE"


def test_chunks():
  assert list(helpers.chunks([1,2,3,4,5,6,7,8,9,10], 3)) == [[1,2,3], [4,5,6], [7,8,9], [10]]

