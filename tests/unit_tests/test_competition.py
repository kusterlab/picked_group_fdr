import pytest

import picked_group_fdr.competition as competition
import picked_group_fdr.protein_groups as protein_groups
from picked_group_fdr.scoring import BestPEPScore


class TestPickedStrategy:
  def test_do_competition(self, proteinGroupsSimple, proteinGroupScoresSimple):
    score_type = BestPEPScore()
    picked_strategy = competition.PickedStrategy()
    filtered_protein_groups, filtered_protein_group_scores = picked_strategy.do_competition(proteinGroupsSimple, proteinGroupScoresSimple, score_type)
    
    assert filtered_protein_groups.get_all_proteins() == set(['proteinA', 'proteinE', 'REV__proteinA', 'REV__proteinC'])
    assert filtered_protein_group_scores == [[(0.001, "", [])], [(0.01, "", [])]]
  
  def test_do_competition_no_grouping(self, proteinGroupsNoGrouping, proteinGroupScoresNoGrouping):
    score_type = BestPEPScore()
    picked_strategy = competition.PickedStrategy()
    filtered_protein_groups, filtered_protein_group_scores = picked_strategy.do_competition(proteinGroupsNoGrouping, proteinGroupScoresNoGrouping, score_type)
    
    assert filtered_protein_groups.get_all_proteins() == set(['proteinA', 'proteinE', 'REV__proteinC'])
    assert filtered_protein_group_scores == [[(0.001, "", [])], [(0.001, "", [])], [(0.01, "", [])]]
  
  def test_do_competition_unordered(self, proteinGroupsUnordered, proteinGroupScoresUnordered):
    score_type = BestPEPScore()
    picked_strategy = competition.PickedStrategy()
    filtered_protein_groups, filtered_protein_group_scores = picked_strategy.do_competition(proteinGroupsUnordered, proteinGroupScoresUnordered, score_type)
    
    assert filtered_protein_groups.get_all_proteins() == set(['proteinA', 'proteinE', 'REV__proteinA', 'REV__proteinC'])
    assert filtered_protein_group_scores == [[(0.001, "", [])], [(0.01, "", [])]]


class TestPickedGroupStrategy:
  def test_do_competition(self, proteinGroupsSimple, proteinGroupScoresSimple):
    score_type = BestPEPScore()
    picked_strategy = competition.PickedGroupStrategy()
    filtered_protein_groups, filtered_protein_group_scores = picked_strategy.do_competition(proteinGroupsSimple, proteinGroupScoresSimple, score_type)
    
    assert filtered_protein_groups.get_all_proteins() == set(['proteinA', 'proteinE'])
    assert filtered_protein_group_scores == [[(0.001, "", [])]]
  
  def test_do_competition_no_grouping(self, proteinGroupsNoGrouping, proteinGroupScoresNoGrouping):
    score_type = BestPEPScore()
    picked_strategy = competition.PickedGroupStrategy()
    filtered_protein_groups, filtered_protein_group_scores = picked_strategy.do_competition(proteinGroupsNoGrouping, proteinGroupScoresNoGrouping, score_type)
    
    assert filtered_protein_groups.get_all_proteins() == set(['proteinA', 'proteinE', 'REV__proteinC'])
    assert filtered_protein_group_scores == [[(0.001, "", [])], [(0.001, "", [])], [(0.01, "", [])]]

  def test_do_competition_unordered(self, proteinGroupsUnordered, proteinGroupScoresUnordered):
    score_type = BestPEPScore()
    picked_strategy = competition.PickedGroupStrategy()
    filtered_protein_groups, filtered_protein_group_scores = picked_strategy.do_competition(proteinGroupsUnordered, proteinGroupScoresUnordered, score_type)
    
    assert filtered_protein_groups.get_all_proteins() == set(['REV__proteinA', 'REV__proteinC'])
    assert filtered_protein_group_scores == [[(0.001, "", [])]]
    

class TestClassicStrategy:
  def test_do_competition(self, proteinGroupsSimple, proteinGroupScoresSimple):
    score_type = BestPEPScore()
    classic_strategy = competition.ClassicStrategy()
    filtered_protein_groups, filtered_protein_group_scores = classic_strategy.do_competition(proteinGroupsSimple, proteinGroupScoresSimple, score_type)
    
    assert filtered_protein_groups.get_all_proteins() == set(['proteinA', 'proteinE', 'REV__proteinA', 'REV__proteinC'])
    assert filtered_protein_group_scores == [[(0.001, "", [])], [(0.01, "", [])]]
  
  def test_do_competition_no_grouping(self, proteinGroupsNoGrouping, proteinGroupScoresNoGrouping):
    score_type = BestPEPScore()
    picked_strategy = competition.ClassicStrategy()
    filtered_protein_groups, filtered_protein_group_scores = picked_strategy.do_competition(proteinGroupsNoGrouping, proteinGroupScoresNoGrouping, score_type)
    
    assert filtered_protein_groups.get_all_proteins() == set(['proteinA', 'proteinE', 'REV__proteinA', 'REV__proteinC'])
    assert filtered_protein_group_scores == [[(0.001, "", [])], [(0.001, "", [])], [(0.01, "", [])], [(0.01, "", [])]]

  def test_do_competition_unordered(self, proteinGroupsUnordered, proteinGroupScoresUnordered):
    score_type = BestPEPScore()
    classic_strategy = competition.ClassicStrategy()
    filtered_protein_groups, filtered_protein_group_scores = classic_strategy.do_competition(proteinGroupsUnordered, proteinGroupScoresUnordered, score_type)
    
    assert filtered_protein_groups.get_all_proteins() == set(['proteinA', 'proteinE', 'REV__proteinA', 'REV__proteinC'])
    assert filtered_protein_group_scores == [[(0.001, "", [])], [(0.01, "", [])]]


@pytest.fixture
def proteinGroupsSimple():
  proteinGroups = protein_groups.ProteinGroups([
          ['proteinA', 'proteinE'], 
          ['REV__proteinA', 'REV__proteinC']])
  proteinGroups.create_index()
  return proteinGroups


@pytest.fixture
def proteinGroupScoresSimple():
  # tuples of (score, peptide, proteins)
  return [[(0.001, "", [])],
          [(0.01, "", [])]]


@pytest.fixture
def proteinGroupsNoGrouping():
  proteinGroups = protein_groups.ProteinGroups([
          ['proteinA'], 
          ['proteinE'],
          ['REV__proteinA'],
          ['REV__proteinC']])
  proteinGroups.create_index()
  return proteinGroups


@pytest.fixture
def proteinGroupScoresNoGrouping():
  # tuples of (score, peptide, proteins)
  return [[(0.001, "", [])],
          [(0.001, "", [])],
          [(0.01, "", [])],
          [(0.01, "", [])]]


@pytest.fixture
def proteinGroupsUnordered():
  proteinGroups = protein_groups.ProteinGroups([
          ['proteinA', 'proteinE'], 
          ['REV__proteinA', 'REV__proteinC']])
  proteinGroups.create_index()
  return proteinGroups


@pytest.fixture
def proteinGroupScoresUnordered():
  # tuples of (score, peptide, proteins)
  return [[(0.01, "", [])],
          [(0.001, "", [])]]
