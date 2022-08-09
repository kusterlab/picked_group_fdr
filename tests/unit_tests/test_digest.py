import pytest

import picked_group_fdr.digest as digest

def test_swapSpecialAAs():
    seq = 'ABCKDEFRRR'
    specialAAs = ['K', 'R']
    assert digest.swapSpecialAAs(seq, specialAAs) == 'ABKCDERRRF'


def test_hasMiscleavage_true():
    seq = 'ABCDEFKR'
    assert digest.hasMiscleavage(seq) == True


def test_hasMiscleavage_false():
    seq = 'ABCDEFKPR'
    assert digest.hasMiscleavage(seq) == False


def test_hasGeneNames_true():
    proteinAnnotations = dict()
    for p in ['A', 'B', 'C', 'D', 'E']:
        proteinAnnotations[f'protein{p}'] = (f'protein{p}', f'gene{p}', f'fasta{p}')
    assert digest.hasGeneNames(proteinAnnotations, minRatioWithGenes=0.5) == True


def test_hasGeneNames_false():
    proteinAnnotations = dict()
    for p in ['A', 'B', 'C', 'D', 'E']:
        proteinAnnotations[f'protein{p}'] = (f'protein{p}', "", f'fasta{p}')
    assert digest.hasGeneNames(proteinAnnotations, minRatioWithGenes=0.5) == False


def test_hasGeneNames_false_below_half():
    proteinAnnotations = dict()
    for p in ['A', 'B', 'C', 'D']:
        proteinAnnotations[f'protein{p}'] = (f'protein{p}', f'gene{p}', f'fasta{p}')
    for p in ['A', 'B', 'C', 'D', 'E']:
        proteinAnnotations[f'protein_no_gene{p}'] = (f'protein_no_gene{p}', "", f'fasta_no_gene{p}')
    
    assert digest.hasGeneNames(proteinAnnotations, minRatioWithGenes=0.5) == False


class TestNonSpecificDigest:
    def test_nonSpecificDigest(self):
        seq = 'ABCDEFGH'
        min_len = 6
        max_len = 30
        assert set(digest.nonSpecificDigest(seq, min_len, max_len)) == set(['ABCDEF', 'BCDEFG', 'CDEFGH', 'ABCDEFG', 'BCDEFGH', 'ABCDEFGH'])

    def test_nonSpecificDigest_maxLen(self):
        seq = 'ABCDEFGH'
        min_len = 6
        max_len = 7
        assert set(digest.nonSpecificDigest(seq, min_len, max_len)) == set(['ABCDEF', 'BCDEFG', 'CDEFGH', 'ABCDEFG', 'BCDEFGH'])


class TestSemiSpecificDigest:
    def test_semiSpecificDigest_noCleavageSite(self):
        seq = 'ABCDEFGH'
        min_len = 6
        max_len = 30
        pre = ['K', 'R']
        not_post = ['P']
        miscleavages = 2
        methionineCleavage = True
        assert set(digest.semiSpecificDigest(seq, min_len, max_len, pre, not_post, miscleavages, methionineCleavage)) == set(
            ['ABCDEF', 
             'ABCDEFG', 
             'ABCDEFGH', 
                'BCDEFGH', 
                 'CDEFGH'])
    
    def test_semiSpecificDigest_methionineCleavage(self):
        seq = 'MABCDEFGH'
        min_len = 6
        max_len = 30
        pre = ['K', 'R']
        not_post = ['P']
        miscleavages = 2
        methionineCleavage = True
        assert set(digest.semiSpecificDigest(seq, min_len, max_len, pre, not_post, miscleavages, methionineCleavage)) == set(
            ['MABCDE',
             'MABCDEF', 
             'MABCDEFG', 
             'MABCDEFGH', 
                'ABCDEF', 
                'ABCDEFG', 
                'ABCDEFGH', 
                 'BCDEFGH', 
                    'CDEFGH'])
    
    # make sure that the methionine cleavage is not counted as a miscleavage
    def test_semiSpecificDigest_methionineCleavagePlusOneMiscleavage(self):
        seq = 'MABCDEFKKK'
        min_len = 6
        max_len = 30
        pre = ['K', 'R']
        not_post = ['P']
        miscleavages = 1
        methionineCleavage = True
        assert set(digest.semiSpecificDigest(seq, min_len, max_len, pre, not_post, miscleavages, methionineCleavage)) == set(
            ['MABCDE',
             'MABCDEF',
             'MABCDEFK',
             'MABCDEFKK', 
                'ABCDEF',
                'ABCDEFK',
                'ABCDEFKK', 
                 'BCDEFK',
                 'BCDEFKK',
                    'CDEFKK'])
        
    def test_semiSpecificDigest_noCleavageSiteMaxLen(self):
        seq = 'ABCDEFGH'
        min_len = 6
        max_len = 7
        pre = ['K', 'R']
        not_post = ['P']
        miscleavages = 2
        methionineCleavage = True
        assert set(digest.semiSpecificDigest(seq, min_len, max_len, pre, not_post, miscleavages, methionineCleavage)) == set(
            ['ABCDEF', 
             'ABCDEFG', 
                'BCDEFGH', 
                 'CDEFGH'])
    
    def test_semiSpecificDigest_noMiscleavage(self):
        seq = 'ABCDEFGKX'
        min_len = 6
        max_len = 30
        pre = ['K', 'R']
        not_post = ['P']
        miscleavages = 0
        methionineCleavage = True
        assert set(digest.semiSpecificDigest(seq, min_len, max_len, pre, not_post, miscleavages, methionineCleavage)) == set(
            ['ABCDEF', 
             'ABCDEFG', 
             'ABCDEFGK', 
                'BCDEFGK', 
                 'CDEFGK'])
    
    def test_semiSpecificDigest_oneMiscleavage(self):
        seq = 'ABCDEFGKX'
        min_len = 6
        max_len = 30
        pre = ['K', 'R']
        not_post = ['P']
        miscleavages = 1
        methionineCleavage = True
        assert set(digest.semiSpecificDigest(seq, min_len, max_len, pre, not_post, miscleavages, methionineCleavage)) == set(
            ['ABCDEF', 
             'ABCDEFG', 
             'ABCDEFGK',
                'BCDEFGK',
                 'CDEFGK',
             'ABCDEFGKX',
                'BCDEFGKX', 
                 'CDEFGKX',
                  'DEFGKX'])
                 
    def test_semiSpecificDigest_notPost(self):
        seq = 'ABCDEFKPX'
        min_len = 6
        max_len = 30
        pre = ['K', 'R']
        not_post = ['P']
        miscleavages = 0
        methionineCleavage = True
        assert set(digest.semiSpecificDigest(seq, min_len, max_len, pre, not_post, miscleavages, methionineCleavage)) == set(
            ['ABCDEF', 
             'ABCDEFK', 
             'ABCDEFKP', 
             'ABCDEFKPX', 
                'BCDEFKPX', 
                 'CDEFKPX', 
                  'DEFKPX'])


class TestFullDigest:
    def test_fullDigest_noCleavageSite(self):
        seq = 'ABCDEFGH'
        min_len = 6
        max_len = 30
        pre = ['K', 'R']
        not_post = ['P']
        miscleavages = 2
        methionineCleavage = True
        assert set(digest.fullDigest(seq, min_len, max_len, pre, not_post, miscleavages, methionineCleavage)) == set(
            ['ABCDEFGH'])
    
    def test_fullDigest_methionineCleavage(self):
        seq = 'MABCDEFGH'
        min_len = 6
        max_len = 30
        pre = ['K', 'R']
        not_post = ['P']
        miscleavages = 2
        methionineCleavage = True
        assert set(digest.fullDigest(seq, min_len, max_len, pre, not_post, miscleavages, methionineCleavage)) == set(
            ['MABCDEFGH',
              'ABCDEFGH'])
    
    # make sure that the methionine cleavage is not counted as a miscleavage
    def test_fullDigest_methionineCleavageOneMiscleavage(self):
        seq = 'MABCDEFGHKKK'
        min_len = 6
        max_len = 30
        pre = ['K', 'R']
        not_post = ['P']
        miscleavages = 1
        methionineCleavage = True
        assert set(digest.fullDigest(seq, min_len, max_len, pre, not_post, miscleavages, methionineCleavage)) == set(
            ['MABCDEFGHK',
             'MABCDEFGHKK',
              'ABCDEFGHK',
              'ABCDEFGHKK'])
    
    def test_fullDigest_noCleavageSiteMaxLen(self):
        seq = 'ABCDEFGH'
        min_len = 6
        max_len = 7
        pre = ['K', 'R']
        not_post = ['P']
        miscleavages = 2
        methionineCleavage = True
        assert set(digest.fullDigest(seq, min_len, max_len, pre, not_post, miscleavages, methionineCleavage)) == set(
            [])
    
    def test_fullDigest_noMiscleavage(self):
        seq = 'ABCDEFGKX'
        min_len = 6
        max_len = 30
        pre = ['K', 'R']
        not_post = ['P']
        miscleavages = 0
        methionineCleavage = True
        assert set(digest.fullDigest(seq, min_len, max_len, pre, not_post, miscleavages, methionineCleavage)) == set(
            ['ABCDEFGK'])
    
    def test_fullDigest_oneMiscleavage(self):
        seq = 'ABCDEFGKX'
        min_len = 6
        max_len = 30
        pre = ['K', 'R']
        not_post = ['P']
        miscleavages = 1
        methionineCleavage = True
        assert set(digest.fullDigest(seq, min_len, max_len, pre, not_post, miscleavages, methionineCleavage)) == set(
            ['ABCDEFGK',
             'ABCDEFGKX'])
    
    def test_fullDigest_notPost(self):
        seq = 'ABCDEFKPX'
        min_len = 6
        max_len = 30
        pre = ['K', 'R']
        not_post = ['P']
        miscleavages = 0
        methionineCleavage = True
        assert set(digest.fullDigest(seq, min_len, max_len, pre, not_post, miscleavages, methionineCleavage)) == set(
            ['ABCDEFKPX'])

