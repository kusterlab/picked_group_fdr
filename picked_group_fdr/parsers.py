import sys
import csv
import re
from typing import List, Dict, Optional
import logging

logger = logging.getLogger(__name__)


#csv.field_size_limit(sys.maxsize)
csv.field_size_limit(2147483647)


TMT_UNIMOD = "[UNIMOD:737]" # TMT6/10/11
TMTPRO_UNIMOD = "[UNIMOD:2016]"
ITRAQ4_UNIMOD = "[UNIMOD:214]"
ITRAQ8_UNIMOD = "[UNIMOD:730]"



# copied from fundamentals/constants.py
MAXQUANT_VAR_MODS = {
    "(ox)": "[UNIMOD:35]",
    "(Oxidation (M))": "[UNIMOD:35]",
    "(tm)": "[UNIMOD:737]",
    "_(tm)": f"_{TMT_UNIMOD}",
    "K(tm)": f"K{TMT_UNIMOD}",
    "_(TMTPro (N-term))": f"_{TMTPRO_UNIMOD}",
    "K(TMTPro (K))": f"K{TMTPRO_UNIMOD}",
    "_(iTRAQ4plex (N-term))": f"_{ITRAQ4_UNIMOD}",
    "K(iTRAQ4plex (K))": f"K{ITRAQ4_UNIMOD}",
    "_(iTRAQ8plex (N-term))": f"_{ITRAQ8_UNIMOD}",
    "K(iTRAQ8plex (K))": f"K{ITRAQ8_UNIMOD}",
    "(ph)": "[UNIMOD:21]",
    "K(Lys8)": "K[UNIMOD:259]",
    "R(Arg10)": "R[UNIMOD:267]",
    "C(Carbamidomethyl (C))": "C[UNIMOD:4]",
}

DEFAULT_FIXED_MODS = {'C': 'C[UNIMOD:4]'}

SILAC_HEAVY_FIXED_MODS = {'C': 'C[UNIMOD:4]',
                          'K': 'K[UNIMOD:259]', 
                          'R': 'R[UNIMOD:267]'}

TMT_FIXED_MODS = {'C': 'C[UNIMOD:4]',
                  '^_':f"_{TMT_UNIMOD}", 
                  'K': f"K{TMT_UNIMOD}"}

TMTPRO_FIXED_MODS = {'C': 'C[UNIMOD:4]',
                     '^_':f"_{TMTPRO_UNIMOD}", 
                     'K': f"K{TMTPRO_UNIMOD}"}

ITRAQ4_FIXED_MODS = {'C': 'C[UNIMOD:4]',
                     '^_':f"_{ITRAQ4_UNIMOD}", 
                     'K': f"K{ITRAQ4_UNIMOD}"}

ITRAQ8_FIXED_MODS = {'C': 'C[UNIMOD:4]',
                     '^_':f"_{ITRAQ8_UNIMOD}", 
                     'K': f"K{ITRAQ8_UNIMOD}"}

FIXED_MODS_UNIMOD = [TMT_UNIMOD, TMTPRO_UNIMOD, ITRAQ4_UNIMOD, ITRAQ8_UNIMOD]
FIXED_MODS_DICTS = [DEFAULT_FIXED_MODS, TMT_FIXED_MODS, TMTPRO_FIXED_MODS, ITRAQ4_FIXED_MODS, ITRAQ8_FIXED_MODS]

def parseMqProteinGroupsFile(mqProteinGroupsFile):
  if mqProteinGroupsFile.endswith('.csv'):
    delimiter = ','
  else:
    delimiter = '\t'
    
  reader = getTsvReader(mqProteinGroupsFile, delimiter)
  headers = next(reader) # save the header
  
  scoreCol = headers.index('Score')
  proteinCol = headers.index('Protein IDs')
  #proteinCol = headers.index('Majority protein IDs')
  
  logger.info("Parsing MaxQuant proteinGroups.txt file")
  for row in reader:
    if len(row[scoreCol]) == 0:
      yield row[proteinCol].split(";"), -100.0
    else:
      yield row[proteinCol].split(";"), float(row[scoreCol])


def parseMqEvidenceFile(mqEvidenceFile, scoreType, forQuantification = False):
  """
  Reads in approximately 100,000 lines per second with forQuantification=False 
  and 50,000 lines per second with forQuantification=True
  """
  if mqEvidenceFile.endswith('.csv'):
    delimiter = ','
  else:
    delimiter = '\t'
  reader = getTsvReader(mqEvidenceFile, delimiter)
  headers = list(map(lambda x : x.lower(), next(reader))) # save the header
  
  if "psmid" in headers:
    yield from parsePercolatorOutFile(mqEvidenceFile, scoreType)
  else:
    if mqEvidenceFile.endswith('.csv'):
      headers = [x.replace(".", " ") for x in headers]
    getHeaderCol = getHeaderColFunc(headers)
    getHeaderColsStartingWith = getHeaderColsStartingWithFunc(headers)
    
    peptCol = getHeaderCol('modified sequence', required = True)
    
    proteinCol = getHeaderCol('leading proteins', required = True) # all protein groups, each represented by the first protein in the group
    if scoreType.use_razor:
      proteinCol = getHeaderCol('leading razor protein', required = True) # best scoring protein group, represented by the first protein in the group      
    
    scoreCol = getHeaderCol(scoreType.get_score_column(), required = True)
        
    experimentCol = getHeaderCol('experiment')
    chargeCol = getHeaderCol('charge', required = forQuantification)
    
    intensityCol = getHeaderCol('intensity', required = forQuantification)
    tmtCols = getHeaderColsStartingWith('reporter intensity ')
    silacCols = list()
    if 'intensity l' in headers: # SILAC
      silacCols.append(getHeaderCol('intensity l', required = forQuantification))
      if 'intensity m' in headers:
        silacCols.append(getHeaderCol('intensity m', required = forQuantification))
      if 'intensity h' in headers:
        silacCols.append(getHeaderCol('intensity h', required = forQuantification))
    
    rawFileCol = getHeaderCol('raw file', required = forQuantification)
    fractionCol = getHeaderCol('fraction')
    evidenceIdCol = getHeaderCol('id', required = forQuantification)
    
    logger.info("Parsing MaxQuant evidence.txt file")
    for lineIdx, row in enumerate(reader):
      if lineIdx % 500000 == 0:
        logger.info(f"  Reading line {lineIdx}")
      
      peptide = row[peptCol]
      proteins = row[proteinCol]
      if experimentCol >= 0:
        experiment = row[experimentCol]
      else:
        experiment = "Experiment1"
      
      score = float(row[scoreCol])
      
      if forQuantification:
        charge = int(row[chargeCol])
        intensity = float(row[intensityCol]) if len(row[intensityCol]) > 0 else 0.0
        if fractionCol >= 0:
          fraction = row[fractionCol]
        else:
          fraction = -1
        rawFile = row[rawFileCol]
        tmtIntensities = [row[tmtCol] for tmtCol in tmtCols]
        silacIntensities = [row[silacCol] if len(row[silacCol]) > 0 else 0 for silacCol in silacCols]
        evidenceId = int(row[evidenceIdCol])
        yield peptide, charge, rawFile, experiment, fraction, intensity, score, tmtIntensities, silacIntensities, evidenceId
      else:
        yield peptide, proteins.split(";"), experiment, score


def getHeaderColFunc(headers):
  def getHeaderCol(name, required = False):
    if required:
      return headers.index(name)
    else:
      if name in headers:
        return headers.index(name)
      return -1
  return getHeaderCol

def getHeaderColsStartingWithFunc(headers):
  def getHeaderColsStartingWith(name):
    cols = [idx for idx, h in enumerate(headers) if h.startswith(name)]
    return cols
  return getHeaderColsStartingWith


def parsePercolatorOutFile(percOutFile, scoreType = "PEP", razor = False):
  if percOutFile.endswith('.csv'):
    delimiter = ','
  else:
    delimiter = '\t'
  reader = getTsvReader(percOutFile, delimiter)
  headers = next(reader) # save the header
  
  peptCol = headers.index('peptide')
  proteinCol = headers.index('proteinIds')
  scoreCol = headers.index(scoreType.get_score_column())
  
  logger.info("Parsing Percolator output file")
  for lineIdx, row in enumerate(reader):
    if lineIdx % 500000 == 0:
      logger.info(f"  Reading line {lineIdx}")
    
    peptide = row[peptCol][1:-1]
    proteins = row[proteinCol:]
    #proteins = list(map(lambda x : "REV__" + x.split("|")[1] if x.startswith("REV__") else x.split("|")[1], proteins)) # for mouse proteome
    experiment = 1
    score = float(row[scoreCol])
    
    yield peptide, proteins, experiment, score


def parsePercolatorOutFileToDict(percOutFile, resultsDict, inputType = ""):
  if percOutFile.endswith('.csv'):
    delimiter = ','
  else:
    delimiter = '\t'
  reader = getTsvReader(percOutFile, delimiter)
  headers = next(reader) # save the header
  
  if 'PSMId' in headers:
    idCol = headers.index('PSMId')
    peptCol = headers.index('peptide')
    scoreCol = headers.index('score')
    postErrProbCol = headers.index('posterior_error_prob')
  else: # mokapot results
    idCol = headers.index('SpecId')
    peptCol = headers.index('Peptide')
    scoreCol = headers.index('mokapot score')
    postErrProbCol = headers.index('mokapot PEP')
  
  logger.info("Parsing Percolator output file")
  fixed_mod_idx = -1
  first = True
  for lineIdx, row in enumerate(reader):
    if lineIdx % 500000 == 0:
      logger.info(f"  Reading line {lineIdx}")
    
    peptide = row[peptCol][2:-2]
    score = float(row[scoreCol])
    postErrProb = float(row[postErrProbCol])
    
    psmId = row[idCol]
    if inputType == "prosit":
      peptide = peptide.replace('m', 'M(ox)')
      if first:
        for i, fixed_mod in enumerate(FIXED_MODS_UNIMOD):
          if fixed_mod in peptide:
            fixed_mod_idx = i
        first = False
      elif fixed_mod_idx >= 0:
        if FIXED_MODS_UNIMOD[fixed_mod_idx] not in peptide:
          fixed_mod_idx = -1
      rawFile = "-".join(psmId.split("-")[:-4])
      scanNumber = int(float(psmId.split("-")[-4]))
    else:
      peptide = peptide.replace('[42]', '(ac)').replace('M[16]', 'M(ox)')
      rawFile = "_".join(psmId.split("_")[:-3])
      scanNumber = int(psmId.split("_")[-3])

    resultsDict[rawFile][(scanNumber, peptide)] = (score, postErrProb)
  
  return FIXED_MODS_DICTS[fixed_mod_idx+1], resultsDict


# copied from fundamentals/mod_string.py
def maxquant_to_internal(
        sequences: List[str],
        fixed_mods: Optional[Dict[str, str]] = {'C': 'C[UNIMOD:4]'}
) -> List[str]:
    """
    Function to translate a MaxQuant modstring to the Prosit format
    :param sequences: List[str] of sequences
    :param fixed_mods: Optional dictionary of modifications with key aa and value mod, e.g. 'M': 'M(UNIMOD:35)'.
    Fixed modifications must be included in the variable modificatons dictionary throws Assertion error otherwise.
    :return: List[str] of modified sequences.
    """
    err_msg = f"Provided illegal fixed mod, supported modifications are {set(MAXQUANT_VAR_MODS.values())}."
    assert all(x in MAXQUANT_VAR_MODS.values() for x in fixed_mods.values()), err_msg
    
    replacements = {**MAXQUANT_VAR_MODS, **fixed_mods}

    def custom_regex_escape(key: str) -> str:
        """
        Subfunction to escape only normal brackets in the modstring
        :param key: The match to escape.
        :return match with escaped special characters.
        """
        for k, v in {"(": "\(", ")": "\)"}.items():
            key = key.replace(k, v)
        return key

    regex = re.compile("|".join(map(custom_regex_escape, replacements.keys())))

    def find_replacement(match: re) -> str:
        """
        Subfunction to find the corresponding substitution for a match.
        :param match: an re.Match object found by re.sub
        :return substitution string for the given match
        """
        key = match.string[match.start():match.end()]
        if "_" in key:  # If _ is in the match we need to differentiate n and c term
            if match.start() == 0:
                key = f"^{key}"
            else:
                key = f"{key}$"

        return replacements[key]

    return [regex.sub(find_replacement, seq)[1:-1] for seq in sequences]


def getTsvReader(filename, delimiter = '\t'):
  # Python 3
  if sys.version_info[0] >= 3:
    return csv.reader(open(filename, 'r', newline = ''), delimiter = delimiter)
  # Python 2
  else:
    return csv.reader(open(filename, 'rb'), delimiter = delimiter)


def getTsvWriter(filename, delimiter = '\t'):
  # Python 3
  if sys.version_info[0] >= 3:
    return csv.writer(open(filename, 'w', newline = ''), delimiter = delimiter)
  # Python 2
  else:
    return csv.writer(open(filename, 'wb'), delimiter = delimiter)

