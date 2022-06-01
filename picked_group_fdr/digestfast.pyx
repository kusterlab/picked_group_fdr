# cython: linetrace=True
# cython: binding=True
# distutils: language = c++
# distutils: define_macros=CYTHON_TRACE=1
cdef extern from "<algorithm>" namespace "std":
    Iter find[Iter, T](Iter first, Iter last, const T& value)
    

# cythonize -a -i digestfast.pyx

from libcpp.vector cimport vector
from libcpp.set cimport set as c_set

import itertools
import sys

# method for generating the list of peptides
def getDigestedPeptides(seq, int min_len = 6, int max_len = 50, pre2 = ['K', 'R'], not_post2 = ['P'], digestion = 'full', int miscleavages = 0, methionineCleavage = True):
  cdef vector[int] cleavageSites, starts
  cdef int i, lenS, lenP, start, methionineCleaved, one, c
  cdef c_set[int] pre, not_post
  
  one = 1
  
  for p in pre2:
    pre.insert(ord(p))
    
  for p in not_post2:
    not_post.insert(ord(p))
  
  lenS = len(seq)
  starts.push_back(0)
  length_accepted = lambda x : x >= min_len and x <= max_len
  methionineCleavage = methionineCleavage and seq[0] == "M"
  
  if digestion == 'none':
    for i in range(lenS + 1):
      for j in range(i + min_len, i + max_len + 1):
        if j < lenS:
          yield seq[i:j]
  elif digestion == 'semi':
    for i in range(lenS + 1):
      isCleavageSite = (seq[min([lenS-1,i])] in pre2 and seq[min([lenS-1,i+1])] not in not_post2)
      isMethionineCleavageSite = (i == 0 and methionineCleavage)
      if i == lenS or isCleavageSite or isMethionineCleavageSite:
        # peptides with enzymatic C-terminal (both enzymatic and non-enzymatic N-terminal)
        start = starts[0]
        for j in range(start, min([i+1, lenS])):
          lenP = min([i, lenS - 1]) - j + 1
          if length_accepted(lenP):
            yield (seq[j : i + 1])
        starts.push_back(i + 1)
        methionineCleaved = int(methionineCleavage and starts[0] == 0)
        if len(starts) > miscleavages + 1 + methionineCleaved or i == lenS:    
          starts = starts[1 + methionineCleaved:]
      else: # peptides with non enzymatic C-terminal
        for start in starts:
          lenP = i - start + 1
          if length_accepted(lenP) and starts[-1] == i + 1:
            yield (seq[start : i + 1])
  else:
    hasMethionineCleavage = False
    if methionineCleavage and seq[0] == 'M':
      cleavageSites.push_back(0)
      hasMethionineCleavage = True
    for i in range(lenS):
      if pre.find(ord(seq[i])) != pre.end() and not_post.find(ord(seq[min([lenS-1,i+1])])) == not_post.end():
        cleavageSites.push_back(i)
    cleavageSites.push_back(lenS)
    for i in cleavageSites:
      for start in starts:
        lenP = i - start + 1
        if lenP >= min_len and lenP <= max_len:
          yield (seq[start : i + 1])
      starts.push_back(i + 1)
      methionineCleaved = int(starts[0] == 0 and hasMethionineCleavage)
      if starts.size() > miscleavages + one + methionineCleaved:   
        starts = starts[one + methionineCleaved:]

def readFastaMaxQuant(filePath, db = "target", parseId = lambda x : x.split(" ")[0], specialAAs = ['R', 'K']):
  cdef int i
  if db not in ["target", "decoy", "concat"]:
    sys.exit("unknown db mode: %s" % db)
  
  hasSpecialAAs = len(specialAAs) > 0
  name, seq = None, []
  with open(filePath, 'r') as fp:
    for line in itertools.chain(fp, [">"]):
      line = line.rstrip()
      if line.startswith(">"):
        if name: 
          if db in ["target", "concat"]:
            yield (name, "".join(seq))
          
          if db in ["decoy", "concat"]:
            seq = list("".join(seq)[::-1])
            if hasSpecialAAs:
              for i in range(1, len(seq)):
                if seq[i] in specialAAs:
                  swapPositions(seq, i, i-1)
            yield ("REV__" + name, "".join(seq))
          
        if len(line) > 1:
          name, seq = parseId(line[1:]), []
      else: seq.append(line)

def swapPositions(seq, pos1, pos2):     
  seq[pos1], seq[pos2] = seq[pos2], seq[pos1] 
