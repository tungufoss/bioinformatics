import sys
sys.path.insert(0,'../')
from common import TextMatching
sys.path.insert(0,'../week2')
from BurrowsWheelerTransform import InvBurrowsWheelerTransformConstruction, BurrowsWheelerTransformConstruction

def FirstOccurrenceColumn(LastColumn):    
  FirstColumn = sorted(LastColumn)  
  FirstOccurrence = {}
  for symbol in set(LastColumn):
    FirstOccurrence[symbol] = FirstColumn.index(symbol)    
  return FirstOccurrence

def Last2FirstColumn(LastColumn):
  N = len(LastColumn)    
  SL = [(LastColumn[i], Countsymbol(LastColumn[i], i+1, LastColumn)) for i in range(N) ]
  SF = sorted(SL)  
  LastToFirst = [SF.index(SL[i]) for i in range(N)]  
  return LastToFirst

def Countsymbol(symbol, i, LastColumn):
  # Countsymbol(i, LastColumn) returns the number of occurrences of symbol in the first i positions of LastColumn
  return LastColumn[:i].count(symbol)

C = 100 # only store the Count arrays when i is divisible by C, where C is a constant; (C is typically equal to 100 in practice) 
def CountArray(LastColumn):
  symbols = sorted(list(set(LastColumn)))
  N = len(LastColumn)  
  K = len(symbols)
  cnts = {}  # these arrays are called checkpoint arrays.  
  cnts[0] = [0]*K  
  for i in range(C,N,C):    
    tmp = [LastColumn[i-C:i].count(symbol) for symbol in symbols]
    cnts[i] = [cnts[i-C][j]+tmp[j] for j in range(K)]
    
  return cnts, symbols
  
def Countsymbol2(symbol, i, LastColumn, CountSymbols, Symbols):  
  ix = Symbols.index(symbol)  
  row = (i / C)*C
  cnt = CountSymbols[row][ix]
  for j in range(row, row+i%C):
    if LastColumn[j] == symbol : 
      cnt+=1  
  return cnt
  
def BetterBWMatchingAux(FirstOccurrence, LastColumn, Pattern, CountSymbols, Symbols): 
  top = 0 # top <- 0
  bottom = len(LastColumn)-1 # bottom <- |LastColumn| - 1
  while top <= bottom : 
    if len(Pattern) > 0 : #Pattern is nonempty
      symbol = Pattern[-1] #symbol <- last letter in Pattern
      Pattern = Pattern[:-1] # remove last letter from Pattern      
      # if positions from top to bottom in LastColumn contain an occurrence of symbol
      check = LastColumn[top:bottom+1]
      if symbol in check:        
        # top <- FirstOccurrence(symbol) + Countsymbol(top, LastColumn)
        top = FirstOccurrence[symbol] + Countsymbol2(symbol, top, LastColumn, CountSymbols, Symbols) 
        # bottom <- FirstOccurrence(symbol) + Countsymbol(bottom + 1, LastColumn) - 1
        bottom = FirstOccurrence[symbol] + Countsymbol2(symbol, bottom+1, LastColumn, CountSymbols, Symbols) - 1 
      else :
        return []
    else :
      return range(top,bottom+1) # bottom - top + 1

def ApproxFirstLastPropertyMatch(LastColumn, FirstColumn, Last2First, pattern, row, d):
  mismatches = 0
  symbol = pattern[-1]
  pattern = pattern[:-1]
  
  if symbol != FirstColumn[row] : 
    mismatches += 1
  
  while len(pattern)>0:
    symbol = pattern[-1]
    pattern = pattern[:-1]
    next_char = LastColumn[row]
    if next_char == '$' : 
      return False     
    if symbol != next_char : 
      mismatches += 1    
    if mismatches > d : 
      return False      
    row = Last2First[row]

  return mismatches <= d

def SuffixArrayVal(LastColumn, Last2First, row):    
  for ix in range(len(LastColumn)) :
    if LastColumn[row] == '$': break 
    row = Last2First[row]
  return ix
  
def BetterBWMatching(LastColumn, patterns):  
  FirstOccurrence = FirstOccurrenceColumn(LastColumn)  
  CountSymbols, Symbols = CountArray(LastColumn)
  L2F = Last2FirstColumn(LastColumn)
  matches={}
  for pattern in patterns:
    rows = BetterBWMatchingAux(FirstOccurrence, LastColumn, pattern, CountSymbols, Symbols)
    matches[pattern] = [SuffixArrayVal(LastColumn, L2F, i) for i in rows]
    
  return matches
    
def ApproxBWMatching(LastColumn, patterns, d):    
  FirstOccurrence = FirstOccurrenceColumn(LastColumn)  
  CountSymbols, Symbols = CountArray(LastColumn)
  L2F = Last2FirstColumn(LastColumn)
  FirstColumn = sorted(LastColumn)
  N = len(LastColumn)
  matches={}
  for pattern in patterns:
    '''
    n = len(pattern)    
    k = n/(d+1)
    seeds = [pattern[i*k:(i+1)*k] for i in range(d+1)]
    seeds[-1] += pattern[(d+1)*k:n]        
    for seed in seeds:
      # check exact match for each seed
      seed_range = BetterBWMatchingAux(FirstOccurrence, LastColumn, seed, CountSymbols, Symbols)      
      if len(seed_range) == 0 : continue      
      print '{}->{} {}'.format(pattern,seeds, seed)
    '''
    rows = []
    for row in range(N):
      if ApproxFirstLastPropertyMatch(LastColumn, FirstColumn, L2F, pattern, row, d) :
        rows.append(row)          
    matches[pattern] = [SuffixArrayVal(LastColumn, L2F, i)-len(pattern)+1 for i in rows]   
  return matches
    
def main_BetterBWMatching(myfile):    
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  
  with open(inputFile) as inFile:        
    transform = inFile.readline().strip()
    patterns = [x for x in inFile.readline().strip().split(' ')]
  
  matches = BetterBWMatching(transform, patterns)
    
  with open(outputFile, 'w') as outFile:        
    txt = ' '.join([str(len(matches[pattern])) for pattern in patterns])
    outFile.write(txt)  
    print txt    

def main_MultiplePatternMatching(myfile):    
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  
  with open(inputFile) as inFile:        
    text = inFile.readline().strip()
    patterns=[]
    for line in inFile :
      patterns.append(line.strip())
  
  if text[-1] != '$' :
    text += '$'
  
  transform = BurrowsWheelerTransformConstruction(text)
  matches = BetterBWMatching(transform, patterns)
    
  with open(outputFile, 'w') as outFile:
    indices = []
    for pattern in patterns:
      indices += matches[pattern]  
    txt = ' '.join([str(ix) for ix in sorted(indices)])
    outFile.write(txt)  
    print txt    

def main_MultipleApproxPatternMatching(myfile):    
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  
  with open(inputFile) as inFile:        
    text = inFile.readline().strip()
    patterns = [x for x in inFile.readline().strip().split(' ')]
    d = int(inFile.readline().strip())      
  
  if text[-1] != '$' :
    text += '$'
  
  transform = BurrowsWheelerTransformConstruction(text)
  matches = ApproxBWMatching(transform, patterns, d)
    
  with open(outputFile, 'w') as outFile:
    indices = []
    for pattern in patterns:
      indices += matches[pattern]  
    txt = ' '.join([str(ix) for ix in sorted(indices)])
    outFile.write(txt)    
    print txt
    
'''    
main_BetterBWMatching('sample_BetterBWMatching')
main_BetterBWMatching('dataset_301_7')
main_MultiplePatternMatching('sample_MultiplePatternMatching')
main_MultiplePatternMatching('dataset_303_4')
'''
main_MultipleApproxPatternMatching('sample_MultipleApproxPatternMatching')
main_MultipleApproxPatternMatching('dataset_304_6')