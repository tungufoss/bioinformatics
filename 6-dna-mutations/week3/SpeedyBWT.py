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
  
def BetterBWMatchingAux(FirstOccurrence, LastColumn, Pattern):  
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
        top = FirstOccurrence[symbol] + Countsymbol(symbol, top, LastColumn) 
        # bottom <- FirstOccurrence(symbol) + Countsymbol(bottom + 1, LastColumn) - 1
        bottom = FirstOccurrence[symbol] + Countsymbol(symbol, bottom+1, LastColumn)-1 
      else :
        return 0
    else :
      return bottom - top + 1   

def rindex(lst, item):  
  rlst = list(reversed(lst))  
  return len(lst) - rlst.index(item) - 1
      
def BetterBWMatching(LastColumn, patterns):  
  FirstOccurrence = FirstOccurrenceColumn(LastColumn)  
  matches={}
  for pattern in patterns:
    matches[pattern] = BetterBWMatchingAux(FirstOccurrence, LastColumn, pattern)    
  return matches

def main_BetterBWMatching(myfile):    
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  
  with open(inputFile) as inFile:        
    transform = inFile.readline().strip()
    patterns = [x for x in inFile.readline().strip().split(' ')]
  
  matches = BetterBWMatching(transform, patterns)
    
  with open(outputFile, 'w') as outFile:        
    txt = ' '.join([str(matches[pattern]) for pattern in patterns])
    outFile.write(txt)  
    print txt    

main_BetterBWMatching('sample_BetterBWMatching')
main_BetterBWMatching('dataset_301_7')    
