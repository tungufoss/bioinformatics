import sys
sys.path.insert(0,'../')
from common import TextMatching
sys.path.insert(0,'../week2')
from BurrowsWheelerTransform import InvBurrowsWheelerTransformConstruction, BurrowsWheelerTransformConstruction

def BWTarrays(LastColumn):  
  CharSet = sorted(set(LastColumn))
  CharSet = {}
  FirstColumn = sorted(LastColumn)  
  for symbol in set(LastColumn):
    FO = FirstColumn.index(symbol)
    CO = FirstColumn.count(symbol)    
    CharSet[symbol] = (FO, CO)
      
  CntCharSet = []    
  for i in range(len(LastColumn)):
    CntCharSet.append([LastColumn[:i].count(symbol) for symbol in sorted(CharSet)])    
  CntCharSet.append([LastColumn.count(symbol) for symbol in sorted(CharSet)])
  
  if False: 
    N=len(LastColumn)
    print 'i\tFC\tLC\t(FO,CO)\t{}\n'.format(sorted(CharSet)),'-'*100
    for i in range(N):
      print '{}\t{}\t{}\t{}\t{}'.format(i, FirstColumn[i], LastColumn[i], CharSet[FirstColumn[i]], CntCharSet[i]) 
    print '{}\t\t\t\t{}'.format(N, CntCharSet[i]) 
    
  return FirstColumn, CharSet, CntCharSet

def Last2FirstColumn(LastColumn):
  N = len(LastColumn)    
  SL = [(LastColumn[i], LastColumn[:i+1].count(LastColumn[i])) for i in range(N) ]
  SF = sorted(SL)  
  LastToFirst = [SF.index(SL[i]) for i in range(N)]  
  return LastToFirst
  
def BetterBWMatchingAux(FirstColumn, LastColumn, Pattern, LastToFirst):  
  top = 0 # top <- 0
  bottom = len(LastColumn)-1 # bottom <- |LastColumn| - 1
  while top <= bottom : 
    if len(Pattern) > 0 : #Pattern is nonempty
      symbol = Pattern[-1] #symbol <- last letter in Pattern
      Pattern = Pattern[:-1] # remove last letter from Pattern      
      # if positions from top to bottom in LastColumn contain an occurrence of symbol
      check = LastColumn[top:bottom+1]
      if symbol in check:
        tix = top+check.index(symbol)   #topIndex <- first position where symbol occurs among positions where symbol occurs from top to bottom in LastColumn
        bix = top+rindex(check,symbol)  #bottomIndex <- last position of symbol among positions from top to bottom in LastColumn
        top = LastToFirst[tix]    #top <- LastToFirst(topIndex)
        bottom = LastToFirst[bix] #bottom <- LastToFirst(bottomIndex)
      else :
        return 0
    else :
      return bottom - top + 1   

def rindex(lst, item):  
  rlst = list(reversed(lst))  
  return len(lst) - rlst.index(item) - 1
      
def BetterBWMatching(LastColumn, patterns):
  FirstColumn = sorted(LastColumn)
  LastToFirst = Last2FirstColumn(LastColumn)
  matches={}
  for pattern in patterns:
    matches[pattern] = BetterBWMatchingAux(FirstColumn, LastColumn, pattern, LastToFirst)
    
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
