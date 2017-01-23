import sys
sys.path.insert(0,'../')
from common import TextMatching

def SuffixArrayConstruction(text):
  if text[-1] != '$':
    text+='$'
  
  suffixes = [text[i:] for i in range(len(text))]
  starting_pos=[]
  N=len(text)
  for s in sorted(suffixes):
    starting_pos.append(N-len(s))
  return starting_pos

def BurrowsWheelerTransformConstruction(text):  
  N=len(text)
  cyclic_rotations=[]  
  # First, form all possible cyclic rotations of Text;   
  for i in range(N):
    cyclic_rotations.append(text)
    # a cyclic rotation is defined by chopping off a suffix from the end of Text and appending this suffix to the beginning of Text. 
    text = text[-1]+text[:-1]
  
  # Next - similarly to suffix arrays - order all the cyclic rotations of Text lexicographically to form a |Text| x |Text| matrix of symbols that we call the Burrows-Wheeler matrix and denote by M(Text). 
  M = sorted(cyclic_rotations)  
  # We are interested in the last column of M(Text), called the Burrows-Wheeler transform of Text, or BWT(Text)
  BWT = ''.join([m[-1] for m in M])  
  return BWT

def InvBurrowsWheelerTransformConstruction(transform):
  N = len(transform)  
  # add subscripts  
  LastColumn = []
  for i in range(N):    
    subscript = transform[:i+1].count(transform[i])
    LastColumn.append((transform[i],subscript))
        
  FirstColumn = sorted(LastColumn)  
    
  row = LastColumn.index(('$',1))    
  if False : # full M matrix
    M = InvBurrowsWheelerMatrix(FirstColumn,LastColumn,N)    
    text= ''.join([ M[row][i][0] for i in range(N) ])    
  else : # only the row that is needed is computed 
    M = [FirstColumn[row]]+[('?',-1)]*(N-2)+[LastColumn[row]]
    M = FirstLastProperty(M, FirstColumn, LastColumn, N)
    text= ''.join([ m[0] for m in M ])    
    
  return text

def InvBurrowsWheelerMatrix(FirstColumn,LastColumn,N):
  M = {}
  for i in range(N):
    M[i] = [FirstColumn[i]]+[('?',-1)]*(N-2)+[LastColumn[i]]
  
  # Using the First-Last Property for inverting the Burrows-Wheeler transform    
  for row in range(N):    
    M[row]=FirstLastProperty(M[row],FirstColumn,LastColumn,N)

  #for m in sorted(M): print m, [x[0] for x in M[m] ] # remove subscripts
  return M
  
def FirstLastProperty(Mrow,FirstColumn,LastColumn,N):
  for col in range(1,N-1):    
    char = Mrow[col-1]     
    nrow = LastColumn.index(char) 
    Mrow[col] = FirstColumn[nrow]
  return Mrow

def BWMatching(text, patterns):
  matches = {}
  for pattern in patterns:
    matches[pattern] = len(TextMatching(text, pattern))

  return matches 
  
def main_SuffixArrayConstruction(myfile):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  
  with open(inputFile) as inFile:        
    text = inFile.readline().strip()
  
  suffix_arr = SuffixArrayConstruction(text)
  
  with open(outputFile, 'w') as outFile:    
    txt = ', '.join([str(x) for x in suffix_arr])
    outFile.write(txt)  
    print txt

def main_BurrowsWheelerTransformConstruction(myfile):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  
  with open(inputFile) as inFile:        
    text = inFile.readline().strip()
  
  BWT = BurrowsWheelerTransformConstruction(text)
  
  with open(outputFile, 'w') as outFile:        
    outFile.write(BWT)  
    print BWT  

def main_InverseBurrowsWheelerTransform(myfile):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  
  with open(inputFile) as inFile:        
    transform = inFile.readline().strip()
    
  text = InvBurrowsWheelerTransformConstruction(transform)
  print 'Sanity check: {}'.format(transform == BurrowsWheelerTransformConstruction(text))
  
  with open(outputFile, 'w') as outFile:        
    outFile.write(text)  
    print text

def main_BWMatching(myfile):    
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  
  with open(inputFile) as inFile:        
    transform = inFile.readline().strip()
    patterns = [x for x in inFile.readline().strip().split(' ')]
  
  text = InvBurrowsWheelerTransformConstruction(transform)
  matches = BWMatching(text, patterns)
  
  with open(outputFile, 'w') as outFile:        
    txt = ' '.join([str(matches[pattern]) for pattern in patterns])
    outFile.write(txt)  
    print txt    
    
'''  
main_SuffixArrayConstruction('sample_SuffixArrayConstruction')
main_SuffixArrayConstruction('dataset_310_2')
main_BurrowsWheelerTransformConstruction('sample_BurrowsWheelerTransformConstruction')
main_BurrowsWheelerTransformConstruction('dataset_297_5')
main_InverseBurrowsWheelerTransform('sample_InverseBurrowsWheelerTransformConstruction')
main_InverseBurrowsWheelerTransform('exbreak_InverseBurrowsWheelerTransformConstruction')
main_InverseBurrowsWheelerTransform('dataset_299_10')
'''
main_BWMatching('sample_BWMatching')
main_BWMatching('dataset_300_8')
