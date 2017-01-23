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
  print text
  # Next - similarly to suffix arrays - order all the cyclic rotations of Text lexicographically to form a |Text| x |Text| matrix of symbols that we call the Burrows-Wheeler matrix and denote by M(Text). 
  M = sorted(cyclic_rotations)  
  # We are interested in the last column of M(Text), called the Burrows-Wheeler transform of Text, or BWT(Text)
  BWT = ''.join([m[-1] for m in M])  
  return BWT

def InvBurrowsWheelerTransformConstruction(transform):
  transform = list(transform)
  N = len(transform)
  
  # add subscripts
  for i in reversed(range(N)):    
    transform[i] += str(transform[:i+1].count(transform[i]))
  
  FirstColumn = sorted(transform)
  LastColumn = transform    
  
  M = {}
  for i in range(N):
    M[i] = [FirstColumn[i]]+['?']*(N-2)+[LastColumn[i]]
  
  # Using the First-Last Property for inverting the Burrows-Wheeler transform    
  for row in range(N):    
    for col in range(1,N-1):    
      char = M[row][col-1]     
      nrow = LastColumn.index(char) 
      M[row][col] = M[nrow][0] # First column
    
  # for m in sorted(M): print m, [x[0] for x in M[m] ] # remove subscripts
  
  row = LastColumn.index('$1')
  text= ''.join([ M[row][i][0] for i in range(N) ])  
  return text
  
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

'''  
main_SuffixArrayConstruction('sample_SuffixArrayConstruction')
main_SuffixArrayConstruction('dataset_310_2')
main_BurrowsWheelerTransformConstruction('sample_BurrowsWheelerTransformConstruction')
main_BurrowsWheelerTransformConstruction('dataset_297_5')
'''
main_InverseBurrowsWheelerTransform('sample_InverseBurrowsWheelerTransformConstruction')
main_InverseBurrowsWheelerTransform('exbreak_InverseBurrowsWheelerTransformConstruction')
main_InverseBurrowsWheelerTransform('dataset_299_10')

