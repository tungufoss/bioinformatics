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
    
'''  
main_SuffixArrayConstruction('sample_SuffixArrayConstruction')
main_SuffixArrayConstruction('dataset_310_2')
'''
main_BurrowsWheelerTransformConstruction('sample_BurrowsWheelerTransformConstruction')
main_BurrowsWheelerTransformConstruction('dataset_297_5')
