def SuffixArrayConstruction(text):
  if text[-1] != '$':
    text+='$'
  
  suffixes = [text[i:] for i in range(len(text))]
  starting_pos=[]
  N=len(text)
  for s in sorted(suffixes):
    starting_pos.append(N-len(s))
  return starting_pos
  
  for s in suffixes:
    print s 
  return []

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
  
main_SuffixArrayConstruction('sample_SuffixArrayConstruction')
main_SuffixArrayConstruction('dataset_310_2')