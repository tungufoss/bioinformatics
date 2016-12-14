def Composition(k,txt):
  kmers = []
  for i in range(len(txt)-k+1):
    kmers.append(txt[i:(i+k)])
  return kmers

def main_composition(myfile):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  with open(inputFile,'r') as inFile:
    k = int(inFile.readline().strip())
    txt = inFile.readline().strip()
  kmers = Composition(k,txt)
  with open(outputFile,'w') as outFile:
    outFile.write('\n'.join(kmers))

main_composition('sample_composition')
main_composition('dataset_197_3')
