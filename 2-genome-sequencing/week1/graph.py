def ReadKmers(inputFile):
  kmers = []
  with open(inputFile) as inFile:
    for line in inFile:
      line = line.strip()
      kmers.append(line)
  
  k = len(kmers[0])
  return k,kmers

def Composition(k,txt):
  kmers = []
  for i in range(len(txt)-k+1):
    kmers.append(txt[i:(i+k)])
  return kmers

def StringSpelledByGenomePath(kmers):
  txt = kmers[0]
  for i in range(1,len(kmers)):
    txt += kmers[i][-1]
  return txt

def main_composition(myfile):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  with open(inputFile,'r') as inFile:
    k = int(inFile.readline().strip())
    txt = inFile.readline().strip()
  kmers = Composition(k,txt)
  with open(outputFile,'w') as outFile:
    outFile.write('\n'.join(kmers))

def main_genomepath(myfile):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  k,kmers = ReadKmers(inputFile)
  txt = StringSpelledByGenomePath(kmers)
  with open(outputFile,'w') as outFile:
    outFile.write(txt)

'''
main_composition('sample_composition')
main_composition('dataset_197_3')
'''
main_genomepath('sample_genomepath')
main_genomepath('dataset_198_3')
