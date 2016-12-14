def ReadKmers(inputFile):
  kmers = []
  with open(inputFile) as inFile:
    for line in inFile:
      line = line.strip()
      kmers.append(line)
  
  k = len(kmers[0])
  return k,kmers

def Composition(k,txt):
  printout = len(txt)<25
  kmers = []
  for i in range(len(txt)-k+1):
    kmer = txt[i:(i+k)]
    if printout : print '{}{}'.format(' '*i,kmer)
    kmers.append(kmer)
  if printout : print txt
  return kmers

def StringSpelledByGenomePath(kmers):
  txt = kmers[0]
  for i in range(1,len(kmers)):
    txt += kmers[i][-1]
  return txt

def Fixes(kmers):
  prefix = []
  suffix = []
  for kmer in kmers:
    prefix.append(kmer[:-1])
    suffix.append(kmer[1:])
  return prefix,suffix

def Overlap(kmers):
  prefix,suffix = Fixes(kmers)
  adj = {}
  n = len(prefix)
  for i in range(n):
    for j in range(n):
      if i == j : continue
      if suffix[i] != prefix[j] : continue
      if i not in adj : 
        adj[i] = [j]
      else :
        adj[i].append(j)

  return adj 

def DeBruijnGraph(kmers):
  prefix,suffix = Fixes(kmers)  
  nNodes = len(kmers)
  nEdges = 0
  adj_dict = {}
  for i in range(nNodes):
    for j in range(nNodes):
      if suffix[i] != prefix[j] : continue       
      if i not in adj_dict:
        adj_dict[i] = []
      if j not in adj_dict[i]:
        adj_dict[i].append(j)
        nEdges += 1
        
  print '|(V,E)|=({},{})'.format(nNodes, nEdges)   
  return adj_dict
    
def AdjencyListToString(adj,kmers):
  txt = []
  for i in adj:
    from_node = kmers[i]
    to_nodes = ','.join([kmers[j] for j in adj[i]])
    txt.append('{} -> {}'.format(from_node,to_nodes))    
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

def main_overlap(myfile):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  k,kmers = ReadKmers(inputFile)
  adj = Overlap(kmers)
  txt = AdjencyListToString(adj,kmers)
  with open(outputFile,'w') as outFile:
    outFile.write('\n'.join(sorted(txt)))

def main_debruijn_fromstring(myfile):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  with open(inputFile,'r') as inFile:
    k = int(inFile.readline().strip())
    txt = inFile.readline().strip()
  kmers = Composition(k-1,txt)
  kmers = reduce(lambda l, x: l.append(x) or l if x not in l else l, kmers, []) # unique kmers
  adj = DeBruijnGraph(kmers)
  txt = AdjencyListToString(adj,kmers) 
  with open(outputFile,'w') as outFile:
    txt='\n'.join(sorted(txt))
    outFile.write(txt)

def main_debruijn_fromkmers(myfile):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  k,kmers = ReadKmers(inputFile)
  kmers = reduce(lambda l, x: l.append(x) or l if x not in l else l, kmers, []) # unique kmers
  adj = DeBruijnGraph(kmers)
  txt = AdjencyListToString(adj,kmers) 
  with open(outputFile,'w') as outFile:
    txt='\n'.join(sorted(txt))
    outFile.write(txt)

'''
main_composition('sample_composition')
main_composition('dataset_197_3')
main_genomepath('sample_genomepath')
main_genomepath('dataset_198_3')
main_overlap('sample_overlap')
main_overlap('dataset_198_10')
'''
main_debruijn_fromstring('sample_debruijn_fromstring')
main_debruijn_fromstring('dataset_199_6') # <--- needs fixing!
main_debruijn_fromkmers('sample_debruijn_fromkmers') # <--- needs fixing!
main_debruijn_fromkmers('dataset_200_8') # <--- needs fixing!
