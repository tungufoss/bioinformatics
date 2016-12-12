'''
EulerianCycle(Graph)
  form a cycle Cycle by randomly walkin gin Graph (don't visit the same edge twice!) 
  while there are unexplored edges in Graph 
    select a node newStart in Cycle with still unexplored edges
    form Cycle' by traversing Cycle (starting at newStart) and then randomly walking
    Cycle <- Cycle'
  return Cycle

'''

def BinaryStrings(k):
  import itertools
  return ["".join(seq) for seq in itertools.product("01", repeat=k)]

def ReadKmers(inputFile):
  kmers = []
  with open(inputFile,'r') as inFile:
    k = inFile.readline().strip()    
    for line in inFile:
      line = line.strip()
      kmers.append(line)
  return k,kmers

def ReadKmerReads(inputFile):
  paired = []
  read1 = []
  read2 = []
  with open(inputFile,'r') as inFile:
    line = inFile.readline().strip().split(' ')
    k = int(line[0])
    d = int(line[1])
    for line in inFile:
      line = line.strip()
      paired.append(line)
      line = line.split('|')
      read1.append(line[0])
      read2.append(line[1])

  return k,d,paired,read1,read2

def KDmerComposition(k,d,txt):
  n = len(txt)
  print '({},{})-mer composition for |{}|={}'.format(k,d,txt,n)
  reads = []
  read1 = []
  read2 = []
  for i in range(n-k*2-d+1):
    read1.append(txt[i:i+k])
    read2.append(txt[i+k+d:i+k+d+k])
    reads.append('{}|{}'.format(read1[-1],read2[-1]))
  
  return reads,read1,read2

def DeBruijnGraphAuxFun(kmers,suffix,prefix):
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
  
  return adj_dict,nNodes,nEdges

def DeBruijnGraph(kmers):
  prefix = []
  suffix = []
  for i in range(len(kmers)):
    prefix.append(kmers[i][:-1])
    suffix.append(kmers[i][1:])

  return DeBruijnGraphAuxFun(kmers,suffix,prefix)

def PairedDeBruijnGraph(paired,read1,read2):
  prefix = []
  suffix = []
  for i in range(len(paired)):
    prefix.append('{}|{}'.format(read1[i][:-1],read2[i][:-1]))
    suffix.append('{}|{}'.format(read1[i][1:],read2[i][1:]))
  
  adj_dict,nNodes,nEdges = DeBruijnGraphAuxFun(paired,suffix,prefix)
  

  return adj_dict,nNodes,nEdges

def ReadAdjacencyList(inFile):

  adj_dict = {}  
  n = 0

  with open(inFile, 'r') as inputFile:
    for line in inputFile:
      line = line.strip()
      line = line.split(' -> ')
      from_node = line[0]
      to_nodes = line[1].split(',')
      adj_dict[int(from_node)] = [int(x) for x in to_nodes]
      n += len(to_nodes)

  return (adj_dict,n)

def rotate(l, n):
  return l[n:] + l[:n]

def EulerianCycle(adj,n):
  #print '#{} adj:{}'.format(n,adj)

  # inital start (could be random, let's use first node)
  cycle = []
  target = 0
  
  cycle.append(adj[0].pop()) # pop removes item from adj list

  while len(cycle) < n:

    #print "begin", cycle, target
    while cycle[-1] != target:
      cycle.append(adj[cycle[-1]].pop(0))

    if len(cycle) == n : break 

    # Find someone with a remaining edge
    while len(adj[cycle[-1]]) == 0:
      cycle = rotate(cycle, 1)
    cycle.append(adj[cycle[-1]].pop(0))

    target = cycle[-1]
    #print "end", cycle
    
  return cycle 

def EulerianPath(adj, nEdges, sink, source):
    print 'Adding dummy edge from sink to source: {}->{}'.format(sink,source)    
    if sink in adj:
      adj[sink].append(source)
    else: 
      adj[sink] = [source]

    nEdges += 1 # note we added dummy edge
   
    cycle = EulerianCycle(adj,nEdges)

    for i in range(nEdges):
      from_node = cycle[i]
      to_node = cycle[(i+1)%nEdges]
      if to_node != source or from_node != sink:
        continue

      path = rotate(cycle,(i+1)%nEdges)
      return path
    
def IsEulerian(adj,nNodes):

  outdegree = [0]*nNodes
  indegree = [0]*nNodes
  for from_node in adj:    
    outdegree[from_node]+=len(adj[from_node])
    for to_node in adj[from_node]:      
      indegree[to_node]+=1
  #print 'outdegree:{}'.format(outdegree)
  #print 'indegree: {}'.format(indegree)

  sink = -1
  source = -1
  
  for i in range(nNodes):
    if outdegree[i] > indegree[i] and source == -1:
      source = i
    elif outdegree[i] < indegree[i] and sink == -1:
      sink = i
    elif outdegree[i] != indegree[i]:
      return False, sink, source

  if (sink == -1 and source == -1) :
    return True, sink, source
  else :
    return False, sink, source

def PrintEulerGraph(cycle,Eulerian):
  txt = '{}'.format(cycle[0])
  for k in range(1,len(cycle)):
    txt += '->{}'.format(cycle[k])
  if Eulerian: 
    txt += '->{}'.format(str(cycle[0]))
  return txt

def PrintEulerPathKmer(path,k,kmers,Eulerian):
  txt = kmers[path[0]]
  for i in range(1,len(path)):
    txt += kmers[path[i]][-1]

  if Eulerian:
    txt = txt[k-1:]

  return txt
  
def Euler(adj,nNodes,nEdges):
  Eulerian, sink, source = IsEulerian(adj,nNodes)

  if Eulerian == True:
    cycle = EulerianCycle(adj,nEdges)
  else :
    cycle = EulerianPath(adj, nEdges, sink, source)

  return cycle, Eulerian
      
def main_Euler(myfile):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'

  adj,n = ReadAdjacencyList(inputFile)
  
  cycle, Eulerian = Euler(adj,n,n)
  
  with open(outputFile, 'w') as outFile:
    txt = PrintEulerGraph(cycle,Eulerian)
    outFile.write(txt)

def EulerGraphForKmers(k,kmers):
  print 'Input {}-mers: {}'.format(k,kmers)
  adj,nNodes,nEdges = DeBruijnGraph(kmers)    
  path,Eulerian = Euler(adj,nNodes,nEdges)  
  txt = PrintEulerPathKmer(path,k,kmers,Eulerian)
  return txt

def RemoveDuplicatePairs(paired,read1,read2):
  nNodes=len(paired)
  unique_pairs = []
  for i in reversed(range(nNodes)):
    if paired[i] not in unique_pairs:
      unique_pairs.append(paired[i])
    else:
      paired.pop(i)
      read1.pop(i)
      read2.pop(i)
  return paired,read1,read2

def EulerGraphForPairedKmers(k,d,paired,read1,read2):
  print 'Input ({},{})-mer composition'.format(k,d)

  paired,read1,read2=RemoveDuplicatePairs(paired,read1,read2)
  adj, nNodes, nEdges = PairedDeBruijnGraph(paired,read1,read2)  
  path, Eulerian = Euler(adj,nNodes,nEdges)

  txt1 = PrintEulerPathKmer(path,k,read1,Eulerian)
  txt2 = PrintEulerPathKmer(path,k,read2,Eulerian)
  txt = txt1[:k+d]+txt2
  return txt  

def main_kmer(myfile):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  
  k,kmers = ReadKmers(inputFile)
  txt = EulerGraphForKmers(k,kmers)
  
  with open(outputFile, 'w') as outFile:
    outFile.write(txt)
    print 'Solution:{}'.format(txt)

def main_univ(myfile):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  
  with open(inputFile,'r') as inFile:
    k = int(inFile.readline().strip())

  kmers = BinaryStrings(k-1)
  txt = EulerGraphForKmers(k-1,kmers)  
  print '{}-universal string: {}'.format(k,txt)
  print 'Number of such strings:',2**(2**(k-1)-k)

  with open(outputFile, 'w') as outFile:
    outFile.write(txt)

def main_reads(myfile):
  '''
  k=3
  d=2
  txt='TAATGCCATGGGATGTT'
  paired,read1,read2 = KDmerComposition(k,d,txt)
  #with open('dummy.out','w') as outFile: outFile.write(' '.join(sorted(paired)))
  print ' '.join(paired)
  print txt == EulerGraphForPairedKmers(k,d,paired,read1,read2)
  '''
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'  
  k,d,paired,read1,read2 = ReadKmerReads(inputFile)
  txt = EulerGraphForPairedKmers(k,d,paired,read1,read2)
  with open(outputFile, 'w') as outFile:
    outFile.write(txt)

'''
main_Euler('sample_cycle') #main_Euler('dataset_203_3')
main_Euler('sample_path') #main_Euler('dataset_203_6')
main_kmer('sample_kmers') #main_kmer('dataset_203_7')
main_univ('dataset_203_11') 
'''
main_reads('sample_reads')
main_reads('dataset_204_15')

