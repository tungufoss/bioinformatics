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

def ReadKmers(inputFile,startsWithK):
  kmers = []
  with open(inputFile,'r') as inFile:
    if startsWithK==True:
      k = int(inFile.readline().strip())
    for line in inFile:
      line = line.strip()
      kmers.append(line)
  
  if startsWithK==False:
    k = len(kmers[0])

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


def EulerianPath(adj, nEdges, nNodes, sink, source):
  nNeeded,indegree,outdegree = MinEdgesNeededForBalancing(adj,nNodes)
  if nNeeded != 1 : 
    print 'Not Eulerian, need at least {} more edges'.format(nNeeded)
    return []

  print 'Adding dummy edge from sink to source: {}->{}'.format(sink,source)    
  if sink in adj:
    adj[sink].append(source)
  else: 
    adj[sink] = [source]   
  nEdges+=1 # added a dummy edge

  cycle = EulerianCycle(adj,nEdges)

  for i in range(nEdges):
    from_node = cycle[i]
    to_node = cycle[(i+1)%nEdges]
    if to_node != source or from_node != sink:
      continue
    path = rotate(cycle,(i+1)%nEdges)

  return path
    
def DegreeOfNodes(adj,nNodes):
  outdegree = [0]*nNodes
  indegree = [0]*nNodes
  for from_node in adj:    
    outdegree[from_node]+=len(adj[from_node])
    for to_node in adj[from_node]:      
      indegree[to_node]+=1
  #print 'outdegree:{}'.format(outdegree)
  #print 'indegree: {}'.format(indegree)
  return indegree, outdegree

def Cyclic(adj,nNodes):
  indegree,outdegree = DegreeOfNodes(adj,nNodes)
  sink = -1
  source = -1  
  for i in range(nNodes):
    if outdegree[i] > indegree[i] and source == -1:
      source = i
    elif outdegree[i] < indegree[i] and sink == -1:
      sink = i
    
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
  cyclic, sink, source = Cyclic(adj,nNodes)
  if cyclic == True:
    return EulerianCycle(adj, nEdges), True
  else :
    return EulerianPath(adj, nEdges, nNodes, sink, source), False
      
def main_Euler(myfile):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'

  adj,n = ReadAdjacencyList(inputFile)
 
  path, cyclic = Euler(adj,n,n)
  
  with open(outputFile, 'w') as outFile:
    txt = PrintEulerGraph(path,cyclic)
    outFile.write(txt)

def EulerGraphForKmers(k,kmers):
  #print 'Input {}-mers: {}'.format(k,kmers)
  adj,nNodes,nEdges = DeBruijnGraph(kmers) 
  path,Eulerian = Euler(adj,nNodes,nEdges)  
  if len(path)==0 : return '<No path found>'
  txt = PrintEulerPathKmer(path,k,kmers,Eulerian)
  return txt

def Contig(k,kmers):
  #print 'Input {}-mers: {}'.format(k,kmers)
  adj,nNodes,nEdges = DeBruijnGraph(kmers)

  indegree,outdegree = DegreeOfNodes(adj,nNodes)
  for i in range(nNodes):
    while outdegree[i] == 1 and indegree[adj[i][0]] == 1:
      from_node = kmers[i]
      j = adj[i][0]
      to_node = kmers[j]
        
      kmer = from_node+to_node[k-1:]
      #print '{}->{}\nA: {}\nB:  {}\nAB:{}'.format(i,j,from_node,to_node,kmer) 
      kmers[i] = kmer
      if outdegree[j] > 0:
        adj[i] = adj[j]
        del adj[j]
        outdegree[i] = outdegree[j]
      else:          
        del adj[i]
        outdegree[i] = 0

      outdegree[j] = -99
      indegree[j] = -99
      nEdges-=1
        
  contig = []
  for i in range(nNodes):
    if indegree[i] > 0 or outdegree[i] > 0 : 
      contig.append(kmers[i])

  '''
  source = indegree.index(0)
  sink = outdegree.index(0)
  print 'Source:{} and Sink:{}'.format(source, sink)
  '''
  return contig,adj,kmers
  
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
  #print 'Input ({},{})-mer composition'.format(k,d)

  paired,read1,read2=RemoveDuplicatePairs(paired,read1,read2)
  adj, nNodes, nEdges = PairedDeBruijnGraph(paired,read1,read2)  
  path, Eulerian = Euler(adj,nNodes,nEdges)
  if len(path)==0 : return '<No path found>'

  txt1 = PrintEulerPathKmer(path,k,read1,Eulerian)
  txt2 = PrintEulerPathKmer(path,k,read2,Eulerian)
  txt = txt1[:k+d]+txt2
  return txt  

def main_kmer(myfile):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  
  k,kmers = ReadKmers(inputFile,True)
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

def main_contigs(myfile):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  k,kmers = ReadKmers(inputFile,False)  
  contigs,adj,kmers_merged = Contig(k,kmers[:])
  with open(outputFile, 'w') as outFile:
    txt = ' '.join(sorted(contigs))
    outFile.write(txt)

def MinEdgesNeededForBalancing(adj,nNodes):  
  indegree,outdegree = DegreeOfNodes(adj,nNodes)
  unbalanced_in = 0
  unbalanced_out = 0
  for i in range(nNodes):
    if indegree[i] < outdegree[i] : 
      unbalanced_in += outdegree[i]-indegree[i]
    elif indegree[i] > outdegree[i] : 
      unbalanced_out += indegree[i]-outdegree[i]
  return max(unbalanced_in,unbalanced_out),indegree,outdegree

def main_test():
  print '\n\nQ1 Give a linear string having the following 4-mer composition.'
  k,kmers = ReadKmers('test_q1_kmers.txt',True)
  contig, adj, kmers_merged = Contig(k, kmers[:])
  print adj
  for key in adj:
    print '{} {}'.format(key,kmers_merged[key]) 
  print 'Q1 sol: <read from solution>'
  
  print '\n\nQ2 Below is the adjacency list of a graph. What is the minimum number of edges we must add to this graph in order to make each node balanced? (You may add duplicate edges connecting the same two nodes, but do not add new nodes.)'
  adj,nNodes = ReadAdjacencyList('test_q2_adjlist.txt')
  nNeeded,indegree,outdegree = MinEdgesNeededForBalancing(adj,nNodes)
  print 'Q2 sol: Need at least {} edges to balance'.format(nNeeded) 
      
  print '\n\nQ3 There is a single (linear) string with the following (3,1)-mer composition. Find it.'
  k,d,paired,read1,read2 = ReadKmerReads('test_q3_reads.txt')
  txt = EulerGraphForPairedKmers(k,d,paired,read1,read2)
  print 'Q3 sol:',txt

def main_univstr(n):  
  kmers = BinaryStrings(n-1)
  txt = EulerGraphForKmers(n-1,kmers)
  txt += txt[0:n-1]
  print '{}-universal string: {}'.format(n,txt)
  N = 2**(2**(n-1)-n)
  print 'Number of distinct De Bruijn sequences B(2,{}): {}'.format(n,N)
  print 'Number of such binary strings:',N**2

'''
main_Euler('sample_cycle') #main_Euler('dataset_203_3')
main_Euler('sample_path') #main_Euler('dataset_203_6')
main_kmer('sample_kmers') #main_kmer('dataset_203_7')
#main_univ('dataset_203_11') 
main_reads('sample_reads') #main_reads('dataset_204_15')
main_contigs('sample_contigs') #main_contigs('dataset_205_5')
main_reads('sample_gappedpatterns') #main_reads('dataset_6206_7')
'''
#main_test()
main_univstr(4)
