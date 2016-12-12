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

def ReadKmers(inFile):
  kmers = []
  with open(inFile,'r') as inputFile:
    k = inputFile.readline().strip()    
    for line in inputFile:
      line = line.strip()
      kmers.append(line)
  return k,kmers

def DeBruijnGraph(k,kmers):
  nNodes = len(kmers)
  match = []
  adj_dict = {}

  for i in range(nNodes):
    for j in range(nNodes):
      if i == j: continue
      left_str = kmers[i]
      right_str = kmers[j]

      if left_str[1:] != right_str[:-1] : continue 
      
      left = kmers.index(left_str)
      right = kmers.index(right_str)

      if left not in adj_dict:
        adj_dict[left] = []
      if right not in adj_dict[left]:
        adj_dict[left].append(right)

  nEdges = 0
  for key in adj_dict:
    nEdges += len(adj_dict[key])

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
      
def main_Euler():
  #myfile = 'sample_cycle' #myfile = dataset_203_3
  #myfile = 'sample_path'  #myfile = dataset_203_6
  myfile = input("Enter filename (without file ending):")
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'

  adj,n = ReadAdjacencyList(inputFile)
  
  cycle, Eulerian = Euler(adj,n,n)
  
  with open(outputFile, 'w') as outFile:
    txt = PrintEulerGraph(cycle,Eulerian)
    outFile.write(txt)

def EulerGraphForKmers(k,kmers,isUnivStr):
  print 'Input {}-mers: {}'.format(k,kmers)
  adj,nNodes,nEdges = DeBruijnGraph(k,kmers)  
  
  if isUnivStr:
    for i in [0,1]:
      ix = kmers.index(str(i)*k)
      adj[ix].append(ix)
      nEdges+=1

  print '|(V,E)|=({},{})'.format(nNodes, nEdges)
  path,Eulerian = Euler(adj,nNodes,nEdges)  
  txt = PrintEulerPathKmer(path,k,kmers,Eulerian)
  return txt

def main_kmer():
  #myfile = 'sample_kmers'
  myfile = 'dataset_203_7'
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  
  k,kmers = ReadKmers(inputFile)
  txt = EulerGraphForKmers(k,kmers,False)
  
  with open(outputFile, 'w') as outFile:
    outFile.write(txt)
    print 'Solution:{}'.format(txt)

def main_univ():
  
  myfile = 'dataset_203_11'
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  
  with open(inputFile,'r') as inFile:
    k = int(inFile.readline().strip())

  kmers = BinaryStrings(k-1)
  txt = EulerGraphForKmers(k-1,kmers,True)  
  print '{}-universal string: {}'.format(k,txt)

  with open(outputFile, 'w') as outFile:
    outFile.write(txt)


main_univ()
