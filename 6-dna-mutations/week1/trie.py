import sys
sys.path.insert(0,'../')
from common import TextMatching

Inf = -1

def TrieConstructionProblem(patterns):
  edges = {}
  nNodes = 0
  
  for pattern in patterns:
    currentNode = 0 # start at root
    for symbol in pattern:
      if (currentNode, symbol) in edges: # match prev. pattern
        currentNode = edges[(currentNode, symbol)] # traverse to subtree
      else: # no prev match, add node as substree
        nNodes+=1
        edges[(currentNode,symbol)] = nNodes
        currentNode = nNodes

  return edges

def TrieMatching(text,patterns):
  N = len(text)
  indices = []
  for pattern in patterns:
    indices += TextMatching(text, pattern)
    
  return indices

def ConsolidateEdges(edges):
  cedges = {}
  connected = {}
  for edge in edges:
    from_edge = edge[0]
    to_edge = edges[edge]
    y = (to_edge, edge[1])
    if from_edge not in connected :
      connected[from_edge]=[y]
    else :
      connected[from_edge].append(y)
      
  InOrderTrie(0, edges, cedges, connected)
  return cedges

def InOrderTrie(from_node, edges, cedges, connected):   
  if from_node not in connected :
    return 0, ''  
    
  children = connected[from_node]  
  for child in children :     
    to_node, symbol = child
    nSubtrees, symbols = InOrderTrie(to_node, edges, cedges, connected)
    
    if nSubtrees <= 1:      
      symbols = symbol+symbols
    else :
      symbols = symbol
      
    if len(children)  > 1 :
      cedges[(from_node,symbols)]=to_node      
    
  return len(children), symbols
  
def SuffixTrie(text):
  if text[len(text)-1] != '$':
    text += '$' # mark the end of text   
  patterns = [text[k:] for k in range(len(text))]
  edges = TrieConstructionProblem(patterns)    
  edges = ConsolidateEdges(edges)  
  return edges

def LongestRepeat(text):  
  for l in reversed(range(1,len(text))):
    cnt = AuxRepeats(text,l)
    longest = max([cnt[c] for c in cnt])
    repeats = [c for c in cnt if cnt[c]==longest]    
    if longest > 1 :
      return repeats    
  return ['']

def AuxRepeats(text,length):
  from collections import Counter  
  patterns = [text[i:i+length] for i in range(len(text)-length+1)]    
  cnt = Counter(patterns) # in descending order
  return cnt
  
def LongestSharedRepeat(text1,text2):  
  for l in reversed(range(1,len(text1))):
    cnt1=AuxRepeats(text1,l)
    cnt2=AuxRepeats(text2,l)
    cnt = {}
    for c in [r for r in cnt1 if r in cnt2] :
      cnt[c] = cnt1[c]+cnt2[c]    
    if cnt != {} :      
      longest = max([cnt[c] for c in cnt])
      repeats = [c for c in cnt if cnt[c]==longest]
      return sorted(repeats)
  return ['']

def ShortestNonSharedSubstringProblem(text1,text2):
  print text1, text2
  for l in range(1,len(text1)):
    cnt1=AuxRepeats(text1,l)
    cnt2=AuxRepeats(text2,l)
    non_shared = [r for r in cnt1 if r not in cnt2] 
    if len(non_shared)>0 :      
      return non_shared
  return ['']
  
def PrintEdges(edges):
  txt=[]
  for node in sorted(edges):
    txt.append('{}->{}:{}'.format(node[0],edges[node],node[1]))
  return '\n'.join(txt)

def AdjencyStr2List(adjacency_str):
  adjacency_list={}
  for txt in adjacency_str :    
    from_node, to_nodes = txt.split(' -> ')
    if to_nodes == '{}' :
      adjacency_list[int(from_node)] = []
    else :
      adjacency_list[int(from_node)] = [int(x) for x in to_nodes.split(',')]
  return adjacency_list

def Degree(adjacency_list):
  indegree = {}
  outdegree = {}
  for from_node in adjacency_list:
    outdegree[from_node] = len(adjacency_list[from_node])
    for to_node in adjacency_list[from_node]:
      if to_node not in indegree:
        indegree[to_node] = 1
      else :
        indegree[to_node] += 1
  for node in [n for n in indegree if n not in outdegree]:
    outdegree[node] = 0
  for node in [n for n in outdegree if n not in indegree]:
    indegree[node] = 0
  return indegree, outdegree
  
def MaximalNonBranchingPaths(adjacency_list):
  import copy
  unused = copy.deepcopy(adjacency_list)  
  
  paths =[] # Paths <- empty list
  indegree, outdegree = Degree(adjacency_list)
  
  # for each node v in Graph
  for v in adjacency_list : 
    # if v is not a 1-in-1-out node
    if indegree[v] != 1 or outdegree[v] != 1 : 
      if outdegree[v]>0 :
        # for each outgoing edge (v, w) from v
        for w in adjacency_list[v] :
          # NonBranchingPath <- the path consisting of the single edge (v, w)
          nonBranchingPath = [v,w]
          unused[v].remove(w)
          # while w is a 1-in-1-out node
          while indegree[w] == 1 and outdegree[w] == 1 :
            u = adjacency_list[w][0]
            #extend NonBranchingPath by the outgoing edge (w, u) from w 
            nonBranchingPath.append(u)
            unused[w].remove(u)
            w = u # w <- u
          # add NonBranchingPath to the set Paths
          paths.append(nonBranchingPath)
  
  #for each isolated cycle Cycle in Graph  
  for v in adjacency_list :      
    if unused[v] == [] : continue
  
    w = adjacency_list[v][0]
    cycle = [v,w]
    unused[v].remove(w)
    
    # add Cycle to Paths    
    while v != w :      
      u = adjacency_list[w][0]
      cycle.append(u)
      unused[w].remove(u)
      w = u      
    paths.append(cycle)
  
  return paths

# TreeColoring: colors the nodes of a suffix tree from the leaves upward. 
def TreeColoring(adjacency_list, color_lbl):
  # This algorithm assumes that the leaves of the suffix tree have been labeled "blue" or "red" 
  # and all other nodes have been labeled "gray". 
  for node in adjacency_list:
    if node not in color_lbl:
      color_lbl[node] = 'gray'

  # A node in a tree is called ripe if it is gray but has no gray children.
  ripe_nodes = [v for v in color_lbl if color_lbl[v] == 'gray' and 'gray' not in [color_lbl[w] for w in adjacency_list[v]]]  
  # while ColoredTree has ripe nodes
  while len(ripe_nodes)>0:
    # for each ripe node v in ColoredTree
    for v in ripe_nodes :
      children = adjacency_list[v]
      children_color = list(set([color_lbl[u] for u in children]))      
      #if there exist differently colored children of v
      if len(children_color) > 1 : 
        color_lbl[v] = 'purple'
      else : 
        color_lbl[v] = children_color[0] # color of all children of v        
    ripe_nodes = [v for v in color_lbl if color_lbl[v] == 'gray']
  
def main_TrieConstructionProblem(myfile):  
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  
  with open(inputFile) as inFile:    
    patterns = [line.strip() for line in inFile]

  edges = TrieConstructionProblem(patterns)
  
  with open(outputFile, 'w') as outFile:
    txt = PrintEdges(edges)      
    outFile.write(txt)

def main_TrieMatchingProblem(myfile):      
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  
  with open(inputFile) as inFile:    
    text = inFile.readline().strip()
    patterns = [line.strip() for line in inFile]

  indices = TrieMatching(text,patterns)
    
  with open(outputFile, 'w') as outFile:
    txt = ' '.join([str(x) for x in sorted(indices)])    
    outFile.write(txt)
    
def main_SuffixTreeConstructionProblem(myfile):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  
  with open(inputFile) as inFile:    
    text = inFile.readline().strip()  
  print text
  edges = SuffixTrie(text)
    
  # output: The edge labels of SuffixTree(Text)
  with open(outputFile, 'w') as outFile:
    txt = '\n'.join(edge[1] for edge in edges)    
    outFile.write(txt)
    
def main_LongestRepeatProblem(myfile):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  
  with open(inputFile) as inFile:    
    text = inFile.readline().strip()  
  
  long = LongestRepeat(text)
  
  with open(outputFile, 'w') as outFile:
    outFile.write(long[0])
    print 'Text {} has solution(s): {}'.format(text,long)

def main_LongestSharedSubstringProblem(myfile):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  
  with open(inputFile) as inFile:    
    text1 = inFile.readline().strip()
    text2 = inFile.readline().strip()
  
  long = LongestSharedRepeat(text1,text2)
  
  with open(outputFile, 'w') as outFile:
    outFile.write(long[0])
    print 'Texts ({},{}) has shared solution(s): {}'.format(text1,text2,long)

def main_ShortestNonSharedSubstringProblem(myfile):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  
  with open(inputFile) as inFile:    
    text1 = inFile.readline().strip()
    text2 = inFile.readline().strip()
  
  short = ShortestNonSharedSubstringProblem(text1,text2)  
  
  with open(outputFile, 'w') as outFile:
    outFile.write(short[0])
    print 'Texts ({},{}) has shared solution(s): {}'.format(text1,text2,short)

def main_MaximalNonBranchingPaths(myfile):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  
  with open(inputFile) as inFile:        
    adjacency_str = [line.strip() for line in inFile]

  adjacency_list = AdjencyStr2List(adjacency_str)  
  paths = MaximalNonBranchingPaths(adjacency_list)

  with open(outputFile, 'w') as outFile:
    for path in paths :
      txt = ' -> '.join([str(x) for x in path])
      outFile.write(txt+'\n')        

def main_TreeColoring(myfile):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  
  with open(inputFile) as inFile:
    input = [line.strip() for line in inFile]
    ix = input.index('-')    
    adjacency_list = AdjencyStr2List(input[:ix])
    color_lbl={}
    for x in input[ix+1:]:
      node,color = x.split(': ')
      color_lbl[int(node)]=color
  
  TreeColoring(adjacency_list, color_lbl)
  
  with open(outputFile, 'w') as outFile:
    for node in sorted(color_lbl):
      txt = '{}: {}'.format(node, color_lbl[node])
      outFile.write(txt+'\n')  

def main_test():

  print '#1\tThe most efficient way to map reads to a reference genome is by using global alignment to align each read to the reference genome.'
  print '\tFalse'
  
  print '#2\tWhat is the runtime of BruteForcePatternMatching?','n = |text| and m = |pattern|'
  print '\tWorst case is when you have m=1, O(nm)'
  print '\tNote: Best case is when you have m=n, O(m)'
  
  print '#3\tWhich patterns could have been used to construct the following trie? (Select all that apply.)'
  print '\t[{}]\tdale'.format(True)
  print '\t[{}]\tbandana'.format(True)
  print '\t[{}]\tbanana'.format(False)
  print '\t[{}]\tdanny'.format(False)
  
  pattern = 'na'
  print '#4\tGiven the following suffix trie/suffix tree of an unknown text, identify all starting positions in the underlying text of\n\tPattern = "{}". Use 0-based indexing and enter the starting positions as space-separated integers in increasing order.'.format(pattern)
  text = 'adnaadnanadadnan$'
  indices = TextMatching(text,pattern)
  print '\tText: "{}" (|Text|={}) has starting positions {} for "{}"'.format(text, len(text[:-1]), indices, pattern)
  
  text='AAATCCCCCCCTCGCGACTTCTAGAGAAGAAGAGTACT$'
  print '#5\tHow many leaves will SuffixTree("{}") have?'.format(text)
  edges = SuffixTrie(text)
  leafs = [x for x in edges if x[1][-1]=='$']
  print '\t#{} leafs'.format(len(leafs))

main_TrieConstructionProblem('sample_TrieConstructionProblem')
main_TrieConstructionProblem('dataset_294_4')
main_TrieMatchingProblem('sample_TrieMatchingProblem')
main_TrieMatchingProblem('dataset_294_8')
main_SuffixTreeConstructionProblem('sample_SuffixTreeConstructionProblem')
main_SuffixTreeConstructionProblem('sample_SuffixTreeConstructionProblem2')
main_SuffixTreeConstructionProblem('dataset_296_4')
main_LongestRepeatProblem('sample_LongestRepeatProblem')
main_LongestRepeatProblem('dataset_296_5')
main_LongestSharedSubstringProblem('sample_LongestSharedSubstringProblem')
main_LongestSharedSubstringProblem('LongestSharedSubstring')
main_LongestSharedSubstringProblem('dataset_296_6')
main_ShortestNonSharedSubstringProblem('sample_ShortestNonSharedSubstring')
main_ShortestNonSharedSubstringProblem('ShortestNonSharedSubstring')
main_ShortestNonSharedSubstringProblem('dataset_296_7')
main_MaximalNonBranchingPaths('sample_MaximalNonBranchingPaths')
main_MaximalNonBranchingPaths('dataset_6207_2')
main_TreeColoring('sample_TreeColoring')
main_TreeColoring('dataset_9665_6')
main_test()