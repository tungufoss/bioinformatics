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
    K = len(pattern)    
    for i in range(0,N-K+1):
      #print ' '*i+text[i:(i+K)]
      if text[i:(i+K)] == pattern:
        indices.append(i)
  
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
      
  InOrder(0, edges, cedges, connected)
  return cedges

def InOrder(from_node, edges, cedges, connected):  
  if from_node not in connected :
    return 0, ''  
  children = connected[from_node]
  nChildren = len(children)
  for child in children :     
    to_node, symbol = child
    nSubtrees, symbols = InOrder(to_node, edges, cedges, connected)
    if nSubtrees <= 1:      
      symbols = symbol+symbols    
    if nChildren>1 :
      if nSubtrees > 1 :
        cedges[(from_node,symbol)]=to_node
      else :
        cedges[(from_node,symbols)]=to_node      
  return len(children), symbols
  
def SuffixTrie(text):
  if text[len(text)-1] != '$':
    text += '$' # mark the end of text   
  patterns = [text[k:] for k in range(len(text))]      
  edges = TrieConstructionProblem(patterns)  
  edges = ConsolidateEdges(edges)
  return edges

def PrintEdges(edges):
  txt=[]
  for node in sorted(edges):
    txt.append('{}->{}:{}'.format(node[0],edges[node],node[1]))
  return '\n'.join(txt)
 
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
'''
main_TrieConstructionProblem('sample_TrieConstructionProblem')
main_TrieConstructionProblem('dataset_294_4')
main_TrieMatchingProblem('sample_TrieMatchingProblem')
main_TrieMatchingProblem('dataset_294_8')
''' 
main_SuffixTreeConstructionProblem('sample_SuffixTreeConstructionProblem')
main_SuffixTreeConstructionProblem('dataset_296_4')