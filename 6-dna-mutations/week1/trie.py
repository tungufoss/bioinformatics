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
  for (from_node,symbol,to_node) in [(x[0],x[1],edges[x]) for x in edges if x[1]=='$']:        
    symbols = symbol
    while True:
      subtrees = [(y[0],y[1],edges[y]) for y in edges if y[0]==from_node]
      connected = [(y[0],y[1],edges[y]) for y in edges if edges[y]==from_node]      
      del edges[(from_node,symbol)]      
      if len(subtrees) == 1 :
        if (from_node,symbol) in edges :
          print 'wtf'      
        from_node, symbol, to_node = connected[0]
        symbols = symbol+symbols                
      else :
        edges[(from_node,symbols)]=to_node        
        break    
  return edges
  
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
    for edge in edges:
      outFile.write(edge[1])
'''
main_TrieConstructionProblem('sample_TrieConstructionProblem')
main_TrieConstructionProblem('dataset_294_4')
main_TrieMatchingProblem('sample_TrieMatchingProblem')
main_TrieMatchingProblem('dataset_294_8')
''' 
main_SuffixTreeConstructionProblem('sample_SuffixTreeConstructionProblem')
main_SuffixTreeConstructionProblem('dataset_296_4')