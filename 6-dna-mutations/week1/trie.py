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

def main_TrieConstructionProblem(myfile):  
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  
  with open(inputFile) as inFile:    
    patterns = [line.strip() for line in inFile]

  edges = TrieConstructionProblem(patterns)
  
  with open(outputFile, 'w') as outFile:
    for node in sorted(edges):
      line = '{}->{}:{}'.format(node[0],edges[node],node[1])
      print line
      outFile.write(line+'\n')

main_TrieConstructionProblem('sample_TrieConstructionProblem')
main_TrieConstructionProblem('dataset_294_4')