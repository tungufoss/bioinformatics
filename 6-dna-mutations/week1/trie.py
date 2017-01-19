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
  for child in children :     
    to_node, symbol = child
    nSubtrees, symbols = InOrder(to_node, edges, cedges, connected)
    
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
  
'''
main_TrieConstructionProblem('sample_TrieConstructionProblem')
main_TrieConstructionProblem('dataset_294_4')
main_TrieMatchingProblem('sample_TrieMatchingProblem')
main_TrieMatchingProblem('dataset_294_8')
main_SuffixTreeConstructionProblem('sample_SuffixTreeConstructionProblem')
main_SuffixTreeConstructionProblem('sample_SuffixTreeConstructionProblem2')
main_SuffixTreeConstructionProblem('SuffixTreeConstruction')
main_SuffixTreeConstructionProblem('dataset_296_4')
main_LongestRepeatProblem('sample_LongestRepeatProblem')
main_LongestRepeatProblem('dataset_296_5')
main_LongestSharedSubstringProblem('sample_LongestSharedSubstringProblem')
main_LongestSharedSubstringProblem('LongestSharedSubstring')
main_LongestSharedSubstringProblem('dataset_296_6')
'''
main_ShortestNonSharedSubstringProblem('sample_ShortestNonSharedSubstring')
main_ShortestNonSharedSubstringProblem('ShortestNonSharedSubstring')
main_ShortestNonSharedSubstringProblem('dataset_296_7')

