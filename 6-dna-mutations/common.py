def TextMatching(text,pattern):
  N = len(text)
  indices = []
  K = len(pattern)
  for i in range(0,N-K+1):
    #print ' '*i+text[i:(i+K)]
    if text[i:(i+K)] == pattern:
      indices.append(i)
  
  return indices