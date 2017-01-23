def TextMatching(text,patterns):
  N = len(text)
  indices = []
  for pattern in patterns:
    K = len(pattern)
    for i in range(0,N-K+1):
      #print ' '*i+text[i:(i+K)]
      if text[i:(i+K)] == pattern:
        indices.append(i)
  
  return indices