def ProbabilityHiddenPathProblem(Sigma, Transition):
  p=0.5 #  an assumption that is modeled by setting Pr(pi0 -> pi1) = 1/2, where pi0 is a "silent" initial state that does not emit any symbols.
  for i in range(1,len(Sigma)):
    p*=Transition[(Sigma[i-1],Sigma[i])]  
  return p

def ProbabilityOutcomeGivenHPP(x, pi, Emission):
  p = 1.0  
  for i in range(len(x)):
    p*=Emission[(pi[i],x[i])] 
  return p

def ViterbiAlgorithm(string, states, transition, emission):    
  dic = {}  
  nStates = len(states)
  nSeq = len(string)
  
  dic[0] = {}
  initprob = 1.0/len(states)
  for i in states:
    dic[0][i] = {}
    dic[0][i]['p'] = initprob*emission[(i,string[0])]
    
  for i in xrange(1,nSeq):
    dic[i] = {}
    for k in states:
      values = [dic[i-1][l]['p']*transition[(l, k)] for l in states]
      maxval = max(values)
      maxvar = states[values.index(maxval)]
      dic[i][k] = {}
      dic[i][k]['p'] = maxval*emission[(k,string[i])]
      dic[i][k]['pre_state'] = maxvar      
  
  maxval = max([dic[nSeq-1][l]['p'] for l in states])  
  maxvar = [l for l in states if dic[nSeq-1][l]['p']==maxval][0]
  pi = [maxvar]
  
  for i in reversed(range(1,nSeq)):
    maxvar = dic[i][maxvar]['pre_state']
    pi.append(maxvar)
  return reversed(pi)
   
def ReadMatrix(inFile, nRows):
  colnames = inFile.readline().strip().split('\t')
  nCol = len(colnames)
  mat = {}
  for ii in range(nRows):       
    line = inFile.readline()    
    line = line.strip().split('\t')    
    for i in range(nCol):
      mat[(line[0],colnames[i])]=float(line[i+1])
  return mat

def PrintMatrix(mat, rownames, colnames, name='\t'):  
  print '{}:\t{}'.format(name, '\t'.join(colnames))
  for s in rownames:
    print '\t{}\t{}'.format(s, '\t'.join([str(mat[(s,t)]) for t in colnames]))
    
def main_ProbabilityHiddenPathProblem(myfile):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  
  with open(inputFile) as inFile:        
    Sigma = inFile.readline().strip()
    inFile.readline() # --------
    States = inFile.readline().strip().split(' ')
    inFile.readline() # --------
    Transition = ReadMatrix(inFile, len(States))
  
  prob = ProbabilityHiddenPathProblem(Sigma, Transition)
  with open(outputFile, 'w') as outFile:    
    outFile.write(str(prob))
    print 'Probability of a Hidden Path Problem for {}: {}'.format(Sigma, prob)
    
def main_ProbabilityOutcomeGivenHPP(myfile):

  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  
  with open(inputFile) as inFile:        
    x = inFile.readline().strip()
    inFile.readline() # --------
    x_states = inFile.readline().strip().split(' ')
    inFile.readline() # --------
    pi = inFile.readline().strip()
    inFile.readline() # --------
    pi_states = inFile.readline().strip().split(' ')
    inFile.readline() # --------
    Emission = ReadMatrix(inFile, len(pi_states))
  
  prob = ProbabilityOutcomeGivenHPP(x, pi, Emission)
  with open(outputFile, 'w') as outFile:    
    outFile.write(str(prob))
    print 'Probability that an HMM will emit a {} given its hidden path {}: {}'.format(x, pi, prob)

def main_ViterbiAlgorithm(myfile):    
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  
  with open(inputFile) as inFile:        
    x = inFile.readline().strip()
    inFile.readline() # --------
    x_states = inFile.readline().strip().split(' ')
    inFile.readline() # --------
    pi_states = inFile.readline().strip().split(' ')
    inFile.readline() # --------
    Transition = ReadMatrix(inFile,len(pi_states))    
    inFile.readline() # --------
    Emission = ReadMatrix(inFile,len(pi_states))

  PrintMatrix(Transition, pi_states, pi_states, 'Transition')
  PrintMatrix(Emission, pi_states, x_states, 'Emission')
  
  pi = ViterbiAlgorithm(x, pi_states, Transition, Emission)    
  HMM = ''.join(pi)
      
  with open(outputFile, 'w') as outFile:    
    outFile.write(HMM)
    print HMM
    
'''
main_ProbabilityHiddenPathProblem('sample_ProbabilityHiddenPathProblem')
main_ProbabilityHiddenPathProblem('dataset_11594_2')
main_ProbabilityOutcomeGivenHPP('sample_ProbabilityOutcomeGivenHPP')
main_ProbabilityOutcomeGivenHPP('extra_ProbabilityOutcomeGivenHPP')
main_ProbabilityOutcomeGivenHPP('dataset_11594_4')
'''
main_ViterbiAlgorithm('sample_ViterbiAlgorithm')
main_ViterbiAlgorithm('extra_ViterbiAlgorithm')
main_ViterbiAlgorithm('dataset_11594_6')