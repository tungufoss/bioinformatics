def ProbabilityHiddenPathProblem(Sigma, States, Transition):
  p=0.5 #  an assumption that is modeled by setting Pr(pi0 -> pi1) = 1/2, where pi0 is a "silent" initial state that does not emit any symbols.
  for i in range(1,len(Sigma)):
    p*=Transition[(Sigma[i-1],Sigma[i])]  
  return p

def ProbabilityOutcomeGivenHPP(x, x_states, pi, pi_states, Emmision):
  p = 1  
  for i in range(len(x)):
    p*=Emmision[(pi[i],x[i])]
  
  return p
  
def ReadMatrix(inFile):
  colnames = inFile.readline().strip().split('\t')
  nCol = len(colnames)
  mat = {}
  for line in inFile:
    line = line.strip().split('\t')
    for i in range(nCol):
      mat[(line[0],colnames[i])]=float(line[i+1])
  return mat
  
def main_ProbabilityHiddenPathProblem(myfile):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  
  with open(inputFile) as inFile:        
    Sigma = inFile.readline().strip()
    inFile.readline() # --------
    States = inFile.readline().strip().split(' ')
    inFile.readline() # --------
    Transition = ReadTransition(inFile)
  
  prob = ProbabilityHiddenPathProblem(Sigma, States, Transition)
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
    Emmision = ReadMatrix(inFile)

  prob = ProbabilityOutcomeGivenHPP(x, x_states, pi, pi_states, Emmision)
  with open(outputFile, 'w') as outFile:    
    outFile.write(str(prob))
    print 'Probability that an HMM will emit a {} given its hidden path {}: {}'.format(x, pi, prob)

'''
main_ProbabilityHiddenPathProblem('sample_ProbabilityHiddenPathProblem')
main_ProbabilityHiddenPathProblem('dataset_11594_2')
'''
main_ProbabilityOutcomeGivenHPP('sample_ProbabilityOutcomeGivenHPP')
main_ProbabilityOutcomeGivenHPP('extra_ProbabilityOutcomeGivenHPP')
main_ProbabilityOutcomeGivenHPP('dataset_11594_4')