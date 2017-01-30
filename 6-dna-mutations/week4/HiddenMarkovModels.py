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
  
def ViterbiAlgorithm(states, sequence, transition, emission):  
  pi = []
  print sequence
  vcurr = {}  
  vprev = {}
  scurr = {}    
  svalues = {}  
  for i in range(len(sequence)+1):    
    for k in states:
      # Find the probabilility, if we are in state l, of choosing the nucleotide at position in the sequence              
      if i == 0 :        
        vcurr[k] = 1.0 * emission[(k,sequence[i])]        
      elif i < len(sequence) :      
        for l in states:
          svalues[l] = vprev[l] * transition[(l,k)] * emission[(k,sequence[i])]        
        scurr[k] = max(svalues, key=svalues.get)
        vcurr[k] = max(svalues.values())
      else : 
        for l in states:
          svalues[l] = vprev[l] * transition[(l,k)]
        scurr[k] = max(svalues, key=svalues.get)
        vcurr[k] = max(svalues.values())
        
    # Go through each of the rows of the matrix v (where each row represents a 
    # position in the DNA sequence), and find out which column has the maximum 
    # value for that row (where each column represents one state of the HMM):  
    if i>0 : 
      most_probable = max(vcurr, key=lambda k: vcurr[k])    
      pi.append(scurr[most_probable])      
    vprev=vcurr.copy()
  
  return pi
  
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
  
  pi = ViterbiAlgorithm(pi_states, x, Transition, Emission)
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
main_ViterbiAlgorithm('dataset_11594_6')