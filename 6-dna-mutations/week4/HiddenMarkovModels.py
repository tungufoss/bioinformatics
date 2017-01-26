def ProbabilityHiddenPathProblem(Sigma, States, Transition):
  p=0.5 #  an assumption that is modeled by setting Pr(pi0 -> pi1) = 1/2, where pi0 is a "silent" initial state that does not emit any symbols.
  for i in range(1,len(Sigma)):
    p*=Transition[(Sigma[i-1],Sigma[i])]  
  return p

def main_ProbabilityHiddenPathProblem(myfile):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  
  with open(inputFile) as inFile:        
    Sigma = inFile.readline().strip()
    inFile.readline() # --------
    States = inFile.readline().strip().split(' ')
    inFile.readline() # --------
    colnames = inFile.readline().strip().split('\t')
    nCol = len(colnames)
    Transition = {}
    for line in inFile:
      line = line.strip().split('\t')
      for i in range(nCol):
        Transition[(line[0],colnames[i])]=float(line[i+1])
  
  prob = ProbabilityHiddenPathProblem(Sigma, States, Transition)
  with open(outputFile, 'w') as outFile:    
    outFile.write(str(prob))
    print 'Probability of a Hidden Path Problem for {}: {}'.format(Sigma, prob)

main_ProbabilityHiddenPathProblem('sample_ProbabilityHiddenPathProblem')
main_ProbabilityHiddenPathProblem('dataset_11594_2')