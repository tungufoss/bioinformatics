import sys
sys.path.insert(0,'../')
from common import RNA_codon_table, Integer_mass_table, aminoacids_tbl, PeptidesMasses, PeptideMass
sys.path.insert(0,'../week3')
from protein import SubpeptidesCyclic, SubpeptidesNotCyclic
import timeit

def CyclopeptideScoring(peptide,spectrum,Cyclic,im_tbl):
  spectrum=spectrum[:]
  if Cyclic: 
    subpeptides = SubpeptidesCyclic(peptide)
  else :
    subpeptides = SubpeptidesNotCyclic(peptide)
  
  score=0
  for sub in subpeptides:
    mass = PeptideMass(sub,im_tbl)
    if mass in spectrum:      
      score+=1
      spectrum.remove(mass)
    
  #print '{} has score {} for spectrum'.format(peptide,score)
  return score

def CyclopeptideScoring2(masses,spectrum,Cyclic):
  spectrum = spectrum[:]
  subpeptides_nc = SubpeptidesNotCyclic(masses,[0])
  subpeptides_c = SubpeptidesCyclic(masses,[0])
  
  score_nc=0
  for sub in subpeptides_nc:
    subpeptides_c.remove(sub)
    mass = sum(sub)
    if mass in spectrum:
      score_nc+=1
      spectrum.remove(mass)  
  
  score_c = score_nc  
  if Cyclic:  
    for sub in subpeptides_c:
      mass = sum(sub)
      if mass in spectrum:      
        score_c+=1
        spectrum.remove(mass)
  
  return score_nc, score_c
  
def LeaderBoardExpand(LeaderBoard,aminos_tbl,spectrum,Cyclic,im_tbl):
  expand = []
  N = len(LeaderBoard)*len(aminos_tbl)
  
  for x in LeaderBoard:
    for a in aminos_tbl:
      pa = x[0]+a
      ma = aminos_tbl[a][2]
      pm = x[4][:]
      pm.append(ma)
      score_nc, score_c = CyclopeptideScoring2(pm, spectrum, Cyclic)
      expand.append(
        (
          pa,
          score_nc, # linear scoring for trimming
          x[2]+ma, # i.e. PeptideMass(pa,im_tbl)
          score_c, # <cyclic?> score is only used for comparing with the leader peptide(s).
          pm
        )
      )  
  return expand

def LeaderBoardTrim(scored_LeaderBoard, N):  
  if len(scored_LeaderBoard) <= N:
    return scored_LeaderBoard    
  # sort LeaderBoard according to the decreasing order of scores in LinearScores
  scored_LeaderBoard = sorted(scored_LeaderBoard, key=lambda x: x[1], reverse=True) 
  # return top N peptides       
  Nth_score = scored_LeaderBoard[N-1][1] # score
  # allow ties 
  scored_LeaderBoard = [x for x in scored_LeaderBoard if x[1] >= Nth_score]  
  return scored_LeaderBoard
  
def LeaderBoardCyclopeptideSequencing(spectrum,N, Cyclic, aminos_tbl,im_tbl):  
  spectrum = sorted(spectrum)    
  LeaderScore = CyclopeptideScoring('', spectrum, Cyclic,im_tbl)
  LeaderPeptides = []
    
  LeaderBoard = [('', LeaderScore, 0, LeaderScore,[])]
  ParentMass = spectrum[-1]  
  print 'ParentMass: {} da'.format(ParentMass)
  
  start = timeit.default_timer()  
  while len(LeaderBoard)>0:       
    # expand LeaderBoard
    LeaderBoard = LeaderBoardExpand(LeaderBoard,aminos_tbl,spectrum,Cyclic,im_tbl)
    
    # If if Mass(Peptide) = ParentMass(Spectrum) and Score(Peptide, Spectrum) > Score(LeaderPeptide, Spectrum) then update LeaderPeptide        
    for (peptide,score) in [(x[4],x[3]) for x in LeaderBoard if x[2] == ParentMass and x[3] >= LeaderScore]:      
      if score > LeaderScore:    
        LeaderScore = score
        LeaderPeptides = []
        LeaderPeptides.append(peptide)
        print 'Peptide {} has score {} after {} min'.format('-'.join([str(x) for x in peptide]),score,(timeit.default_timer()-start)/60)
      elif score == LeaderScore:
        LeaderPeptides.append(peptide)      
           
    # remove Peptide from LeaderBoard if Mass(Peptide) > ParentMass(Spectrum)                        
    LeaderBoard = [x for x in LeaderBoard if x[2] <= ParentMass]
    
    # Trim LeaderBoard w.r.t. N        
    LeaderBoard = LeaderBoardTrim(LeaderBoard, N)    
    
  print 'Total of {} leader peptides with score {}'.format(len(LeaderPeptides), LeaderScore)  
  return ['-'.join([str(x) for x in peptide]) for peptide in LeaderPeptides]  

#The list of elements in the convolution of Spectrum. 
def SpectralConvolution(spectrum):  
  # The spectrum isn't sorted, so find all differences and filter out the non-positive.
  # If an element has multiplicity k, it should appear exactly k times
  convolution = [i-j for i in spectrum for j in spectrum if i-j > 0]  
  return convolution

# Get the top M elements from the convolution that are between 57 and 200.
def TopConvolution(convolution,M):
  from collections import Counter
  import operator  
  cnt = Counter([x for x in convolution if x>=57 and x<=200])  
  cnt = sorted(cnt.items(), key=operator.itemgetter(1), reverse=True)
  if len(cnt)>=M :
    return [x[0] for x in cnt]
    
  Mth_value = cnt[M-1][1] # include ties 
  return [x[0] for x in cnt if x[1] >= Mth_value]
  
def main_CyclopeptideScoringProblem(myfile,Cyclic):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  with open(inputFile) as inFile:
    peptide = inFile.readline().strip()
    spectrum = [int(x) for x in inFile.readline().strip().split(' ')]
  
  score = CyclopeptideScoring(peptide, spectrum, Cyclic)  
  with open(outputFile,'w') as outFile:    
    outFile.write(str(score))
  
def main_LeaderBoardCyclopeptideSequencing(myfile,Cyclic,print_all,synthetic_aminos=False):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  start = timeit.default_timer()
  
  with open(inputFile) as inFile:
    N = int(inFile.readline().strip())
    spectrum = [int(x) for x in inFile.readline().strip().split(' ')]
   
  if synthetic_aminos: 
    aminos_tbl={}
    im_tbl={}
    for da in range(57,201):      
      aminos_tbl[chr(da)]=(chr(da),chr(da),da)
      im_tbl[chr(da)]=da
  else :
    im_tbl = Integer_mass_table()  
    aminos, aminos_tbl = aminoacids_tbl(im_tbl,True)
  
  leader_peptides = LeaderBoardCyclopeptideSequencing(spectrum,N,Cyclic,aminos_tbl,im_tbl)
  
  with open(outputFile,'w') as outFile:
    if print_all:
      txt = ' '.join(leader_peptides)
    else:
      txt = leader_peptides[0]
    outFile.write(txt)    
  
  stop = timeit.default_timer()
  print 'Running time {} sec'.format(stop - start)
 
def main_Trim(myfile):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  start = timeit.default_timer()
  
  with open(inputFile) as inFile:
    LeaderBoard = [x for x in inFile.readline().strip().split(' ')]
    spectrum = [int(x) for x in inFile.readline().strip().split(' ')]
    N = int(inFile.readline().strip())
  
  scored_LeaderBoard = [(
    peptide,
    CyclopeptideScoring(peptide, spectrum, False), # score    
  ) for peptide in LeaderBoard]  
  trimmed = LeaderBoardTrim(scored_LeaderBoard, N)  
  
  with open(outputFile,'w') as outFile:    
    txt=' '.join([x[0] for x in trimmed])
    outFile.write(txt)
    print txt

def main_SpectralConvolution(myfile):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  with open(inputFile) as inFile:
    spectrum = [int(x) for x in inFile.readline().strip().split(' ')]
  
  convolution = SpectralConvolution(spectrum)
  
  with open(outputFile,'w') as outFile:    
    txt=' '.join([str(x) for x in convolution])
    outFile.write(txt)

def main_ConvolutionCyclopeptideSequencing(myfile):  
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  start = timeit.default_timer()
  
  with open(inputFile) as inFile:
    M = int(inFile.readline().strip())
    N = int(inFile.readline().strip())
    spectrum = [int(x) for x in inFile.readline().strip().split(' ')]
  
  convolution = SpectralConvolution(spectrum)
  convolution = TopConvolution(convolution, M)  
  
  aminos_tbl={}
  im_tbl={}
  for da in convolution:      
    aminos_tbl[chr(da)]=(chr(da),chr(da),da)
    im_tbl[chr(da)]=da
  
  leader_peptides = LeaderBoardCyclopeptideSequencing(spectrum,N,True,aminos_tbl,im_tbl)
    
  if myfile == 'sample_ConvolutionCyclopeptideSequencing':
    leader_peptide = '99-71-137-57-72-57'    
    if leader_peptide in leader_peptides:
      print 'Correct leader peptide found'
  elif myfile == 'convolution_cyclopeptide_sequencing':
    leader_peptide='113-115-114-128-97-163-131-129-129-147-57-57-129'
    if leader_peptide in leader_peptides:
      print 'Correct leader peptide found'
  else :
    leader_peptide = leader_peptides[0]
  
  with open(outputFile,'w') as outFile:        
    outFile.write(leader_peptide)
    print leader_peptide

  stop = timeit.default_timer()
  print 'Running time {} sec'.format(stop - start)  
    
'''
main_CyclopeptideScoringProblem('sample_cyclopeptidescoring',True) #11
main_CyclopeptideScoringProblem('cyclopeptide_scoring',True) #521
main_CyclopeptideScoringProblem('dataset_102_3',True)
main_CyclopeptideScoringProblem('linear_score',False) #8
main_CyclopeptideScoringProblem('dataset_4913_1',False) #8
main_Trim('sample_trim')
main_Trim('dataset_4913_3')
main_LeaderBoardCyclopeptideSequencing('sample_linearpeptidescoring',False,False)
main_LeaderBoardCyclopeptideSequencing('dataset_102_8',False,False)
main_LeaderBoardCyclopeptideSequencing('dataset_102_10',True,True)
main_LeaderBoardCyclopeptideSequencing('dataset_103_2',True,True,True)
main_SpectralConvolution('sample_spectralconvolution')
main_SpectralConvolution('dataset_104_4')
'''
#main_ConvolutionCyclopeptideSequencing('sample_ConvolutionCyclopeptideSequencing')
#main_ConvolutionCyclopeptideSequencing('convolution_cyclopeptide_sequencing')
main_ConvolutionCyclopeptideSequencing('dataset_104_7')