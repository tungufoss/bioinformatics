import sys
sys.path.insert(0,'../')
from common import RNA_codon_table, Integer_mass_table, aminoacids_tbl, PeptidesMasses, PeptideMass
sys.path.insert(0,'../week3')
from protein import SubpeptidesCyclic, SubpeptidesNotCyclic

def CyclopeptideScoring(peptide,spectrum,im_tbl,Cyclic):
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
    
  print '{} has score {} for spectrum'.format(peptide,score)
  return score
 
def main_CyclopeptideScoringProblem(myfile,Cyclic):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  with open(inputFile) as inFile:
    peptide = inFile.readline().strip()
    spectrum = [int(x) for x in inFile.readline().strip().split(' ')]
  
  im_tbl = Integer_mass_table()  
  score = CyclopeptideScoring(peptide, spectrum, im_tbl,Cyclic)  
  with open(outputFile,'w') as outFile:    
    outFile.write(str(score))
  
main_CyclopeptideScoringProblem('sample_cyclopeptidescoring',True) #11
main_CyclopeptideScoringProblem('cyclopeptide_scoring',True) #521
main_CyclopeptideScoringProblem('dataset_102_3',True)