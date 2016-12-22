def aminoacids_tbl(im_tbl,unique):
  aminos = []
  weights = []
  tbl = {}
  with open('amino_acids.txt') as inFile:
    for line in inFile: 
      line = line.strip().split(' ')
      char = line[0]
      shrt = line[1]
      name = line[2]      
      mass = PeptideMass(char,im_tbl)
      if unique==False or mass not in weights : 
        aminos.append(char)
        weights.append(mass)
        tbl[char] = (shrt,name,mass)
  
  print '{} amino acids: {}'.format(len(aminos),''.join(aminos))
  return aminos,tbl
  
def RNA_codon_table():
  tbl = {}
  with open('RNA_codon_table_1.txt') as inFile:
    for line in inFile: 
      line = line.strip().split(' ')      
      rna = line[0]
      codon = '' if len(line) == 1 else line[1]
      tbl[rna]=codon
  return tbl

def Integer_mass_table():
  tbl = {}
  with open('integer_mass_table.txt') as inFile:
    for line in inFile: 
      line = line.strip().split(' ')
      integer = line[0]
      mass = int(line[1])
      tbl[integer]=mass
  return tbl
  
def TranslateRNAToCodon(txt,tbl):
  sol = ''
  i = 0  
  while i < len(txt):
    key=txt[i:(i+3)]
    if key in tbl : 
      sol += tbl[key]
    else : 
      sol += ' '
    i+=3
  return sol

def ReverseComplement(str):
  rev = ''
  switcher = {
      'G': 'C',
      'C': 'G',
      'T': 'A',
      'A': 'T',
    }
  for char in reversed(str):
    rev += switcher[char]
  return rev

def Transcribe(str):
  return str.replace("T", "U")
  
def PeptideEncoding(txt,codon,tbl):
  rev_txt = Transcribe(ReverseComplement(txt))
  trs_txt = Transcribe(txt)
  nCodon = 3*len(codon)  
  substrs = []
  n=len(txt)-nCodon
  
  for i in range(n):
    str = txt[i:(i+nCodon)]
    trs_str = trs_txt[i:(i+nCodon)]
    rev_str = rev_txt[(n-i):(n-i+nCodon)]
    rev_str = rev_txt[(n-i):(n-i+nCodon)]
        
    if TranslateRNAToCodon(trs_str,tbl) == codon : 
      substrs.append(str)
    elif TranslateRNAToCodon(rev_str,tbl) == codon: 
      substrs.append(str)
  
  return substrs
 
def SubpeptidesCyclic(peptide):
  subpeptides = ['',peptide]  
  N = len(peptide)
  peptide+=peptide
  for n in range(1,N):
    for i in range(N):      
      subpeptides.append(peptide[i:i+n])  
  return subpeptides

def SubpeptidesNotCyclic(peptide):
  subpeptides = ['',peptide]  
  N = len(peptide)  
  for n in range(1,N):
    for i in range(N-n+1):
      subpeptides.append(peptide[i:i+n])  
  return subpeptides
  
def PeptideMassAux(peptide,im_tbl):  
  mass = []
  for char in peptide:      
    mass.append(im_tbl[char])
  return mass

def PeptideMass(peptide,im_tbl):  
  mass = PeptideMassAux(peptide,im_tbl)  
  return sum(x for x in mass)
  
def Weigths(peptides,im_tbl):
  weights = {}
  for peptide in peptides:
    weights[peptide] = PeptideMass(peptide,im_tbl)
  return weights

def CyclopeptideSequencingExpandAndBound(peptides,aminos,cyclospec,im_tbl):
  expand={}
  for a in aminos: 
    for p in peptides:      
      pa = p+a      
      if ConsistentToSpectrum(pa,cyclospec[:],im_tbl,False):
        expand[pa]='-'.join([str(x) for x in PeptideMassAux(pa,im_tbl)])
  return expand
  
def ConsistentToSpectrum(peptide,spectrum,im_tbl,Cyclic):
  if Cyclic: 
    subpeptides = SubpeptidesCyclic(peptide)
  else :
    subpeptides = SubpeptidesNotCyclic(peptide)
    
  for sub in subpeptides:
    mass = PeptideMass(sub,im_tbl)
    if mass not in spectrum:      
      return False
    spectrum.remove(mass)  
  return True
  
def CyclopeptideSequencing(cyclospec):
  im_tbl = Integer_mass_table()
  L = len(cyclospec)
  print 'Spectrum:',cyclospec

  all_aminos, aminos_tbl = aminoacids_tbl(im_tbl,True)
  aminos=[]
  peptides={}
  for a in all_aminos:
    mass = PeptideMass(a,im_tbl)    
    if mass in cyclospec :       
      peptides[a] = mass
      aminos.append(a)
  
  while len(peptides)>0:    
    sol = peptides.copy()
    peptides = CyclopeptideSequencingExpandAndBound(peptides,aminos,cyclospec,im_tbl)        
    #print '#{}'.format(len(peptides)), ' '.join('{}:{}'.format(x,peptides[x]) for x in peptides)     
  
  return sol
      
def main_translateprotein(myfile):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  with open(inputFile) as inFile:
    txt = inFile.readline().strip()
  tbl = RNA_codon_table()
  sol = TranslateRNAToCodon(txt,tbl)
  with open(outputFile,'w') as outFile:
    outFile.write(sol)  

def main_peptideencoding(myfile):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  with open(inputFile) as inFile:
    txt = inFile.readline().strip()
    codon = inFile.readline().strip()
  tbl = RNA_codon_table()
  sol = PeptideEncoding(txt,codon,tbl)  
  with open(outputFile,'w') as outFile:
    outFile.write('\n'.join(sol))

def main_CountSubpeptides(myfile,IsCyclic):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  with open(inputFile) as inFile:
    n = int(inFile.readline().strip())  
  if IsCyclic: 
    N=n*(n-1)
  else : 
    N=n*(n+1)/2 + 1
  print '{} length peptide has {} {} combinations'.format(n,N,'cyclic' if IsCyclic  else 'linear')
  with open(outputFile,'w') as outFile:
    outFile.write(str(N))

    
def ComputerNumberCompomers(M):
  import numpy as np
  import math  
  
  im_tbl = Integer_mass_table()    
  PEP_NAMES, aminos_tbl = aminoacids_tbl(im_tbl,False)
  
  PEP_MASSES = []  
  for a in PEP_NAMES :    
    PEP_MASSES.append(PeptideMass(a,im_tbl))
  LEN_PEP_MASSES = len(PEP_MASSES)
  NUM_COMB       = 2**LEN_PEP_MASSES-1
  
  # numpy array containing all possible coeficients  
  C = np.array([[int(x) for x in np.binary_repr(K, width=LEN_PEP_MASSES)] for K in xrange(NUM_COMB)])
  # each element is an array of coefficients representing a number between 0 and NUM_COMB in binary form
  print "type(C)      = ",type(C)
  print "type(C[0])   = ",type(C[0])
  print "C.shape      = ",str(C.shape)
  print "C[0].shape   = ",str(C[0].shape)
  print "C[0]         = ",C[0]
  print "C[15]        = ",C[15]
  print "C[255]       = ",C[255]
  
  # Calculate sum of all combinations
  PROD = C.dot(PEP_MASSES)

  # find the ones that match M
  valid_combinations = [(i,x) for i,x in enumerate(PROD) if x == M]
  print 'Found {} possibilities with total mass = {}'.format(len(valid_combinations), M)
  #print valid_combinations
  total=0  
  for comb_index, comb_mass in valid_combinations:
    # work back the combinations in string format
    comb_str = [PEP_NAMES[i] for i,x in enumerate(C[comb_index]) if x==1]    
    total+= math.factorial(len(comb_str))
    #print '{} --> {}'.format(comb_index, ''.join(comb_str))
  
  print 'Total permutations: {}'.format(total)
  return total
    
def main_CountSubpeptidesWithMass(myfile):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  with open(inputFile) as inFile:
    mass = int(inFile.readline().strip())
  
  total = ComputerNumberCompomers(mass)
  with open(outputFile,'w') as outFile:        
    outFile.write(str(total))    
      
def main_spectrum(myfile,Cyclic):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  with open(inputFile) as inFile:
    peptide = inFile.readline().strip()  
  
  if Cyclic: 
    subpeptides = SubpeptidesCyclic(peptide)
  else :
    subpeptides = SubpeptidesNotCyclic(peptide)

  im_tbl = Integer_mass_table()
  masses = Weigths(subpeptides,im_tbl)    
  with open(outputFile,'w') as outFile:    
    txt = ' '.join([str(m) for m in sorted([masses[x] for x in masses])])    
    outFile.write(txt)
  
def main_CyclopeptideSequencing(myfile):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  with open(inputFile) as inFile:
    spectrum = [int(x) for x in inFile.readline().strip().split(' ')]

  sol = CyclopeptideSequencing(spectrum)
  
  with open(outputFile,'w') as outFile:    
    txt=' '.join(sorted('{}'.format(sol[x]) for x in sol))
    print '#{} solutions: {}\n'.format(len(sol),txt)
    outFile.write(txt)

def CountStringsTranscribeToTxt(txt,rna_tbl):
  from numpy import prod
  cnt=[]
  for t in txt:
    cnt.append(sum([rna_tbl[x] == t for x in rna_tbl]))
  total=prod(cnt)
  print '{}: #{}={}'.format(txt,cnt,total)
  return total
    
def main_test():  
  rna_tbl = RNA_codon_table()
  im_tbl = Integer_mass_table()
  all_aminos, aminos_tbl = aminoacids_tbl(im_tbl,False)
  
  print '\n#1 True or False: Tyrocidine B1 is synthesized by NRP synthetase.'
  print True
  
  print '\n#2 Which of the following RNA strings could translate into the amino acid string PRTEIN?'
  for txt in ['CCGAGGACCGAAAUCAAC','CCCCGUACGGAGAUGAAA','CCCAGGACUGAGAUCAAU','CCCAGUACCGAAAUUAAC'] :
    translated=TranslateRNAToCodon(txt,rna_tbl)
    print '{}->{} {}'.format(txt,translated,translated=='PRTEIN')

  print '\n#3 How many DNA strings transcribe and translate into the amino acid string MASS?'  
  txt='MASS'
  total = CountStringsTranscribeToTxt(txt,rna_tbl)
    
  print '\n#4 What is the integer mass of glutamine?'
  # Glutamine (abbreviated as Gln or Q; encoded by the codons CAA and CAG)  
  char = [a for a in all_aminos if aminos_tbl[a][0] == 'Gln'][0]
  print '{}: {} mass'.format(char,PeptideMass(char,im_tbl))
  
  print '\n#5 Which of the following cyclic peptides could have generated the theoretical spectrum'
  theoretical_spectrum=[0,71,101,113,131,184,202,214,232,285,303,315,345,416]
  for peptide in ['MTAI','TMLA','TLAM','TMIA','TAIM','MAIT']:    
    if ConsistentToSpectrum(peptide,theoretical_spectrum[:],im_tbl,True): 
      print '{} matches spectrum'.format(peptide)
    else :
      print '{} does not match spectrum'.format(peptide)
 
  print '\n#6 Which of the following linear peptides is consistent with Spectrum?'
  theoretical_spectrum=[0,71,99,101,103,128,129,199,200,204,227,230,231,298,303,328,330,332,333]
  for peptide in ['QCV','TCQ','AQV','TCE','VAQ','CTV']:    
    if ConsistentToSpectrum(peptide,theoretical_spectrum[:],im_tbl,False): 
      print '{} matches spectrum'.format(peptide)
    else :
      print '{} does not match spectrum'.format(peptide)    
  

def main_cnt():
  rna_tbl = RNA_codon_table()
  im_tbl = Integer_mass_table()
  all_aminos, aminos_tbl = aminoacids_tbl(im_tbl,False)
  amino=''
  for shrt in ['Val','Lys','Leu','Phe','Pro','Trp','Phe','Asn','Gln','Tyr'] :
    amino+=[a for a in all_aminos if aminos_tbl[a][0] == shrt][0]  
  CountStringsTranscribeToTxt(amino,rna_tbl)  
  
  #Bacillus_brevis
  

'''
main_translateprotein('sample_translateprotein')
main_translateprotein('dataset_96_4')
main_peptideencoding('sample_peptideencoding')
main_peptideencoding('dataset_96_7')
main_CountSubpeptides('sample_cntsubpeptides_cyclic',1)
main_CountSubpeptides('dataset_98_3',1)
main_spectrum('sample_spectrum',True)
main_spectrum('dataset_98_4',True)
main_CountSubpeptidesWithMass('sample_mass') # not working properly 
main_CountSubpeptides('sample_cntsubpeptides_path',0)
main_CountSubpeptides('dataset_100_3',0)
main_CyclopeptideSequencing('sample_cyclopeptideseq')
main_CyclopeptideSequencing('sample_cyclopeptideseq2')
main_CyclopeptideSequencing('leaderboard_spectrum')
main_CyclopeptideSequencing('dataset_100_6')
main_spectrum('sample_linearspectrum',False)
main_spectrum('dataset_4912_2',False)
'''
#main_test()
main_cnt()
