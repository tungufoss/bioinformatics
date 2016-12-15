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
        
    if TranslateRNAToCodon(trs_str,tbl) == codon : 
      substrs.append(str)
    elif TranslateRNAToCodon(rev_str,tbl) == codon: 
      substrs.append(str)
  
  return substrs
 
def Subpeptides(peptide):
  subpeptides = ['',peptide]  
  N = len(peptide)
  peptide+=peptide
  for n in range(1,N):
    for i in range(N):      
      subpeptides.append(peptide[i:i+n])  
  return subpeptides

def PeptideMass(peptide,im_tbl):  
  mass = 0
  for char in peptide:      
    mass += im_tbl[char]
  return mass

def Weigths(peptides,im_tbl):
  weights = {}
  for peptide in peptides:
    weights[peptide] = PeptideMass(peptide,im_tbl)
  return weights


def CyclopeptideSequencing(cyclospec):
  import itertools
  from math import sqrt
  parent_mass = cyclospec[-1]
  sol=[]  
  # Let n be the length of a given peptide, and L be the length of its cyclospectrum.  Then L = n(n-1) + 2.
  # Using the quadratic formula to to solve for n:  n = (sqrt(4L-7) + 1)/2
  n = int((sqrt(4*len(cyclospec)-7)+1)/2)  
  print n, len(cyclospec)
  for comb in itertools.combinations(cyclospec[1:-1], n): 
    total_mass = 0
    for mass in comb :
      total_mass += mass      
    if total_mass == parent_mass :  
      perm = itertools.permutations(comb, n)
      for peptide in perm:          
        sol.append('-'.join([str(x) for x in peptide]))
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

def main_CountSubpeptidesWithMass(myfile):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  with open(inputFile) as inFile:
    mass = int(inFile.readline().strip())  
  im_tbl = Integer_mass_table()
  
  aminoacids = 'GASPVTCILNDKQEMHFRYW'  #this are all the aminoacids  
  print mass, aminoacids
  
  maxLength = 100
  import itertools
  total = 0
  pepList = []
  for i in range(maxLength+1):
    for p in itertools.combinations_with_replacement(aminoacids, i): 
    #order matters for the total number of peptides but not for calculating the total mass
      amino = ''.join(p)      
      if PeptideMass(amino, im_tbl) == mass:
        pepList.append(amino)
    print 'Iter: {} |{}|'.format(i,len(pepList))

  newpepList = []
  for i in pepList:
    for p in itertools.permutations(i, r = len(i)): 
    #I use permutations here to get the total number because order matters
      if p not in newpepList:
        newpepList.append(p)
        total +=1

  print (total)
      
def main_spectrum(myfile):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  with open(inputFile) as inFile:
    peptide = inFile.readline().strip()  
  subpeptides = Subpeptides(peptide)
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
    txt= '\n'.join(sol)
    print txt
    outFile.write(txt)
    
'''
main_translateprotein('sample_translateprotein')
main_translateprotein('dataset_96_4')
main_peptideencoding('sample_peptideencoding')
main_peptideencoding('dataset_96_7')
main_CountSubpeptides('sample_cntsubpeptides_cyclic',1)
main_CountSubpeptides('dataset_98_3',1)
main_spectrum('sample_spectrum')
main_spectrum('dataset_98_4')
main_CountSubpeptidesWithMass('sample_mass') # not working properly 
main_CountSubpeptides('sample_cntsubpeptides_path',0)
main_CountSubpeptides('dataset_100_3',0)
'''

main_CyclopeptideSequencing('sample_cyclopeptideseq')
main_CyclopeptideSequencing('dataset_100_6')
