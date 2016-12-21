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
      if ConsistentToSpectrum(pa,cyclospec[:],im_tbl):        
        expand[pa]='-'.join([str(x) for x in PeptideMassAux(pa,im_tbl)])
  return expand
  
def ConsistentToSpectrum(peptide,spectrum,im_tbl):  
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

  all_aminos='GASPVTCINDKEMHFRYW' #this are all the aminoacids
  aminos=[]
  peptides={}
  for a in all_aminos:
    mass = PeptideMass(a,im_tbl)    
    if mass in cyclospec :       
      peptides[a] = mass
      aminos.append(a)
  aminos=all_aminos  
  
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

def main_CountSubpeptidesWithMass(myfile):
  inputFile = myfile + '.txt'
  outputFile = myfile + '.out'
  with open(inputFile) as inFile:
    mass = int(inFile.readline().strip())  
  im_tbl = Integer_mass_table()
  
  aminoacids = 'GASPVTCILNDKQEMHFRYW'  #this are all the aminoacids  
  print mass, aminoacids, len(aminoacids)
  
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
'''

main_spectrum('sample_linearspectrum',False)
main_spectrum('dataset_4912_2',False)
