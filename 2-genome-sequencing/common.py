def RNA_codon_table():
  tbl = {}
  with open('../RNA_codon_table_1.txt') as inFile:
    for line in inFile: 
      line = line.strip().split(' ')      
      rna = line[0]
      codon = '' if len(line) == 1 else line[1]
      tbl[rna]=codon
  return tbl

def Integer_mass_table():
  tbl = {}
  with open('../integer_mass_table.txt') as inFile:
    for line in inFile: 
      line = line.strip().split(' ')
      integer = line[0]
      mass = int(line[1])
      tbl[integer]=mass
  return tbl

def aminoacids_tbl(im_tbl,unique):
  aminos = []
  weights = []
  tbl = {}
  with open('../amino_acids.txt') as inFile:
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

def PeptideMassAux(peptide,im_tbl):  
  mass = []
  for char in peptide:      
    mass.append(im_tbl[char])
  return mass

def PeptideMass(peptide,im_tbl):  
  mass = PeptideMassAux(peptide,im_tbl)  
  return sum(x for x in mass)