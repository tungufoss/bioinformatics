def RNA_codon_table():
  tbl = {}
  with open('RNA_codon_table_1.txt') as inFile:
    for line in inFile: 
      line = line.strip().split(' ')      
      rna = line[0]
      codon = '' if len(line) == 1 else line[1]
      tbl[rna]=codon
  return tbl

def TranslateRNAToCodon(txt,tbl):
  sol = ''
  i = 0
  n = len(txt)-3
  while i < n:
    key=txt[i:(i+3)]
    sol += tbl[key]
    i+=3
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
  
main_translateprotein('sample_translateprotein')
main_translateprotein('dataset_96_4')
