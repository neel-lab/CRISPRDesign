import pandas 
import numpy as np
import azimuth
import azimuth.model_comparison
import sys

infile=sys.argv[1]
prefix=infile.split('.')[0]

#Import the data:
guideTable=pandas.read_csv(infile,sep='\t',header=None)

#Fix Reverse
table=str.maketrans({'A':'T','T':'A','G':'C','C':'G'})
guideTable[1]=[guideTable[1][x].translate(table)[::-1] if guideTable[6][x]=='R' else guideTable[1][x] for x in range(guideTable.shape[0])]

#Define Data Parameters:
seqs=guideTable[1].values
peptide_cut=guideTable[4].values
frac_peptide=guideTable[5].values

#Store the scores:
guideTable[7]=azimuth.model_comparison.predict(seqs,peptide_cut,frac_peptide)

#print(glycoCRISPR)
guideTable.to_csv(prefix+'_potency.txt',sep='\t',header=None,index=None)
