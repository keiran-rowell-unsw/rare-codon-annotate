#This script will generate the codon tables for a local extracted copy of codonstatsdb.unr.edu and put them into a pickle file (codon_tables.pkl) 
from Bio.SeqUtils import seq1, seq3
import copy
import pandas as pd
from pathlib import Path
import pickle 


#def compareM_to_condon_table(compareM_out.csv): # TO DO: if someone has a new species I can provide a processing pathway from the raw frequency output from compareM.

# Steph has given me the data for 'Nelly' the archeon here: https://unsw-my.sharepoint.com/personal/z5062530_ad_unsw_edu_au/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Fz5062530%5Fad%5Funsw%5Fedu%5Fau%2FDocuments%2FResearch%2FCodon%5FUsage%5F4%5FKeiran&g1


# TO DO: 

# May add more codon table templates later: codon code can vary by organism: https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG1, 
# Codon table date: 7 Jan 2019 
# If I was to do that, I could just use 'differences from standard code' and modify a deepcopy of the standard inverse codon table
standard_inverse_codon_table = {
'*': {'TAA': 0.0, 'TAG': 0.0, 'TGA': 0.0}, 
'A': {'GCA': 0.0, 'GCC': 0.0, 'GCG': 0.0, 'GCT': 0.0}, 
'C': {'TGC': 0.0, 'TGT': 0.0}, 
'D': {'GAC': 0.0, 'GAT': 0.0}, 
'E': {'GAA': 0.0, 'GAG': 0.0}, 
'F': {'TTC': 0.0, 'TTT': 0.0}, 
'G': {'GGA': 0.0, 'GGC': 0.0, 'GGG': 0.0, 'GGT': 0.0}, 
'H': {'CAC': 0.0, 'CAT': 0.0}, 
'I': {'ATA': 0.0, 'ATC': 0.0, 'ATT': 0.0},
'K': {'AAA': 0.0, 'AAG': 0.0},
'L': {'CTA': 0.0, 'CTC': 0.0, 'CTG': 0.0, 'CTT': 0.0, 'TTA': 0.0, 'TTG': 0.0},
'M': {'ATG': 0.0},
'N': {'AAC': 0.0, 'AAT': 0.0},
'P': {'CCA': 0.0, 'CCC': 0.0, 'CCG': 0.0, 'CCT': 0.0},
'Q': {'CAA': 0.0, 'CAG': 0.0}, 
'R': {'AGA': 0.0, 'AGG': 0.0, 'CGA': 0.0, 'CGC': 0.0, 'CGG': 0.0, 'CGT': 0.0},
'S': {'AGC': 0.0, 'AGT': 0.0, 'TCA': 0.0, 'TCC': 0.0, 'TCG': 0.0, 'TCT': 0.0}, 
'T': {'ACA': 0.0, 'ACC': 0.0, 'ACG': 0.0, 'ACT': 0.0}, 
'V': {'GTA': 0.0, 'GTC': 0.0, 'GTG': 0.0, 'GTT': 0.0}, 
'W': {'TGG': 0.0}, 
'Y': {'TAC': 0.0, 'TAT': 0.0}
}

taxID_codon_tables = {}

p = Path('codon_data')
tsv_paths = p.glob('**/*.tsv')

for tsv_path in tsv_paths:
    taxID = tsv_path.parts[1]
    print(taxID)
    #print(f'Creating codon table for: {taxID}') # TO DO: consider adding progress bar: https://github.com/rsalmei/alive-progress
    #could add logic here to check which codon table applies based upon taxID
    taxID_codon_table = copy.deepcopy(standard_inverse_codon_table)  # deepcopy for new dict per taxID, not a reference that overwrites 
    df = pd.read_csv(tsv_path, sep='\t')
    codon_frac_zip = zip(df['CODON'], df['Amino acid'], df['Fraction']) 

    for codon, AA_3let, fraction in codon_frac_zip:
        AA_1let = seq1(AA_3let) #This does handle 'Stop'
        if AA_1let == 'X':
            AA_1let = '*' #Used as the termination charactered in python_codon_tables
        frac = round(fraction, 2)  #May have rounding errors so sum != 1. So far, off by 0.01--0.02 from python_codon_tables
        taxID_codon_table[AA_1let][codon] = frac
        taxID_codon_tables[taxID] = taxID_codon_table

print("Here's the table for nelly:")
print(taxID_codon_tables['nelly'])

pickle.dump(taxID_codon_tables, open('codon_tables.pkl', 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
