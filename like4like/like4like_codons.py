from Bio.SeqUtils import seq1, seq3
import pandas as pd
from pathlib import Path
import pickle 
import bz2

# make something that wgets the updated databases from codonstatsdb.unr.edu/codonstatsdb_MONTHYEAR.tar.gz, but this is infrequent
# extract the .tar.gz file, delete the gene, ribosomal, plastid, and mitochrondrial stats, bottom decile, .tsvs
# delete .RData and .Rhistory 
# remove the 7 heading lines before tsv data ---- tail -n +9 file.tsv > file.tsv
# get taxID from directory numbering 

p = Path('codon_data')
tsv_paths = p.glob('**/*.tsv')
for tsv_path in tsv_paths:
    taxID = tsv_path.parts[1]
    df = pd.read_csv(tsv_path, sep='\t')
    print(taxID) 
    print(df['CODON'])
    print(df['Amino acid'])
    #AA_3let = df['Amino acid'] #Need to do this loop better
    #AA_1let = seq1(AA_3let)
    print(df['Fraction'])



#Make nested dictionary, the name of the dictionary is the taxID, then nest under AA (converted with seq1(), then by codon with Fraction

### Format like so to match python_codon_tables 
#taxID# = {
    #'*': {
        #'TAA': 0.47, 
        #'TAG': 0.23, 
        #'TGA': 0.3
    #}, 
    #'A': {
        #'GCA': 0.29, 
        #'GCC': 0.22, 
        #'GCG': 0.11, 
        #'GCT': 0.38}
    #}
#}


#pickle the python dictionaries for easy loading
#compress with bz2 for 1/4 size  
