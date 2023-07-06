import Bio
from Bio.PDB import PDBList

PDB_IDs = ['1NKD','4DO2','3K79','1RPO','7KAE','1RPO','2IJH','1GTO','2IJJ','2IJK','2GHY','2IJI','1F4M','1YO7','3F0O','5DSF','5U7A','5C0T','3FN8','1Q5V','1NEK','2WS3','2WDR','2WP9','1KFY','3CIR','1L0V','2B76']


for PDB_ID in PDB_IDs:
    PDBList().retrieve_pdb_file(PDB_ID, file_format='pdb', pdir='.', overwrite=True)
