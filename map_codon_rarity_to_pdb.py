import argparse
import Bio 
from Bio import SeqIO 
from Bio.PDB import PDBList, PDBParser, PDBIO
from Bio.SeqUtils import seq1, seq3 
from pathlib import Path
import pickle
import numpy
import math

parser = argparse.ArgumentParser(prog='map  rare codons', description='testing for mapping rare codon use in Ebola and HepB viral n    ucleic acid sequences to AF2 .pdb structures')
#parser.add_argument('--debug',  type=bool, help='To print values when something goes wrong') #TO DO, wrap print statements in debug 
parser.add_argument('-nuc', '--nucseq',  type=str, help='A manually entered nucleic acid sequence to split into codons', required=True) 
parser.add_argument('-pdb', '--PDBID', type=str, help='.pdb file', required=True)
parser.add_argument('-tax', '--taxID',  type=str, help='TaxID to load the relevant codon rarity table', required=True)
parser.add_argument('-ct', '--codontables',  type=str, help=f'Pickle file (usually .pkl) of the codon tables, organised by taxID', default='codon_tables.pkl', required=True)  # Might bundle with it later 
parser.add_argument('-o', '--outfile', type=str, help='name out annotate .pdb output file', required=True)
args = parser.parse_args()

def load_codon_table(codon_tables, taxID):  #I can probably have the codon table inbuilt
    with open(codon_tables, 'rb') as pickle_file:
        codon_tables = pickle.load(pickle_file)
    return codon_tables[taxID]

def pdbfile2struct(pdbfile): # TO DO: consider file -> pdbfile variable rename. Can Bio.PDB help here with reading (PDBIO.select pulling out 'REMARK', etc)
    struct_name = str(pdbfile).replace('.pdb', '')
    pdb_struct = PDBParser().get_structure(struct_name, pdbfile)
    return pdb_struct
    
def nucseq2codons(nucseq):
    codons = [nucseq[i:i+3] for i in range(0, len(nucseq), 3)]
    return codons

def replace_b_factor(pdb_struct, codons, codon_table): 
    a_idx = 0  # no neat way to link idx to enumerate
    res_idx = 0
    for model in pdb_struct: #doesn't hurt, but check later if there are pdbs with multiple models
        for chain in model:
            for residue in chain:
                res_idx += 1 
                res_3let = residue.get_resname()
                res_1let = seq1(res_3let)
                print(f'Residue index: {res_idx}, AA is: {res_1let}, Codon is: {codons[res_idx-1]}')
                for atom in residue:
                    codon_rarity_val = codon_table[res_1let][codons[res_idx-1]] 
                    atom.set_bfactor(codon_rarity_val) #this is used in AlphaFold to colour residue position certainty
   
    
pdbfile = Path(args.PDBID)
codon_table = load_codon_table(args.codontables, args.taxID)
print(f'Codon fraction table for {args.taxID} is {codon_table}')
pdb_struct = pdbfile2struct(pdbfile)
if args.nucseq is not None: #Need to decide which overrides
    seq_record = SeqIO.read(args.nucseq, 'fasta') 
    nucseq = str(seq_record.seq).upper()
    print(f'The nucleic acid sequence is: {nucseq}')
codons = nucseq2codons(nucseq)
print(codons)
replace_b_factor(pdb_struct, codons, codon_table) 
io=PDBIO()
io.set_structure(pdb_struct) 
io.save(args.outfile)
