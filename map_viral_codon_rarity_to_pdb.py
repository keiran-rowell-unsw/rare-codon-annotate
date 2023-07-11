import argparse
import Bio 
from Bio import SeqIO 
from Bio.PDB import PDBList, PDBParser
from Bio.SeqUtils import seq1, seq3 
from pathlib import Path
import pickle
import numpy
import math

parser = argparse.ArgumentParser(prog='map viral rare codons', description='testing for mapping rare codon use in Ebola and HepB viral n    ucleic acid sequences to AF2 .pdb structures')
parser.add_argument('-nuc', '--nucseq',  type=str, help='A manually entered nucleic acid sequence to split into codons ') 
parser.add_argument('-pdb', '--PDBID', type=str, help='RSCB Protein Data Bank ID to fetch .pdb file')
parser.add_argument('-tax', '--taxID',  type=str, help='TaxID to load the relevant codon rarity table', required=True)
parser.add_argument('-ct', '--codontables',  type=str, help=f'Pickle file (usually .pkl) of the codon tables, organised by taxID', default='codon_tables.pkl')
parser.add_argument('-o', '--outfile', type=str, help='name out annotate .pdb output file')
args = parser.parse_args()

def load_codon_table(codon_tables, taxID):  #I can probably have the codon table inbuilt
    with open(codon_tables, 'rb') as pickle_file:
        codon_tables = pickle.load(pickle_file)
    return codon_tables[taxID]

def pdbfile2struct(pdbfile): # TO DO: consider file -> pdbfile variable rename. Can Bio.PDB help here with reading (PDBIO.select pulling out 'REMARK', etc)
    parser = PDBParser()
    struct_name = str(pdbfile).replace('.pdb', '')
    pdb_struct = parser.get_structure(struct_name, pdbfile)
    return pdb_struct
    
def nucseq2codons(nucseq):
    codons = [nucseq[i:i+3] for i in range(0, len(nucseq), 3)]
    return codons

def replace_b_factor(pdb_struct, codons, codon_table): 
    new_pdb_contents = [] 
    a_idx = 0  # no neat way to link idx to enumerate
    res_idx = 0
    for model in pdb_struct: #doesn't hurt, but check later if there are pdbs with multiple models
        for chain in model:
            for residue in chain:
                res_idx += 1 
                res_3let = residue.get_resname()
                res_1let = seq1(res_3let)
                #print(res_1let)
                for atom in residue:
                    a_idx +=1
                    codon_rarity_val = codon_table[res_1let][codons[res_idx-1]] 
                    # Set up variables for printing the ATOM line
                    a_name = atom.get_name()
                    chain_id = chain.get_id().replace('<Chain id=','').replace('>','')
                    a_coords = numpy.array2string(atom.get_coord()).replace('[','').replace(']','')
                    a_occ = atom.get_occupancy()
                    #a_bfactor = atom.get_bfactor()  #is replaced by rarity value in the 2nd-to-last column 
                    atom_line = f'ATOM\t {a_idx}\t {a_name}\t {res_3let}\t {chain_id}\t {res_idx}\t {a_coords}\t {a_occ}\t {codon_rarity_val}\t {a_name}\n'
                    new_pdb_contents.append(atom_line)
        new_pdb_contents.append(f'TER\t {a_idx+1}\t  \t {res_3let}\t {chain_id}\t {res_idx}\n') 
        new_pdb_contents.append('END') 
    return new_pdb_contents
    
pdbfile = Path(args.PDBID)
outfile = Path(args.outfile)
codon_table = load_codon_table(args.codontables, args.taxID)
#print(f'Codon fraction table for {args.taxID} is {codon_table}')
pdb_struct = pdbfile2struct(pdbfile)
if args.nucseq is not None: #Need to decide which overrides
    seq_record = SeqIO.read(args.nucseq, 'fasta') 
    nucseq = str(seq_record.seq).upper()
   #print(f'The nucleic acid sequence is: {nucseq}')
codons = nucseq2codons(nucseq)
#print(codons)
new_pdb_contents = replace_b_factor(pdb_struct, codons, codon_table) 
outfile.write_text(''.join(new_pdb_contents))
