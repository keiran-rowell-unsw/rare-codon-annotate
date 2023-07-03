import Bio
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1, seq3
import numpy as np

parser = PDBParser()

PDB_struct = parser.get_structure("test-PDB", "RCSB-PDB-nitrogenase-1n2c.pdb")
EBI_struct = parser.get_structure("test-EBI", "EBI-nitrogenase-AF-P00459-F1-model_v4.pdb")

for model in EBI_struct:
    a_idx = 0  # no neat way to link idx to enumerate 
    for chain in model:
        for res_idx, residue in enumerate(chain):
            for atom in residue:
                a_idx +=1
                a_name = atom.get_name()
                res_name = residue.get_resname()
                chain_id = chain.get_id().replace('<Chain id=','').replace('>','') 
                #a_coords = atom.get_coord()
                a_coords = np.array2string(atom.get_coord()).replace('[','').replace(']','')
                a_occ = atom.get_occupancy() 
                a_bfactor = atom.get_bfactor()
                new_bfactor = 'ZZZ' 
                atom_line = f'ATOM\t {a_idx}\t {a_name}\t {res_name}\t {chain_id}\t {res_idx}\t {a_coords}\t {a_occ}\t {new_bfactor}\t {a_name}'  
                print(atom_line)
#for idx, res in enumerate(EBI_struct.get_residues()):
#    print(seq1(res.get_resname()))
#    for a in res:
#        new_bfactor = 'ZZ'
#       print(atom_line)



