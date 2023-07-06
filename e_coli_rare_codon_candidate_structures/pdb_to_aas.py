import Bio
from Bio.PDB import PDBParser 
from Bio.SeqUtils import seq1, seq3 
from pathlib import Path

p = Path('.')

pdb_files = list(p.glob('*.pdb'))

parser = PDBParser()

for pdb_file in pdb_files:
    AA_string = []
    struct_name = str(pdb_file).replace('.pdb', '')
    print(struct_name)
    pdb_struct = parser.get_structure(struct_name, pdb_file)
    out_fname = str(pdb_file).replace('.pdb','.fna')
    for residue in pdb_struct.get_residues():
        tags = residue.id
        if tags[0] == ' ':  # Standard residues shouldn't be labelled with H for heteroatom and W for water etc 
            AA_string.append(seq1(residue.get_resname()))
    outf = Path(out_fname)
    outf.write_text("".join(AA_string)) 
