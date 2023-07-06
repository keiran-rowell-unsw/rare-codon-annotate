import Bio
from Bio.Seq import Seq 
from pathlib import Path 
    
p = Path('.')

fna_files = list(p.glob('*.fna'))

for fname in fna_files:
        DNA_contents = fname.read_text()
        print(fname)
        DNA_seq = Seq(DNA_contents)
        AA_seq = DNA_seq.translate()
        out_fname = str(fname).replace('.fna','.faa')
        outf = Path(out_fname)
        outf.write_text(str(AA_seq))
