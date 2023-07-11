import Bio
from Bio import Entrez 
from Bio import SeqIO 

seq_IDs = ['NZ_STEB01000014.1']
#seq_IDs = ['NC_008488.1','NC_000913.3','NZ_STEB01000004.1','NZ_SSZK01000033.1','NZ_STEB01000035.1','NZ_SSZK01000018.1','NZ_STEB01000014.1']

Entrez.email = "k.rowell@unsw.edu.au" 

for seq_ID in seq_IDs:
    print(seq_ID)
    RefSeq_handle = Entrez.efetch(db='nucleotide', id=seq_ID, rettype='fasta', retmode='text')
    RefSeq_record = SeqIO.read(RefSeq_handle, 'fasta')
    RefSeq_nuc_seq = str(RefSeq_record.seq)
    fname = seq_ID + '.fasta'
    with open(fname, 'w') as f:
        f.write(RefSeq_nuc_seq)
RefSeq_handle.close()
