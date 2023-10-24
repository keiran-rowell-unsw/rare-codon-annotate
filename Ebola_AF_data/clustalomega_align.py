from Bio import AlignIO, SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline
from tqdm import tqdm
import subprocess

def perform_alignment(input_file, output_file):
    # Define the ClustalOmega command line
    clustalo_cline = ClustalOmegaCommandline(infile=input_file, outfile=output_file, verbose=True, auto=True)
    
    # Start the ClustalOmega alignment as a subprocess
    process = subprocess.Popen(str(clustalo_cline), shell=True)
    
    # Use tqdm to track the progress
    with tqdm(total=100, desc="Aligning", dynamic_ncols=True) as pbar:
        while process.poll() is None:
            pbar.update(1)
            pbar.refresh()
    
    # Convert the output to Clustal format
    temp_output_file = output_file + "_temp"
    with open(output_file) as input_handle:
        sequences = list(SeqIO.parse(input_handle, "fasta"))
    with open(temp_output_file, "w") as output_handle:
        SeqIO.write(sequences, output_handle, "clustal")
    
    # Load the alignment
    alignment = AlignIO.read(open(temp_output_file), "clustal")
    
    # Print the alignment
    print(alignment)

# Prompt the user for the input and output file names
input_file = input("Enter the input file name: ")
output_file = input("Enter the output file name: ")

perform_alignment(input_file, output_file)
