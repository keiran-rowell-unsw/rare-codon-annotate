from Bio import AlignIO, SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline
from tqdm import tqdm
import subprocess
import os

def perform_alignment(input_file):
    # Ensure the input file has the ".fasta" extension
    if not input_file.endswith(".fasta"):
        input_file += ".fasta"
    
    # Create the output filenames by appending "_output" to the input filename
    output_file = os.path.splitext(input_file)[0] + "_output.aln"
    maf_output_file = os.path.splitext(input_file)[0] + "_output.maf"
    
    # Define the ClustalOmega command line with the --outfmt option
    clustalo_cline = ClustalOmegaCommandline(infile=input_file, outfile=output_file, outfmt="clu", verbose=True, auto=True)
    
    # Start the ClustalOmega alignment as a subprocess
    process = subprocess.Popen(str(clustalo_cline), shell=True)
    
    # Use tqdm to track the progress
    with tqdm(total=100, desc="Aligning", dynamic_ncols=True) as pbar:
        while process.poll() is None:
            pbar.update(1)
            pbar.refresh()
    
    # Load the Clustalomega .aln output
    clustal_alignment = AlignIO.read(open(output_file), "clustal")
    
    # Print the Clustalomega .aln alignment
    print("Clustalomega .aln alignment:")
    print(clustal_alignment)
    
    # Convert the output to MAF format
    with open(maf_output_file, "w") as output_handle:
        SeqIO.write(clustal_alignment, output_handle, "maf")
    
    # Load the MAF alignment
    maf_alignment = AlignIO.read(open(maf_output_file), "maf")
    
    # Print the MAF alignment
    print("\nMAF alignment:")
    print(maf_alignment)

# Prompt the user for the input file name (without .fasta extension)
input_file = input("Enter the input file name (without .fasta extension): ")

perform_alignment(input_file)
