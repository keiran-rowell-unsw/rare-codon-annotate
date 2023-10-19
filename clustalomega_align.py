from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline
from tqdm import tqdm

def perform_alignment(input_file, output_file):
    # Define the ClustalOmega command line
    clustalo_cline = ClustalOmegaCommandline(infile=input_file, outfile=output_file, verbose=True, auto=True)
    
    # Start the ClustalOmega alignment
    process = clustalo_cline()
    
    # Use tqdm to track the progress
    with tqdm(total=100, desc="Aligning", dynamic_ncols=True) as pbar:
        for i in range(101):
            pbar.update(1)
            pbar.refresh()
            process.poll()
            if process.returncode is not None:
                break
    
    # Load the alignment
    alignment = AlignIO.read(open(output_file), "clustal")
    
    # Print the alignment
    print(alignment)

# Prompt the user for the input and output file names
input_file = input("Enter the input file name: ")
output_file = input("Enter the output file name: ")

perform_alignment(input_file, output_file)
