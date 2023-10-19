from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline

def perform_alignment(input_file, output_file):
    # Define the ClustalW command line with the full path to the executable
    clustalw_cline = ClustalwCommandline("/usr/local/bin/clustalw2", infile=input_file, outfile=output_file)
    
    # Run the ClustalW alignment
    clustalw_cline()
    
    # Load the alignment
    alignment = AlignIO.read(open(output_file), "clustal")
    
    # Print the alignment
    print(alignment)

# Prompt the user for the input and output file names
input_file = input("Enter the input file name: ")
output_file = input("Enter the output file name: ")

perform_alignment(input_file, output_file)

# Print that the job has finished
print("Alignment job completed successfully.")
