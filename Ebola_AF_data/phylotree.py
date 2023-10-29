from Bio import AlignIO, SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo

def construct_tree(input_file, output_file):
    # Load the Clustalomega .aln output
    clustal_alignment = AlignIO.read(input_file, "clustal")
    
    # Convert the output to MAF format
    maf_output_file = output_file + ".maf"
    with open(maf_output_file, "w") as output_handle:
        SeqIO.write(clustal_alignment, output_handle, "maf")
    
    # Load the MAF alignment
    maf_alignment = AlignIO.read(maf_output_file, "maf")
    
    # Calculate distances
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(maf_alignment)
    
    # Construct the tree
    constructor = DistanceTreeConstructor(calculator, "upgma")
    tree = constructor.build_tree(maf_alignment)
    
    # Save the tree in Newick format
    tree_output_file = output_file + ".mod"
    Phylo.write(tree, tree_output_file, "newick")
    print("\nPhylogenetic tree saved in", tree_output_file)

# Prompt the user for the input and output file names
input_file = input("Enter the input file name: ")
output_file = input("Enter the output file name: ")

construct_tree(input_file, output_file)
