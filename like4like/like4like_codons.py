import argparse
import pandas as pd
import pickle 

desc_string = '''Matches codon tables of two species in terms of the ordering within each species of fraction of the time codon encodes for an  amino acid.'''
epilog_string = '''This is *not* codon optimisation, but instead an attempt to match like4like in terms of codon rarity.

Files generated:
  - ..._codons_table.tsv: Inverse codon table of the two organism listed in side-by-side columns.
  _ ..._all_changes.tsv: Lists every occurrence of codons changing in rarity order between species.
  - ..._critical_changes_X%.tsv: Lists only codon changes where fractions differ by some critical percentage (specified with -c)'''

parser = argparse.ArgumentParser(prog='like4like_codons', description=desc_string, epilog=epilog_string, formatter_class=argparse.RawDescriptionHelpFormatter) #Don't line wrap the description
parser.add_argument('-i', '--inputtaxID', type=str, help='Taxonomy ID number of the input species', required=True)
parser.add_argument('-o', '--outputtaxID', type=str, help='Taxonomy ID number of the output species', required=True)
parser.add_argument('-c', '--criticalpercent', type=int, default=25, help='Threshold for the percentage codon use fractions can vary between species in the codon change assignment before it is listed in ..._critical_changes.tsvea')
args = parser.parse_args()

AA_alphabet = ['*','A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

with open('codon_tables.pkl', 'rb') as pickle_file:
    codon_tables = pickle.load(pickle_file)

in_codon_table = codon_tables[args.inputtaxID]
out_codon_table = codon_tables[args.outputtaxID]

print(f'The input codon table is {in_codon_table}')
print(f'The output codon table is {out_codon_table}')

in_org_df = pd.DataFrame(columns=['Amino acid', args.inputtaxID])    
out_org_df = pd.DataFrame(columns=['Amino acid', args.outputtaxID])    

codon_changes = []
codon_changes_w_frac = []

for AA in AA_alphabet:
    inp_codons_sorted = {codon: frac for codon, frac in sorted(in_codon_table[AA].items(), key=lambda item: item[1])}
    out_codons_sorted = {codon: frac for codon, frac in sorted(out_codon_table[AA].items(), key=lambda item: item[1])}
    i_cod_sort_l = list(inp_codons_sorted) 
    o_cod_sort_l = list(out_codons_sorted) 
     
    if(len(i_cod_sort_l) != len(o_cod_sort_l)): #Check number of codons per AA don't match. 
        print(f'At amino acid {AA} the number of codons do not match between the organisms')

    for idx in range(len(inp_codons_sorted)):  #Python > 3.7: regular dictionaries are ordered
        # Add the rows to the dfs with concat. Too few rows to care about performance, as of pandas 2.0 .append() is depreceated
        i_row_df = pd.DataFrame({'Amino acid': AA, args.inputtaxID: i_cod_sort_l[idx]}, index=[0])
        o_row_df = pd.DataFrame({'Amino acid': AA, args.outputtaxID: o_cod_sort_l[idx]}, index=[0])
        in_org_df = pd.concat([in_org_df, i_row_df], ignore_index=True)  
        out_org_df = pd.concat([out_org_df, o_row_df], ignore_index=True)  
        
        # Record any codon frequency ordering changes into a list 
        if(i_cod_sort_l[idx] != o_cod_sort_l[idx]):
            codon_changes.append(f'{AA}: {i_cod_sort_l[idx]} -> {o_cod_sort_l[idx]}')  
            i_codon_frac = in_codon_table[AA][i_cod_sort_l[idx]] 
            o_codon_frac = out_codon_table[AA][o_cod_sort_l[idx]]
            if float(i_codon_frac) == 0:
                percent_diff = 0  # TO DO: better divide by zero handling 
            else:  
                percent_diff = ( abs(float(i_codon_frac) - float(o_codon_frac)) / float(i_codon_frac)  ) * 100 
            if (round(percent_diff) > args.criticalpercent):
                codon_changes_w_frac.append(f'{AA}: {i_cod_sort_l[idx]} ({in_codon_table[AA][i_cod_sort_l[idx]]}) -> {o_cod_sort_l[idx]} ({out_codon_table[AA][o_cod_sort_l[idx]]})')  # Long, unopt, line -- but will do   

#Write out the files, like4like most used, but the changes could be useful
codon_tables_comp = pd.merge(in_org_df, out_org_df[args.outputtaxID], left_index=True, right_index=True, how='outer')
codon_tables_comp.to_csv(f'{args.inputtaxID}_{args.outputtaxID}_like4like_codons_table.tsv', sep='\t', index=False)

with open(f'{args.inputtaxID}_{args.outputtaxID}_like4like_codons_all_changes.tsv', 'w') as changes_file:
    for change in codon_changes: #change this if you want the fractions or not 
        changes_file.write(change + '\n')

with open(f'{args.inputtaxID}_{args.outputtaxID}_like4like_codons_critical_changes_{args.criticalpercent}%.tsv', 'w') as crit_changes_file:
    for change in codon_changes_w_frac: #change this if you want the fractions or not 
        crit_changes_file.write(change + '\n')
