nelly_codon_table = {'*': {'TAA': 0.01, 'TAG': 0.0, 'TGA': 0.0}, 'A': {'GCA': 0.35, 'GCC': 0.18, 'GCG': 0.13, 'GCT': 0.33}, 'C': {'TGC': 0.23, 'TGT': 0.77}, 'D': {'GAC': 0.12, 'GAT': 0.88}, 'E': {'GAA': 0.85, 'GAG': 0.15}, 'F': {'TTC': 0.23, 'TTT': 0.77}, 'G': {'GGA': 0.44, 'GGC': 0.1, 'GGG': 0.2, 'GGT': 0.27}, 'H': {'CAC': 0.19, 'CAT': 0.81}, 'I': {'ATA': 0.22, 'ATC': 0.22, 'ATT': 0.56}, 'K': {'AAA': 0.82, 'AAG': 0.18}, 'L': {'CTA': 0.09, 'CTC': 0.08, 'CTG': 0.04, 'CTT': 0.15, 'TTA': 0.48, 'TTG': 0.15}, 'M': {'ATG': 1.0}, 'N': {'AAC': 0.17, 'AAT': 0.83}, 'P': {'CCA': 0.29, 'CCC': 0.22, 'CCG': 0.14, 'CCT': 0.34}, 'Q': {'CAA': 0.87, 'CAG': 0.13}, 'R': {'AGA': 0.18, 'AGG': 0.05, 'CGA': 0.31, 'CGC': 0.11, 'CGG': 0.16, 'CGT': 0.18}, 'S': {'AGC': 0.06, 'AGT': 0.2, 'TCA': 0.23, 'TCC': 0.13, 'TCG': 0.1, 'TCT': 0.27}, 'T': {'ACA': 0.28, 'ACC': 0.24, 'ACG': 0.11, 'ACT': 0.36}, 'V': {'GTA': 0.31, 'GTC': 0.15, 'GTG': 0.21, 'GTT': 0.34}, 'W': {'TGG': 1.0}, 'Y': {'TAC': 0.14, 'TAT': 0.86}}

#TaxID 83333
e_coli_codon_table = {'*': {'TAA': 0.64, 'TAG': 0.0, 'TGA': 0.36}, 'A': {'GCA': 0.21, 'GCC': 0.31, 'GCG': 0.38, 'GCT': 0.11}, 'C': {'TGC': 0.58, 'TGT': 0.42}, 'D': {'GAC': 0.35, 'GAT': 0.65}, 'E': {'GAA': 0.7, 'GAG': 0.3}, 'F': {'TTC': 0.43, 'TTT': 0.57}, 'G': {'GGA': 0.13, 'GGC': 0.46, 'GGG': 0.12, 'GGT': 0.29}, 'H': {'CAC': 0.45, 'CAT': 0.55}, 'I': {'ATA': 0.07, 'ATC': 0.35, 'ATT': 0.58}, 'K': {'AAA': 0.73, 'AAG': 0.27}, 'L': {'CTA': 0.05, 'CTC': 0.1, 'CTG': 0.46, 'CTT': 0.12, 'TTA': 0.15, 'TTG': 0.12}, 'M': {'ATG': 1.0}, 'N': {'AAC': 0.53, 'AAT': 0.47}, 'P': {'CCA': 0.14, 'CCC': 0.13, 'CCG': 0.55, 'CCT': 0.17}, 'Q': {'CAA': 0.3, 'CAG': 0.7}, 'R': {'AGA': 0.02, 'AGG': 0.03, 'CGA': 0.07, 'CGC': 0.44, 'CGG': 0.07, 'CGT': 0.36}, 'S': {'AGC': 0.33, 'AGT': 0.14, 'TCA': 0.15, 'TCC': 0.11, 'TCG': 0.16, 'TCT': 0.11}, 'T': {'ACA': 0.13, 'ACC': 0.47, 'ACG': 0.24, 'ACT': 0.16}, 'V': {'GTA': 0.17, 'GTC': 0.18, 'GTG': 0.4, 'GTT': 0.25}, 'W': {'TGG': 1.0}, 'Y': {'TAC': 0.47, 'TAT': 0.53}}

AA_alphabet = ['*','A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

input_codon_table = nelly_codon_table
output_codon_table = e_coli_codon_table

codon_changes = []
codon_changes_w_frac = []
for AA in AA_alphabet:
    inp_codons_sorted = {codon: frac for codon, frac in sorted(input_codon_table[AA].items(), key=lambda item: item[1])}
    out_codons_sorted = {codon: frac for codon, frac in sorted(output_codon_table[AA].items(), key=lambda item: item[1])}
    i_cod_sort_l = list(inp_codons_sorted) 
    o_cod_sort_l = list(out_codons_sorted) 
    #print(inp_codons_sorted)
    #print(i_cod_sort_l)
    #print(out_codons_sorted)
    #print(o_cod_sort_l)
     
    if(len(i_cod_sort_l) != len(o_cod_sort_l)):
        print(f'At amino acid {AA} the number of codons do not match between the organisms')
    for idx in range(len(inp_codons_sorted)):  #Python > 3.7: regular dictionaries are ordered
        if(i_cod_sort_l[idx] != o_cod_sort_l[idx]):
            codon_changes.append(f'{AA}: {i_cod_sort_l[idx]} -> {o_cod_sort_l[idx]}')  
            codon_changes_w_frac.append(f'{AA}: {i_cod_sort_l[idx]} ({input_codon_table[AA][i_cod_sort_l[idx]]}) -> {o_cod_sort_l[idx]} ({output_codon_table[AA][o_cod_sort_l[idx]]})')  # Long, unopt, line -- but will do   

#print(*codon_changes, sep='\n') 
prin(*codon_changes_w_frac, sep='\n') 
