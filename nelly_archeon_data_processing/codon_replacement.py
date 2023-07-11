
import argparse



def translate(args):
    with open(args.DNAsequence, "r") as f:
        seq = f.read()
        seq = seq.replace("\n", "")
        seq = seq.replace("\r", "")

    outputfile = args.outfolder + "/new_codon_sequence.txt"
    outputhere = open(outputfile, "w")

    with open(args.dictionary, encoding="utf8", errors='ignore') as d:
        dictionary = d.read().splitlines()[2:]
    bin214codonlist = []
    ecolicodonlist  = []

    for line in dictionary:
        col = line.split("\t")
        bin214codonlist.append( col[1] )
        ecolicodonlist.append(col[2])
    print(bin214codonlist)
    print(ecolicodonlist)

    table = {bin214codonlist[i]: ecolicodonlist[i] for i in range(len(bin214codonlist))}
    print(table)

    new_codon = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            new_codon += table[codon]
    print(new_codon)
    outputhere.write(new_codon)
    outputhere.close()

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-r","--DNAsequence", help="original DNA sequence ",required=True)
    parser.add_argument("-d","--dictionary", help="The codon dictionary of ecoli and new organism",required=True)
    parser.add_argument("-o","--outfolder", help="folder to save final codon sequence",required=True)

    args = parser.parse_args()
    return args


def main():
    args = parseargs()
    translate(args=args)



if __name__ == '__main__':
    main()
