from tinyfasta import FastaParser
import argparse

#usage: python ~/github/Centromere/script/AnchorRatio.py --fasta Zygaena_filipendulae.fasta --out temp.txt

def anchor_ratio(fasta_file, out_file):
    """
    out_fie:
    Species	Total	Chromosome	Scaffold	AnchorRatio
    Zygaena_filipendulae	365,946,273	356,282,259	9,664,014	0.974
    """
    name = fasta_file.split('/')[-1].strip("fasta").strip(".")
    print(name)

    TotalLen = 0
    ChromLen = 0
    ScaffoldLen = 0
    for record in FastaParser(fasta_file):
        desc = str(record.description)
        seq = str(record.sequence)
        seqLen = len(seq)
        if "chromosome" in desc.lower():
            ChromLen += seqLen
        else:
            ScaffoldLen += seqLen
        
    if ChromLen == 0:
        print("Please check whether the string 'chromosome' is in the assembly file.")

        TotalLen = 0
        ChromLen = 0
        ScaffoldLen = 0
        for record in FastaParser(fasta_file):
            desc = str(record.description)
            seq = str(record.sequence)
            seqLen = len(seq)
            if seqLen > 5000000:
                ChromLen += seqLen
            else:
                ScaffoldLen += seqLen


    TotalLen = ChromLen + ScaffoldLen
    AnchorRatio = ChromLen / TotalLen
    AnchorRatio = "%.3f" % AnchorRatio
    out_h = open(out_file, "w")
    out_h.write("Species\tTotal\tChromosome\tScaffold\tAnchorRatio\n")
    out_h.write("%s\t%s\t%s\t%s\t%s\n" % (name, format(TotalLen, ","), format(ChromLen, ","), format(ScaffoldLen, ","), AnchorRatio))


def main():
    parser = argparse.ArgumentParser(description='Get the anchored ratio of the assembly.')
    parser.add_argument('-f', '--fasta', help='The input fasta file.')
    parser.add_argument('-o', '--out', help='The output prefix name.')
    args = parser.parse_args()
    anchor_ratio(args.fasta, args.out)



if __name__ == '__main__':
    main()
