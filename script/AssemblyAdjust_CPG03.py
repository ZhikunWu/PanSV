from tinyfasta import FastaParser
import argparse

#usage: python ~/github/PanSV/script/AssemblyAdjust_CPG03.py --fasta ragtag.scaffold.fasta --out ragtag.scaffold.revised.fasta


def revComp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    seq = seq.upper()
    bases = list(seq)
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    return bases


def adjust_assembly(fa_file, out_file):
    split = "N" * 100
    ChrSeq = {}
    out_h = open(out_file, "w")
    for record in FastaParser(fa_file):
        desc = str(record.description)
        desc = desc.strip(">")
        seq = str(record.sequence)
        if desc == "Vu03_RagTag":
            seq1 = seq[:33650000]
            seq2 = seq[33650000:55100000]
            seq3 = seq[55100000:]
            ChrSeq[desc] = seq1 + split +  seq3
        else:
            ChrSeq[desc] = seq

    for i in ChrSeq:
        Seq = ChrSeq[i]
        if i == "Vu11_RagTag":
            Seq = Seq + split + revComp(seq2)
        out_h.write(">%s\n%s\n" % (i, Seq))
    out_h.close()

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-f", "--fasta", help="The input fasta file.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    adjust_assembly(args.fasta, args.out)

if __name__ == "__main__":
    main()
