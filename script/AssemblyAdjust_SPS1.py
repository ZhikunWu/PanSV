from tinyfasta import FastaParser
import argparse

#usage: python ~/github/PanSV/script/AssemblyAdjust_SPS1.py --fasta ragtag.scaffold.fasta --out ragtag.scaffold.revised.fasta


# def revComp(sequence):
#     complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'} 

#     reverse_complement = ""
#     sequence = sequence.upper()

#     for base in sequence:
#         reverse_complement = complement[base] + reverse_complement
#     return reverse_complement


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
        if desc == "Vu01_RagTag":
            seq1 = seq[:22300000]
            seq2 = seq[22300000:33200000]
            seq3 = seq[33200000:]
            ChrSeq[desc] = seq3
        else:
            ChrSeq[desc] = seq

    for i in ChrSeq:
        Seq = ChrSeq[i]
        if i == "Vu05_RagTag":
            Seq = Seq + split + revComp(seq1) 
        elif i == "Vu06_RagTag":
            Seq = revComp(seq2) + split + Seq
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
