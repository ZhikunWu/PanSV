from tinyfasta import FastaParser
import argparse
import sys

#usage: python ~/github/PanSV/script/RenameScaffold.py --assembly /home/wuzhikun/Project/PanVigna/Assembly/RagTag/SPS1/ragtag.scaffold.fasta --out temp.fasta

def rename_scaffold(fa_file, out_file, lengthThreshold):
    lengthThreshold = int(lengthThreshold)

    samplePath = fa_file.split("/")
    if samplePath[-1] == "ragtag.scaffold.fasta":
        sample = samplePath[-2]
    else:
        sample = samplePath[-1].split(".")[0]

    out_h = open(out_file, "w")
    for record in FastaParser(fa_file):
        desc = str(record.description)
        desc = desc.strip(">")
        seq = str(record.sequence)
        seqLen = len(seq)
        if seqLen > lengthThreshold:
            Chr = desc.split("_")[0]
            newChr = "%s_%s" % (sample, Chr)
            out_h.write(">%s\n%s\n" % (newChr, seq))
    out_h.close()



def main():
    parser = argparse.ArgumentParser(description="Get aligned values for seq of BAM files.")
    parser.add_argument("-a", "--assembly", help="The input assembly file.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-t", "--lengthThreshold", default=5000000, help="The length threshold.")
    args = parser.parse_args()
    rename_scaffold(args.assembly, args.out, args.lengthThreshold)


if __name__ == "__main__":
    main()
