from tinyfasta import FastaParser
import argparse


#usage: python ~/github/Centromere/script/AssemblyGap.py --genome /home/wuzhikun/Project/Centromere/genome/Chromosome/Yponomeuta_sedellus.fasta --out temp.txt

def identify_genome_gap(genome_file, out_file):
    ChrSeq = {}
    for record in FastaParser(genome_file):
        sequence = str(record.sequence)
        desc = str(record.description)
        desc = desc.strip(">").split()[0]
        ChrSeq[desc] = sequence

    gaps = []
    for c in ChrSeq:
        sequence = ChrSeq[c]
        gap_start = None
        for i, base in enumerate(sequence):
            if base.upper() == "N":
                if gap_start is None:
                    gap_start = i
            elif gap_start is not None:
                gap_end = i - 1
                gaps.append((c, gap_start, gap_end))
                gap_start = None
        if gap_start is not None:
            gap_end = len(sequence) - 1
            gaps.append((c, gap_start, gap_end))
    out_h = open(out_file, "w")
    for c, gap_start, gap_end in gaps:
        out_h.write(f"{c}\t{gap_start}\t{gap_end}\n")
    out_h.close()



def main():
    parser = argparse.ArgumentParser(description='Get the gap of the assembly')
    parser.add_argument("-g", "--genome", help="The file contain genome chromosome length.")
    parser.add_argument('-o', '--out', help='The output prefix name.')
    args = parser.parse_args()
    identify_genome_gap(args.genome, args.out)


if __name__ == '__main__':
    main()   
