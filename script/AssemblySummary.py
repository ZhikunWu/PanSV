import argparse
import os
import collections

#usage: python ~/github/PanSV/script/AssemblySummary.py --file CPG01/03.ctg_graph/nd.asm.fasta.stat --out temp.txt

def assembly_summary(in_file, out_file):
    """
    in_file:
    Type           Length (bp)            Count (#)
    N10             24770484                   2
    N20             23973100                   4
    N30             21462069                   6
    N40             19328709                   9
    N50             16423765                  12
    N60             15270408                  15
    N70             12252099                  19
    N80              9404406                  23
    N90              4866499                  31

    Min.               32534                   -
    Max.            35179134                   -
    Ave.             4497194                   -
    Total          508182972                 113

    """
    Samples = []
    FileSummary = collections.defaultdict(dict)
    files = in_file.split(',')
    for f in files:
        name = f.strip().split('/')[-3]
        Samples.append(name)
        in_h = open(f, "r")
        for line in in_h:
            line =line.strip()
            lines = line.split()
            if line.startswith("N50"):
                FileSummary[name]["N50_length"] = lines[1]
                FileSummary[name]["N50_count"] = lines[2]
            elif line.startswith("Min"):
                FileSummary[name]["Min"] = lines[1]
            elif line.startswith("Max"):
                FileSummary[name]["Max"] = lines[1]
            elif line.startswith("Ave"):
                FileSummary[name]["Ave"] = lines[1]
            elif line.startswith("Total"):
                FileSummary[name]["Bases"] = lines[1]
                FileSummary[name]["Contigs"] = lines[2]
        in_h.close()

    out_h = open(out_file, "w")
    out_h.write("Contig_number\tTotal_length\tN50_count\tN50_length\tMin\tAve\tMax\n")
    for s in sorted(Samples):
        Contigs = FileSummary[s]["Contigs"]
        Bases = FileSummary[s]["Bases"]
        N50_count = FileSummary[s]["N50_count"]
        N50_length = FileSummary[s]["N50_length"]
        Min = FileSummary[s]["Min"]
        Ave = FileSummary[s]["Ave"]
        Max = FileSummary[s]["Max"]
        out_h.write("%s\t%s\n" % (s,"\t".join([Contigs, Bases, N50_count, N50_length, Min, Ave, Max])))
    out_h.close()





def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-f", "--file", help="The file of the genome stats.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    assembly_summary(args.file, args.out)

if __name__ == "__main__":
    main()
