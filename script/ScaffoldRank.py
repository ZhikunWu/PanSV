import argparse
import operator

#usage: python ~/github/PanSV/script/ScaffoldRank.py --fai MH63RS3.fasta.fai --out temp.txt

def scaffold_rank(fai_file, out_file):
    """
    fai_file:
    GWHBCKY00000004	37319239	123444041	100	101
    GWHBCKY00000005	31307418	161136518	100	101
    GWHBCKY00000006	31921180	192757056	100	101

    out_file:
    GWHBCKY00000001	45027022	0.114
    GWHBCKY00000003	39893253	0.215
    GWHBCKY00000004	37319239	0.309
    GWHBCKY00000002	37301368	0.403
    """
    ChrLength = {}
    Sum = 0
    in_h = open(fai_file, "r")
    for line in in_h:
        lines = line.strip().split("\t")
        Chr, length = lines[:2]
        length = int(length)
        ChrLength[Chr] = length
        Sum += length
    in_h.close()

    ### sort the length of chromosome(from high to low)
    sortLength = sorted(ChrLength.items(), key=operator.itemgetter(1), reverse=True)

    out_h = open(out_file, "w")
    ratioSum = 0
    for i in sortLength:
        c, l = i
        ratio = l / Sum
        ratioSum += ratio
        ratioS = "%.3f" % ratioSum
        out_h.write(f"{c}\t{l}\t{ratioS}\n")
    out_h.close()

def main():
    parser = argparse.ArgumentParser(description="sort the length of chromosomes")
    parser.add_argument("-f", "--fai", help="The fai file of the genome.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    scaffold_rank(args.fai, args.out)

if __name__ == "__main__":
    main()
