import argparse
import sys
import os
import collections

#usage: python ~/github/PanSV/script/snpEffStats.py --annotation Sample_jointcall.twoallel.maf005.filtmiss.recode.annotation_test.vcf --out temp.txt

def snpEff_anno_stats(in_file, out_file):
    FeatureNumber = collections.Counter()
    in_h = open(in_file, "r")
    out_h = open(out_file, "w")
    out_h.write("Feature\tNumber\n")
    for line in in_h:
        line = line.strip()
        if line.startswith("#"):
            continue
        else:
            lines = line.split("\t")
            Ref, Alt = lines[3:5]
            Infor = lines[7]
            Infors = Infor.split("|")
            feature = Infors[1]
            FeatureNumber[feature] += 1
    in_h.close()
    
    for c in FeatureNumber:
        n = FeatureNumber[c]
        out_h.write("%s\t%d\n" % (c, n))
    out_h.close()    



def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-a", "--annotation", help="The input annotation file.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    snpEff_anno_stats(args.annotation, args.out)

if __name__ == "__main__":
    main()
