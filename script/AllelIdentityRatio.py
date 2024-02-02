import argparse
import os
import sys
import collections

#usage: python ~/github/PanSV/script/AllelIdentityRatio.py --vcf Sample_jointcall.twoallel.maf005.GT_test.vcf --out temp.txt > temp

def allele_identity_ratio(vcf_file, out_file):
    PairCounts = collections.Counter()
    in_h = open(vcf_file, "r")
    out_h = open(out_file, "w")
    allNum = 0
    for line in in_h:
        line = line.strip()
        lines = line.split("\t")
        if line.startswith("##"):
            continue
        elif line.startswith("#CHR"):
            samples = lines[9:]
            samNum = len(samples)
        else:
            if line.startswith("Chr"):
                allNum += 1
                genos = lines[9:]
                genotypes = [g.replace("|", "/") for g in genos]
                for i in range(samNum-1):
                    for j in range(i+1, samNum):
                        s1 = samples[i]
                        s2 = samples[j]
                        g1 = genotypes[i]
                        g2 = genotypes[j]
                        if g1 == g2:
                            PairCounts[(s1, s2)] += 1
    in_h.close()
                    
    keys = sorted(list(PairCounts.keys()))
    for k in keys:
        kk = "\t".join(list(k))
        count = PairCounts[k]
        ratio = "%.3f" % (count / allNum)
        out_h.write("%s\t%s\n" % (kk, ratio))
    out_h.close()


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-v", "--vcf", help="The input vcf file.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    allele_identity_ratio(args.vcf, args.out)

if __name__ == "__main__":
    main()
