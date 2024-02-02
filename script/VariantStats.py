import statistics
import collections
import argparse
import os
import sys

#usage: python ~/github/PanSV/script/VariantStats.py --vcf Sample_jointcall.twoallel.maf005.GT.vcf --stats stats.txt --heter heter.txt


def varaint_stats(vcf_file, stat_file, heter_file):
    SNVChrCount = collections.Counter()
    INDELChrCount = collections.Counter()
    SampleGenos = collections.defaultdict(list)
    in_h = open(vcf_file, "r")
    for line in in_h:
        line = line.strip()
        lines = line.split("\t")
        if line.startswith("##"):
            continue
        elif line.startswith("#CHR"):
            samples = lines[9:]
        else:
            Chr, Pos, Id, Ref, Alt = lines[:5]
            genotypes = lines[9:]

            ### sample genotype
            for s, g in zip(samples, genotypes):
                g = g.replace("|", "/")
                if g != "./." and g !="0/0":
                    SampleGenos[s].append(g)


            ### variant stats 
            if len(Ref) == len(Alt):
                if Chr.startswith("Chr"):
                    SNVChrCount[Chr] += 1
                else:
                    SNVChrCount["Other"] += 1
            else:
                if Chr.startswith("Chr"):
                    INDELChrCount[Chr] += 1
                else:
                    INDELChrCount["Other"] += 1   
    in_h.close()

    # SNVCount = []
    # INDELCount = []
    stat_h = open(stat_file, "w")
    stat_h.write("Chr\tSNP\tInDel\tSum\n")
    Chrs = sorted(list(SNVChrCount.keys()))
    for Chr in SNVChrCount:
        cc = SNVChrCount[Chr]
        ci = INDELChrCount[Chr]
        stat_h.write("%s\t%d\t%d\t%d\n" % (Chr, cc, ci, cc+ci))
        # SNVChrCount.append(cc)
        # INDELChrCount.append(ci)
    stat_h.close()
    

    heter_h = open(heter_file, "w")
    heter_h.write("Sample\tHeterogeneity\n")
    for s in samples:
        ge = SampleGenos
        heterRatio = heter_ratio(ge)
        heter_h.write("%s\t%s\n" % (s, heterRatio))
    heter_h.close()


def heter_ratio(genos):
    heterCount = 0
    for g in genos:
        if g == "1/0" or g == "0/1":
            heterCount += 1
    heterRatio = "%.4f" % (heterCount / 632552108 * 100)
    return heterRatio

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-v", "--vcf", help="The input vcf file.")
    parser.add_argument("-s", "--stats", help="The output stats file.")
    parser.add_argument("-e", "--heter", help="The output heter file.")
    args = parser.parse_args()
    varaint_stats(args.vcf,  args.stats, args.heter)

if __name__ == "__main__":
    main()



                