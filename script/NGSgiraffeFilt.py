from BaseFunc import parse_genotype_format, Infor_target_values
import argparse
import sys
import os

#usage: python ~/github/PanSV/script/NGSgiraffeFilt.py --vcf   /home/wuzhikun/Project/PanVigna/SVCall/NGS/CPG32.giraffe.vcf --out out.vcf --altThreshold 10

def NGS_giraffe_filt(vcf_file, out_file, AltThreshold):
    AltThreshold = int(AltThreshold)
    out_h = open(out_file, "w")
    vcf_h = open(vcf_file, "r")
    for line in vcf_h:
        line = line.strip()
        lines = line.split("\t")
        if line.startswith("#"):
            out_h.write("%s\n" % line)
        else:
            Format = lines[8]
            Geno = lines[9]
            GT = parse_genotype_format(Format, Geno, "GT")
            AD = parse_genotype_format(Format, Geno, "AD")
            if GT != "./." and GT != "0/0":
                ADS = AD.split(",")
                ADLen = len(ADS)
                if ADLen == 2:
                    RefDepth, AltDepth = ADS
                    RefDepth = int(RefDepth)
                    AltDepth = int(AltDepth)
                    Sum = RefDepth + AltDepth
                    if Sum != 0:
                        ratio = AltDepth / (RefDepth + AltDepth)
                        if AltDepth >= AltThreshold:
                            out_h.write("%s\n" % line)
                elif ADLen == 3:
                    RefDepth, AltDepth1, AltDepth2 = ADS
                    AltDepth1 = int(AltDepth1)
                    AltDepth2 = int(AltDepth2)
                    if AltDepth1 >= AltThreshold or AltDepth2 >= AltThreshold:
                        out_h.write("%s\n" % line)
                else:
                    print(line)
    vcf_h.close()
    out_h.close()


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-v", "--vcf", help="The input vcf file.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-a", "--altThreshold", default=8, help="The altnative depth threshold.")
    args = parser.parse_args()
    NGS_giraffe_filt(args.vcf, args.out, args.altThreshold)

if __name__ == "__main__":
    main()


