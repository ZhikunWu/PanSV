from BaseFunc import parse_genotype_format, Infor_target_values
import argparse
import sys
import os

#usage: python ~/github/PanSV/script/PanSVExtract.py --vcf Samples_SV_merge_vgformat.vcf --out temp.vcf --sample CPG03

def PanSV_extract(vcf_file, out_file, target):
    out_h = open(out_file, "w")
    in_h = open(vcf_file, "r")
    for line in in_h:
        line = line.strip()
        lines = line.split("\t")
        if line.startswith("##"):
            out_h.write("%s\n" % line)
        elif line.startswith("#"):
            samples = lines[9:]
            if target in samples:
                targetIndex = samples.index(target)
                out_h.write("%s\t%s\n" % ("\t".join(lines[:9]), target))
            else:
                print("Please check whether the sample %s is in vcf file %s." % (target, vcf_file))
                sys.exit(1)
        else:
            Format = lines[8]
            genos = lines[9:]
            targetGeno = genos[targetIndex]
            GT = parse_genotype_format(Format, targetGeno, "GT")
            DR = parse_genotype_format(Format, targetGeno, "DR")
            DV = parse_genotype_format(Format, targetGeno, "DV")
            if GT == "./." or GT == "0/0":
                continue
            else:
                out_h.write("%s\t%s\n" % ("\t".join(lines[:8]), targetGeno))
    in_h.close()
    out_h.close()

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-v", "--vcf", help="The input vcf file.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-s", "--sample", default=8, help="The target sample name.")
    args = parser.parse_args()
    PanSV_extract(args.vcf, args.out, args.sample)

if __name__ == "__main__":
    main()

