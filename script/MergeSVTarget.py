from BaseFunc import Infor_target_values
from BaseFunc import parse_genotype_format
import operator
import argparse

#usage: python /home/wuzhikun/github/SCZSV/script/MergeSVTarget.py --vcf /home/wuzhikun/Project/Schizophrenia/SVCall3/Merge/CN028_merge.vcf --out temp.vcf

def target_geno(gts, drs, dvs, supRead):
    TargetRatio = {}
    for i in range(len(gts)):
        g = gts[i]
        if g != "./.":
            r = int(drs[i])
            v = int(dvs[i])
            total = r + v
            if total != 0 :
                ratio = v / total
                if ratio >= 0.2 and total >= 5 and v >= supRead:
                    TargetRatio[(g, r, v)] = ratio
    ks = list(TargetRatio.keys())
    if len(ks) == 0:
        target = ""
    elif len(ks) == 1:
        target = ks[0]
    else:
        sortedTargetRatio = sorted(TargetRatio.items(), key=operator.itemgetter(1), reverse=True)
       
        target = sortedTargetRatio[0][0]
    return target


def MergeSV_filter(vcf_file, out_file, supThreshold, supRead):
    supRead = int(supRead)
    supThreshold = int(supThreshold)
    targetTypes = ["DEL", "INS", "DUP", "INV"]
    #ChrList = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"]
    ChrList = ["Vu01", "Vu02", "Vu03", "Vu04", "Vu05", "Vu06", "Vu07", "Vu08", "Vu09", "Vu10", "Vu11"]
    in_h = open(vcf_file, "r")
    out_h = open(out_file, "w")
    for line in in_h:
        line = line.strip()
        lines = line.split("\t")
        if line.startswith("##"):
            out_h.write("%s\n" % line)
        elif line.startswith("#CHR"):
            sample = lines[9].lstrip("0_")
            out_h.write("%s\t%s\n" % ("\t".join(lines[:9]), sample))
        else:
            Chr = lines[0]
            Infor = lines[7]
            SVType = Infor_target_values(Infor, "SVTYPE")
            Supp = int(Infor_target_values(Infor, "SUPP"))
            if Chr in ChrList and  SVType in targetTypes and Supp >= supThreshold:
                Format = lines[8]
                Genos = lines[9:]
                gts = [parse_genotype_format(Format, g, "GT") for g in Genos]
                drs = [parse_genotype_format(Format, g, "DR") for g in Genos]
                dvs = [parse_genotype_format(Format, g, "DV") for g in Genos]
                target = target_geno(gts, drs, dvs, supRead)
                if target != "":
                    gt, dr, dv = target
                    if gt == "0/0":
                        gt = "0/1"
                    Target = "%s:%d:%d" % (gt, dr, dv)
                    out_h.write("%s\tGT:DR:DV\t%s\n" % ("\t".join(lines[:8]), Target))
    in_h.close()
    out_h.close()

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-v", "--vcf", help="The input vcf file.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-s", "--supThreshold", default=2, help="The threshold of support tool number.")
    parser.add_argument("-r", "--supRead", default=3, help="The threshold of support reads.")
    args = parser.parse_args()
    MergeSV_filter(args.vcf, args.out, args.supThreshold, args.supRead)

if __name__ == "__main__":
    main()



