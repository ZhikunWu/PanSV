from BaseFunc import Infor_target_values
import argparse
import sys

#usage: python ~/github/PanSV/script/HomoSVRecord.py --vcf SPS1.sniffles.vcf --out temp.bed


def homo_SV_record(vcf_file, out_file, REThreshold, AFThreshold):
    REThreshold = int(REThreshold)
    AFThreshold = float(AFThreshold)
    in_h = open(vcf_file, "r")
    out_h = open(out_file, "w")
    for line in in_h:
        line = line.strip()
        if line.startswith("#") or line == "":
            continue
        else:
            lines = line.split("\t")
            Chr, Start = lines[:2]
            Ref, Alt = lines[3:5]
            Infor = lines[7]
            SVLen = Infor_target_values(Infor, "SVLEN")
            SVLen = int(SVLen.strip("-"))
            End = Infor_target_values(Infor, "END")
            SVType = Infor_target_values(Infor, "SVTYPE")
            RE = int(Infor_target_values(Infor, "RE"))
            AF = float(Infor_target_values(Infor, "AF"))
            if RE >= REThreshold and AF >= AFThreshold and "/" not in SVType and SVType != "TRA":
                if SVType == "DEL":
                    out_h.write("%s\t%s\t%s\t%s\t%s\n" % (Chr, Start, End, SVType, Ref))
                elif SVType == "INS":
                    out_h.write("%s\t%s\t%s\t%s\t%s\n" % (Chr, Start, End, SVType, Alt))
                else:
                    out_h.write("%s\t%s\t%s\t%s\n" % (Chr, Start, End, SVType))
    in_h.close()
    out_h.close()


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-v", "--vcf", help="The vcf file.")
    parser.add_argument("-a", "--AFThreshold", default=0.9, help="The threshold of AF.")
    parser.add_argument("-r", "--REThreshold", default=20, help="The threshold of SV supported reads.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    homo_SV_record(args.vcf, args.out, args.REThreshold, args.AFThreshold)

if __name__ == "__main__":
    main()





            
