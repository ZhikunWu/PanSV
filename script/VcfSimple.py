
from BaseFunc import Infor_target_values
import argparse

def vcf_simple(vcf_file, out_file):
    in_h = open(vcf_file, "r")
    out_h = open(out_file, "w")
    for line in in_h:
        line = line.strip()
        lines = line.split("\t")
        if line.startswith("##"):
            out_h.write("%s\n" % line)
        elif line.startswith("#CHR"):
            samples = lines[9:]
            newSamples = [i.split("_")[-1] for i in samples]
            out_h.write("%s\t%s\n" % ("\t".join(lines[:9]), "\t".join(newSamples)))
        else:
            Infor = lines[7]
            svtype  = Infor_target_values(Infor, "SVTYPE")
            svlen = Infor_target_values(Infor, "SVLEN")
            end = Infor_target_values(Infor, "END")
            newInfor = "SVTYPE=%s;SVLEN=%s;END=%s" % (svtype, svlen, end)
            genos = lines[9:]
            newGenos = [i.split(":")[0] for i in genos]
            out_h.write("%s\t%s\tGT\t%s\n" % ("\t".join(lines[:7]), newInfor, "\t".join(newGenos)))
    in_h.close()
    out_h.close()


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-v", "--vcf", help="The input gff file.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    vcf_simple(args.vcf, args.out)

if __name__ == "__main__":
    main()

