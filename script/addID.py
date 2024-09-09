from BaseFunc import Infor_target_values
import argparse


def add_ID(vcf_file, out_file, name):
    targetTypes = ["DEL", "INS", "DUP", "INV"]
    in_h = open(vcf_file, "r")
    out_h = open(out_file, "w")
    for line in in_h:
        line = line.strip()
        if line.startswith("##"):
            out_h.write("%s\n" % line)
        elif line.startswith("#CHR"):
            lines = line.split("\t")
            sample = lines[-1].split("/")[-1].strip(".bam")
            out_h.write("%s\t%s\n" % ("\t".join(lines[:-1]), sample))
        else:
            lines = line.split("\t")
            ID = lines[2]
            Infor = lines[7]
            SVType = Infor_target_values(Infor, "SVTYPE")
            if SVType in targetTypes:
                if name != "None":
                    ID = "%s.%s.%s" % (name, SVType, ID)
                    out_h.write("%s\t%s\t%s\n" % ("\t".join(lines[:2]), ID, "\t".join(lines[3:])))
                else:
                    out_h.write("%s\n" % line)
    in_h.close()
    out_h.close()


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-v", "--vcf", help="The input vcf file.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-n", "--name", help="The name of tool.")
    args = parser.parse_args()
    add_ID(args.vcf, args.out, args.name)

if __name__ == "__main__":
    main()

