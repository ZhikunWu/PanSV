import argparse
import collections
import sys


def mix_base(Ref, Alt):
    baseDict = {"A/G": "R", "G/A": "R", "C/T": "Y", "T/C": "Y", "A/C": "M", "C/A": "M", "G/T": "K", "T/G": "K", "G/C": "S", "C/G": "S", "A/T": "W", "T/A" : "W"}
    ge = "%s/%s" % (Ref, Alt)
    g = baseDict[ge]
    return g



def geno_to_base(vcf_file, out_file):
    in_h = open(vcf_file, "r")
    out_h = open(out_file, "w")
    for line in in_h:
        line = line.strip()
        lines = line.split("\t")
        if line.startswith("##"):
            continue
        elif line.startswith("#"):
            samples = lines[9:]
            out_h.write("Chr\tPos\tRef\t%s\n" % "\t".join(samples))
        else:
            Chr, Pos = lines[:2]
            Ref, Alt = lines[3:5]
            if len(Ref) == 1 and len(Alt) == 1:
                Ref = Ref.upper()
                Alt = Alt.upper()
                genos = lines[9:]
                newGeno = []
                for ge in genos:
                    ge = ge.replace("|", "/")
                    if ge == "0/0":
                        g = Ref
                    elif ge == "1/1":
                        g = Alt
                    elif ge == "./." or ge == ".":
                        g = "N"
                    elif ge == "0/1" or ge == '1/0':
                        g = mix_base(Ref, Alt)
                    else:
                        print("Please check the genotype %s in line %s" % (ge, line))
                        sys.exit(1)
                    newGeno.append(g)
                if len(newGeno) == len(samples):
                    out_h.write("%s\t%s\t%s\t%s\n" % (Chr, Pos, Ref, "\t".join(newGeno)))
                else:
                    print("please check the genotype number of line %s." % line)
                    sys.exit(1)
    in_h.close()
    out_h.close()



def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', '--vcf', help='The input vcf file.')
    parser.add_argument('-o', '--out', help='The out file.')
    args = parser.parse_args()
    geno_to_base(args.vcf, args.out)



if __name__ == '__main__':
    main()

