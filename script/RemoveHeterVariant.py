
import argparse


def remove_heter_variant(vcf_file, out_file, heterThreshold):
    heterThreshold = float(heterThreshold)
    in_h = open(vcf_file, "r")
    out_h = open(out_file, "w")
    for line in in_h:
        line = line.strip()
        lines = line.split("\t")
        if line.startswith("##"):
            out_h.write("%s\n" % line)
        elif line.startswith("#CH"):
            out_h.write("%s\n" % line)
            samples = lines[9:]
            sampleNum = len(samples)
        else:
            genos = lines[9:]
            genos = [g.replace("|", "/") for g in genos]
            genoCount = 0
            heter = 0
            for g in genos:
                if g != "./.":
                    genoCount += 1
                    if g == "0/1" or g == "1/0":
                        heter += 1
            heterRatio = heter / genoCount
            if heterRatio <= heterThreshold:
                out_h.write("%s\n" % line)
    in_h.close()
    out_h.close()


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-v", "--vcf", help="The input gff file.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-t", "--heterThreshold", default=0.2, help="The threshold of heterzygosity")
    args = parser.parse_args()
    remove_heter_variant(args.vcf, args.out, args.heterThreshold)

if __name__ == "__main__":
    main()


