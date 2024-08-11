import argparse
from tinyfasta import FastaParser
import sys


def rename_value(Value, sample):
    newValue = []
    values = Value.split(";")
    for v in values:
        i, j = v.split("=")
        newV = i + "=" + sample + "_" + j
        newValue.append(newV)
    NewValues = ";".join(newValue)
    return NewValues




def add_sample_name(in_file, out_file, Format, sample):
    out_h = open(out_file, "w")
    if Format == "gff":
        in_h = open(in_file, "r")
        for line in in_h:
            line = line.strip()
            lines = line.split("\t")
            if line.startswith("#"):
                out_h.write("%s\n" % line)
            elif line == "":
                continue
            else:
                Value = lines[8]
                NewValues = rename_value(Value, sample)
                out_h.write("%s\t%s\n" % ("\t".join(lines[:8]), NewValues))
        out_h.close()
    elif Format == "fasta":
        for record in FastaParser(in_file):
            name = sample + "_Vu"
            desc = str(record.description).strip(">")
            newDesc = name + "_" + desc
            seq = str(record.sequence)
            out_h.write("%s\n%s\n" % (newDesc, seq))
        out_h.close()
    else:
        print("Please make sure that the format is fasta or gff")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-i", "--input", help="The input file.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-n", "--name", help="The name of sample.")
    parser.add_argument("-f", "--format", help="The file format.")
    args = parser.parse_args()
    add_sample_name(args.input, args.out, args.format, args.name)

if __name__ == "__main__":
    main()


