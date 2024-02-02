#!/usr/bin/python
import argparse
import os
import sys
import re

#usage: python ~/github/NanoHub/src/NanoHub/bamFlagSum.py --file temp --out temp.xls

def flag_stat_summary(fileStr, out_file):
    files = fileStr.split(",")
    files = [f.strip() for f in files]

    out_h = open(out_file, "w")
    out_h.write("Sample\tTotal_reads\tMapped_reads\tRatio(%)\n")

    for f in files:
        fileBase = f.split("/")[-1].split(".")[0].split("_")[0]
        in_h = open(f, "r")
        for line in in_h:
            line = line.strip()
            lines = line.split()
            if "total" in line:
                Total = lines[0]
            elif "mapped (" in line:
                Mapped = lines[0]
                match = re.findall("(\d+.\d+)%", line)
                Ratio = match[0]
        out_h.write("%s\t%s\t%s\t%s\n" % (fileBase, Total, Mapped, Ratio))
    out_h.close()


def main():
    parser = argparse.ArgumentParser(description="Summary for bam flag statistics.")
    parser.add_argument("-f", "--file", help="The input flag stats files.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    flag_stat_summary(args.file, args.out)



if __name__ == "__main__":
    main()


