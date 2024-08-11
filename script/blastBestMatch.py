#!/usr/bin/python
import collections
import argparse
import sys
import os
import operator

#usage:  python ~/github/GenomeAssembly/src/GenomeAssembly/blastBestMatch.py --blast nr_diamond_blast.txt --out nr_diamond_blast_best.txt


def blast_best_match(blast_file, out_file):
    """
    balst_file:
    Vu11G000020	XP_027912833.1	100	179	0	0	1	179	1	179	2.43e-115	338
    Vu11G000020	XP_027912832.1	100	179	0	0	1	179	1	179	3.71e-115	338
    Vu11G000020	XP_017437079.1	99.4	179	1	0	1	179	1	179	1.98e-112	338
    """
    BlastScores = collections.defaultdict(dict)
    in_h = open(blast_file, "r")
    for line in in_h:
        line = line.strip()
        lines = line.split("\t")
        query, subject = lines[:2]
        score = float(lines[-1])
        BlastScores[query][subject] = score
    in_h.close()

    out_h = open(out_file, "w")
    for q in BlastScores:
        subjects = BlastScores[q]
        sortScores = sorted(subjects.items(), key=operator.itemgetter(1), reverse=True)
        targetSub = sortScores[0][0]
        out_h.write("%s\t%s\n" % (q, targetSub))
    out_h.close()



def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-b', '--blast', help='The input blast file.')
    parser.add_argument('-o', '--out', help='The output file.')
    args = parser.parse_args()
    blast_best_match(args.blast, args.out)


if __name__ == '__main__':
    main()
