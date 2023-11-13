#!/usr/bin/python
from BaseFunc import parse_genotype_format, Infor_target_values
import argparse
import sys
import os
import collections


#usage: python /home/wuzhikun/github/NanoHub/src/NanoHub/PopGenoStats.py --vcf /home/wuzhikun/Project/HuaXi2021/Schizophrenia/AllSamples/Samples_SV_refine_case_filt.vcf --stat1 /home/wuzhikun/Project/HuaXi2021/Schizophrenia/AllSamples/Samples_SV_refine_case_filt_stat1.txt --stat2 /home/wuzhikun/Project/HuaXi2021/Schizophrenia/AllSamples/Samples_SV_refine_case_filt_stat2.txt


def population_genotype_stats(vcf_file, out_file, stat_file):
    SampleTypeCount = collections.defaultdict(lambda: collections.Counter())
    in_h = open(vcf_file, "r")
    stat_h = open(stat_file, "w")
    Types = set()
    for line in in_h:
        line = line.strip()
        lines = line.split("\t")
        if line.startswith("##"):
            continue
        elif line.startswith("#"):
            samples = lines[9:]
        else:
            Chr, Start = lines[:2]
            genotypes = lines[9:]
            Infor = lines[7]
            svtype = Infor_target_values(Infor, "SVTYPE")
            End = Infor_target_values(Infor, "END")

            if "SVLEN" in Infor:
                SVLength = Infor_target_values(Infor, "SVLEN")
            elif "AVGLEN" in Infor:
                SVLength = Infor_target_values(Infor, "AVGLEN")
            else:
                print("Please check whether the tag name 'SVLEN' or 'AVGLEN' in the record %s." % line)
                sys.exit(1)
            # if svtype == "INV":
            #     SVLength = int(End) - int(Start)


            Region = "\t".join([Chr, Start, End])
            Tag = "_".join([Chr, Start, End, str(SVLength), svtype])



            Types.add(svtype)
            
            SVCount = 0
            Samples = []
            Format = lines[8]
            if len(samples) == len(genotypes):
                for s, gg in zip(samples, genotypes):
                    g = parse_genotype_format(Format, gg, "GT")
                    if g == "0/1" or g == "1/0" or g == "1/1":
                        SampleTypeCount[s][svtype] += 1
                        SVCount += 1
                        Samples.append(s)
                stat_h.write("%s\t%s\t%d\t%s\n" % (Region, Tag, SVCount, ",".join(Samples)))
            else:
                print("Please make sure that the number of samples and genotypes is identical.")
                sys.exit(1)
    in_h.close()

    out_h = open(out_file, "w")
    typeList = sorted(list(Types))
    out_h.write("Sample\t%s\tAll\n" % "\t".join(typeList))

    sampleList = sorted(list(SampleTypeCount.keys()))
    for s in sampleList:
        count = 0
        typeCount = []
        for t in  typeList:
            if t in SampleTypeCount[s]:
                c = SampleTypeCount[s][t]
            else:
                c = 0
            count += c
            typeCount.append(c)
        countStr = [str(i) for i in typeCount]
        out_h.write("%s\t%s\t%d\n" % (s, "\t".join(countStr), count))
    out_h.close()

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-v", "--vcf", help="The input vcf file.")
    parser.add_argument("-s1", "--stat1", help="The output file.")
    parser.add_argument("-s2", "--stat2", help="The output file.")
    args = parser.parse_args()
    population_genotype_stats(args.vcf,  args.stat1, args.stat2)

if __name__ == "__main__":
    main()


