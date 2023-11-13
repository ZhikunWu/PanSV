from BaseFunc import parse_genotype_format, Infor_target_values
import argparse
import sys
import os

#usage: python ~/github/PanSV/script/MergeVCFFilt.py --vcf /home/wuzhikun/Project/PanVigna/SVCall/Sniffles2/Samples_SV_merge.vcf --out filt_sv.vcf > temp.tx

def genotype_revised(Format, Genos):
    miss = 0
    fc = 0
    ac = 0
    AltCount = 0
    NewGenos = []
    for i in Genos:
        gt = parse_genotype_format(Format, i, "GT")
        dr = parse_genotype_format(Format, i, "DR")
        dv = parse_genotype_format(Format, i, "DV")
        dr = int(dr)
        dv = int(dv)
        if gt == "./.":
            miss += 1
        else:
            af = dv / (dr + dv)
            if af < 0.2:
                gt = "0/0"
                fc += 2
            elif af > 0.8:
                gt = "1/1"
                ac += 2
                AltCount += 1
            elif af >= 0.2 and af <= 0.8:
                gt = "0/1"
                fc += 1
                ac += 1
                AltCount += 1
            else:
                print("Please check the value of %s" % af)
                sys.exit(1)
        newgeno = "%s:%d:%d" %  (gt, dr, dv)
        NewGenos.append(newgeno)
    nomiss = len(Genos) - miss
    missRatio = miss / len(Genos)
    PopAF = ac / (fc + ac)
    return missRatio, PopAF, AltCount, "\t".join(NewGenos)




def population_vcf_filt(vcf_file, out_file):
    SVTYPES = ["INS", "DEL", "INV", "DUP"]

    in_h = open(vcf_file, "r")
    out_h = open(out_file, "w")
    for line in in_h:
        line = line.strip()
        lines = line.split("\t")
        if line.startswith("#"):
            out_h.write("%s\n" % line)
        else:
            Chr, Start = lines[:2]
            Infor = lines[7]
            Format = lines[8]
            Genos = lines[9:]
            missRatio, PopAF, AltCount, NewGenos = genotype_revised(Format, Genos)

            SVType = Infor_target_values(Infor, "SVTYPE")
            if SVType in SVTYPES:
                SVLen = Infor_target_values(Infor, "SVLEN")
                SVLen = int(SVLen.strip("-"))
                End = Infor_target_values(Infor, "END")
                if AltCount == 0:
                    print("No alternative genotype:%s\n" % line)
                else:
                    if SVLen >= 5000000 and AltCount <= 1:
                        print("Too long and rare: %s" % line)
                    else:
                        # if missRatio > 0.2:
                        #     print("Too many genotypes missed: %s" % line )
                        # else:
                        out_h.write("%s\t%s\t%s\n" % ("\t".join(lines[:8]), "GT:DR:DV", NewGenos))
            else:
                print("Not target type: %s" % line)
    in_h.close()
    out_h.close()



def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-v", "--vcf", help="The input vcf file.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    population_vcf_filt(args.vcf, args.out)

if __name__ == "__main__":
    main()

