import argparse
import sys

#usage : python ~/github/PanSV/script/targetPheno.py  --phenotype /home/wuzhikun/Project/PanSV/phenotypes.txt --genome /home/wuzhikun/Project/PanSV/NGS/Variant/Sample_jointcall.twoallel.maf005.fam --out Sample_jointcall.twoallel.maf005.phenotype.txt

def target_pheno(pheno_file, geno_file, out_file):
    # GenoPheno = {}
    # cor_h = open(cor_file, "r")
    # for line in cor_h:
    #     lines = line.strip().split("\t")
    #     geno, pheno = lines
    #     GenoPheno[geno] = pheno
    # cor_h.close()

    PhenoValues = {}
    pheno_h = open(pheno_file, "r")
    header = pheno_h.readline().strip()
    for line in pheno_h:
        line = line.strip()
        lines = line.split("\t")
        Values = []
        for i in lines[1:]:
            if i == "" or i == "-":
                v = "NA"
            else:
                v = i
            Values.append(v)
        ID = lines[0]
        PhenoValues[ID] = "\t".join(Values)
    pheno_h.close()

    geno_h = open(geno_file, "r")
    out_h = open(out_file, "w")
    out_h.write("%s\n" % header)
    headers = header.split("\t")[1:]
    for line in geno_h:
        lines = line.strip().split()
        geno = lines[0]
        genoID = geno.split("_")[0]
        if genoID in PhenoValues:
            Values = PhenoValues[genoID]
            out_h.write("%s\t%s\n" % (geno, Values))
        else:
            print("Please check the geno ID %s" % genoID)
            phenoLen = len(headers)
            target = ["NA" for i in range(phenoLen)]
            out_h.write("%s\t%s\n" % (geno, "\t".join(target)))
            # sys.exit(1)
    geno_h.close()
    out_h.close()


def main():
    parser = argparse.ArgumentParser(description='')
    # parser.add_argument('-c', '--correspond', help='The input corresponding file.')
    parser.add_argument('-g', '--genome', help='The input genome file.')
    parser.add_argument('-p', '--phenotype', help='The input phenotype file.')
    parser.add_argument('-o', '--out', help='The output file.')
    args = parser.parse_args()
    target_pheno( args.phenotype, args.genome, args.out)


if __name__ == '__main__':
    main()


