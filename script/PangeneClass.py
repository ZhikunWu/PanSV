import argparse
import collections
import sys


def pangene_class(count_file, out_file, cat_out):
    """
    count_file"
    Orthogroup      CPG01.pep       CPG02.pep       CPG03.pep       CPG05.pep       CPG06.pep       CPG08.pep       CPG09.pep       CPG10>
    OG0000000       31      44      49      41      29      52      51      63      56      44      51      46      50      43      47   >
    OG0000001       27      37      46      38      33      42      49      38      30      41      35      41      34      35      42   >
    OG0000002       22      29      22      39      28      26      20      27      33      32      22      23      34      28      21   
    """
    Genes = []
    GeneCount = collections.Counter()
    in_h = open(count_file, "r")
    headers = in_h.readline().strip().split("\t")
    sampleNum = len(headers[1:-1])
    core = sampleNum
    softcore = core * 0.9
    print(core)
    print(softcore)
    for line in in_h:
        lines = line.strip().split("\t")
        gene = lines[0]
        Genes.append(gene)
        samples = lines[1:-1]
        for c in samples:
            c = int(c)
            if c != 0:
                GeneCount[gene] += 1
    in_h.close()

    CategoryCount = collections.Counter()
    out_h = open(out_file, "w")
    out_h.write("Orthogroup\tCategory\n")
    for s in Genes:
        count = GeneCount[s]
        if count == core:
            cat = "Core"
        elif count > softcore and count < core:
            cat = "Softcore"
        elif count >= 2 and count < softcore:
            cat = "Dispensable"
        elif count == 1:
            cat = "Private"
        else:
            print("Please check the sample number %d of gene %s" % (count, s))
            sys.exit(1)
        out_h.write("%s\t%s\n" % (s, cat))

        CategoryCount[cat] += 1
    out_h.close()

    geneNum = len(Genes)
    cat_h = open(cat_out, "w")
    cat_h.write("Category\tNumber\tRatio\n")
    for i in CategoryCount:
        c = CategoryCount[i]
        ratio = c / geneNum
        ratio = "%.3f" % ratio
        cat_h.write("%s\t%d\t%s\n" % (i, c, ratio))
    cat_h.close()


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-i", "--input", help="The input orthogroup file.")
    parser.add_argument("-c", "--category", help="The out category file.")
    parser.add_argument("-s", "--stats", help="The output summary file.")
    args = parser.parse_args()
    pangene_class(args.input, args.category, args.stats)

if __name__ == "__main__":
    main()

