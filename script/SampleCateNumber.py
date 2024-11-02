
import collections
import argparse

def sample_category_stats(cat_file, ortho_file, out_file):
    """
    ortho_file:
    Orthogroup      CPG01.pep       CPG02.pep       CPG03.pep       CPG05.pep       CPG06.pep       CPG08.pep       CPG09.pep       CPG10>
    OG0000000       31      44      49      41      29      52      51      63      56      44      51      46      50      43      47   >
    OG0000001       27      37      46      38      33      42      49      38      30      41      35      41      34      35      42   >
    OG0000002       22      29      22      39      28      26      20      27      33      32      22      23      34      28      21  

    cat_file:
    Orthogroup      Category
    OG0000000       Core
    OG0000001       Core
    OG0000002       Core
    """
    GeneCategory = {}
    Categories = set()
    cat_h = open(cat_file, "r")
    header = cat_h.readline()
    for line in cat_h:
        lines = line.strip().split('\t')
        gene, category = lines
        GeneCategory[gene] = category
        Categories.add(category)
    cat_h.close()

    SampleCateCount = collections.defaultdict(lambda: collections.Counter())
    in_h = open(ortho_file, "r")
    out_h = open(out_file, "w")
    headers = in_h.readline().strip().split("\t")
    samples = headers[1:-1]
    for line in in_h:
        lines = line.strip().split("\t")
        gene = lines[0]
        cat = GeneCategory[gene]
        counts = lines[1:-1]
        for c, s in zip(counts, samples):
            c = int(c)
            if c != 0:
                SampleCateCount[s][cat] += 1
    in_h.close()

    Cates = sorted(list(Categories))
    out_h.write("Sample\t%s\n" % "\t".join(Cates))
    for sample in samples:
        catCount = []
        for c in Cates:
            if c in SampleCateCount[sample]:
                n = SampleCateCount[sample][c]
            else:
                n = 0
            catCount.append(n)
        catCount = [str(n) for n in catCount]
        out_h.write("%s\t%s\n" % (sample, "\t".join(catCount)))
    out_h.close()


            

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-c", "--category", help="The input category file.")
    parser.add_argument("-g", "--orthogroup", help="The input orthogroup file.")
    args = parser.parse_args()
    sample_category_stats(args.category, args.orthogroup,  args.out)

if __name__ == "__main__":
    main()


