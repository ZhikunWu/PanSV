import argparse

#usage: python /home/wuzhikun/github/PanSV/script/miniProtIdenFilter.py --gff aln.gff --out aln_filter.gff

def target_value(Infor):
    KeyValues = {}
    infors = Infor.split(";")
    for i in infors:
        item, value = i.split("=")
        KeyValues[item] = value
    return KeyValues

def miniprot_identity_filter(gff_file, out_file, idenThreshold):
    idenThreshold = float(idenThreshold)
    in_h = open(gff_file, "r")
    out_h = open(out_file, "w")
    Remove = ""
    for line in in_h:
        line = line.strip()
        if not line.startswith("#"):
            lines = line.split("\t")
            feature = lines[2]
            Infor = lines[8]
            InforValues = target_value(Infor)
            if feature == "mRNA":
                ID = InforValues["ID"]
                ident = float(InforValues["Identity"])
                if ident <= idenThreshold:
                    Remove = ID
                else:
                    Remove = ""
                    out_h.write("%s\n" % line)
            elif feature == "CDS":
                Parent = InforValues["Parent"]
                if Parent != Remove:
                    out_h.write("%s\n" % line)
    in_h.close()
    out_h.close()
    



def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-g", "--gff", help="The input gff file.")
    parser.add_argument("-i", "--idenThreshold", default=0.9, help="The identity thrshold.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    miniprot_identity_filter(args.gff, args.out, args.idenThreshold)

if __name__ == "__main__":
    main()

