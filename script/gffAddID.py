import argparse

#usage: python /home/wuzhikun/github/PanSV/script/gffAddID.py --gff CPG11_align.gff --out temp.gff


def target_value(target):
    KeyValues = {}
    tt = target.split(";")
    for t in tt:
        k, v = t.split("=")
        KeyValues[k] = v
    return KeyValues

def add_ID(record):
    featureCount = {}
    newRecord = []
    for r in record:
        lines = r.split("\t")
        feature = lines[2]
        Infor = lines[8]
        KeyValues = target_value(Infor)
        if feature == "mRNA":
            ID = KeyValues["ID"]
            newInfor = "ID=%s" % ID
        else:
            Parent = KeyValues["Parent"]
            if feature not in featureCount:
                feCount = 0
                featureCount[feature] = feCount
                ID = "%s.%s%d" % (Parent, feature, feCount)
                newInfor = "ID=%s;Parent=%s" % (ID, Parent)
            else:
                feCount = featureCount[feature] + 1
                featureCount[feature] = feCount
                ID = "%s.%s%d" % (Parent, feature, feCount)
                newInfor = "ID=%s;Parent=%s" % (ID, Parent)
        rr = "%s\t%s" % ("\t".join(lines[:8]), newInfor)
        newRecord.append(rr)
    return newRecord



def gff_add_ID(gff_file, out_file):
    gff_h = open(gff_file, "r")
    out_h = open(out_file, "w")
    RECORDS = []
    Record = []
    for line in gff_h:
        line = line.strip()
        if line.startswith("##gff"):
            out_h.write("%s\n" % line)
        elif line.startswith("##PAF"):
            RECORDS.append(Record)
            Record = []
        else:
            Record.append(line)
    RECORDS.append(Record)
    gff_h.close()

    for record in RECORDS:
        if record != []:
            newRecord = add_ID(record)
            out_h.write("%s\n" % "\n".join(newRecord))
    out_h.close()



def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-g", "--gff", help="The input gff file.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    gff_add_ID(args.gff, args.out)

if __name__ == "__main__":
    main()

