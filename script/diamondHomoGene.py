
from tinyfasta import FastaParser
import argparse

def best_homo(prot_file, blast_file, out_file, out_pro, idenThreshold):
    idenThreshold = float(idenThreshold)
    ProDic = {}
    for record in FastaParser(prot_file):
        desc = str(record.description).strip(">").split()[0]
        seq = str(record.sequence)
        ProDic[desc] = seq

    ProcTarget = {}
    ProcValue = {}
    in_h = open(blast_file, "r")
    for line in in_h:
        lines = line.strip().split("\t")
        prot, subject, ident = lines[:3]
        ident = float(ident)
        if ident > idenThreshold:
            if prot not in ProcValue:
                ProcValue[prot] = ident
                ProcTarget[prot] = subject
            else:
                i = ProcValue[prot]
                if ident > i:
                    ProcValue[prot] = ident
                    ProcTarget[prot] = subject
    in_h.close()

    out_h = open(out_file, "w")
    for i in ProcTarget:
        t = ProcTarget[i]
        out_h.write("%s\t%s\n" % (i, t))
    out_h.close()

    prot_h = open(out_pro, "w")
    for j in ProDic:
        if j not in ProcTarget:
            seq = ProDic[j]
            prot_h.write(">%s\n%s\n" % (j, seq))
    prot_h.close()



def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-p", "--protein", help="The input protein file with fasta format.")
    parser.add_argument("-b", "--blast", help="The blast file.")
    parser.add_argument("-ho", "--homologous", help="The homologous pair.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-d", "--idenThreshold", default=90, help="The threshold of identity")
    args = parser.parse_args()
    best_homo(args.protein, args.blast, args.homologous, args.out,  args.idenThreshold)

if __name__ == "__main__":
    main()





