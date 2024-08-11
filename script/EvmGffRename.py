
import argparse

#python /home/wuzhikun/github/GenomeAssembly/src/GenomeAssembly/EvmGffRename.py  --gff Boleracea.EVM.gff3 --out temp.txt --out2 temp2.txt


def evm_rename_gene(gff_file, out_file, rename_file):
    """
    gff_file:
    Chr7    EVM gene    20  936 .   +   .   ID=evm.TU.Chr7.1;Name=EVM%20prediction%20Chr7.1
    Chr7    EVM mRNA    20  936 .   +   .   ID=evm.model.Chr7.1;Parent=evm.TU.Chr7.1;Name=EVM%20prediction%20Chr7.1
    Chr7    EVM exon    20  172 .   +   .   ID=evm.model.Chr7.1.exon1;Parent=evm.model.Chr7.1
    Chr7    EVM CDS 20  172 .   +   0   ID=cds.evm.model.Chr7.1;Parent=evm.model.Chr7.1
    Chr7    EVM exon    264 542 .   +   .   ID=evm.model.Chr7.1.exon2;Parent=evm.model.Chr7.1
    Chr7    EVM CDS 264 542 .   +   0   ID=cds.evm.model.Chr7.1;Parent=evm.model.Chr7.1
    Chr7    EVM exon    682 936 .   +   .   ID=evm.model.Chr7.1.exon3;Parent=evm.model.Chr7.1
    Chr7    EVM CDS 682 936 .   +   0   ID=cds.evm.model.Chr7.1;Parent=evm.model.Chr7.1
    """
    GeneName = {}
    in_h = open(gff_file, "r")
    out_h1 = open(rename_file, "w")
    out_h2 = open(out_file, "w")
    for line in in_h:
        line = line.strip()
        if not line.startswith("#") and line != "":
            lines = line.split("\t") 
            feature = lines[2]
            Value = lines[8]
            if feature == "mRNA":
                values = Value.split(";")
                ID = values[0].split("=")[-1]
                ids = ID.split(".")
                newID = ids[2] + "g" + ids[3].zfill(5) + "0.1"
                GeneName[ID] = newID
                out_h1.write("%s\t%s\n" % (ID, newID))
                Parent = values[1].split("=")[-1]
                parents = Parent.split(".")
                newParent = parents[2] + "g" + parents[3].zfill(5) + "0"
                GeneName[Parent] = newParent
                out_h1.write("%s\t%s\n" % (Parent, newParent))
    in_h.close()
    out_h1.close()

    in_h = open(gff_file, "r")
    for line in in_h:
        line = line.strip()
        if line.startswith("#"):
            out_h2.write("%s\n" % line)
        elif line == "":
            continue
        else:
            lines = line.split("\t")
            feature = lines[2]
            if feature == "gene":
                ID = lines[8].split(";")[0].split("=")[1]
                gene = GeneName[ID]
                out_h2.write("%s\tID=%s\n" % ("\t".join(lines[:8]), gene))
            elif feature == "mRNA":
                values = lines[8].split(";")
                ID = values[0].split("=")[1]
                Parent = values[1].split("=")[1]
                newID = GeneName[ID]
                newParent = GeneName[Parent]
                out_h2.write("%s\tID=%s;Parent=%s\n" % ("\t".join(lines[:8]), newID, newParent))
            elif feature == "exon":
                values = lines[8].split(";")
                ID = values[0].split("=")[1]
                IDPre = ".".join(ID.split(".")[:-1])
                exon = ID.split(".")[-1]
                newID = GeneName[IDPre] + "." + exon
                Parent = values[1].split("=")[1]
                newParent = GeneName[Parent]
                out_h2.write("%s\tID=%s;Parent=%s\n" % ("\t".join(lines[:8]), newID, newParent))
            elif feature == "CDS":
                values = lines[8].split(";")
                ID = values[0].split("=")[1]
                IDlast = ".".join(ID.split(".")[1:])
                cds = ID.split(".")[0]
                newID = cds + "." + GeneName[IDlast]
                Parent = values[1].split("=")[1]
                newParent = GeneName[Parent]
                out_h2.write("%s\tID=%s;Parent=%s\n" % ("\t".join(lines[:8]), newID, newParent))
    in_h.close()
    out_h2.close()



def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-g', '--gff', help='The input gff file.')
    parser.add_argument('-o', '--out', help='The output file.')
    parser.add_argument('-o2', '--out2', help='The out file 2.')
    args = parser.parse_args()
    evm_rename_gene(args.gff, args.out, args.out2)



if __name__ == '__main__':
    main()


