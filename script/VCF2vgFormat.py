from BaseFunc import Infor_target_values
# from BaseFunc import Infor_substution_tag
from tinyfasta import FastaParser
import argparse

#usage: python ~/github/PanSV/script/VCF2vgFormat.py --vcf /home/wuzhikun/Project/PanVigna/SVCall/Sniffles2/Samples_SV_merge.vcf --fasta /home/wuzhikun/Project/Vigna/Final/Final/Vigna_unguiculata_assembly.fasta --out /home/wuzhikun/Project/PanVigna/SVCall/Sniffles2/Samples_SV_merge_filt.vcf


def Infor_substution_tag(Infor, targetTag, NewValue):
    TagValues = {}
    Infors = Infor.split(";")

    Tags = []
    for i in Infors[1:]:
        print(i)
        tag, value = i.split("=")
        Tags.append(tag)
        TagValues[tag] = value

    if targetTag in Tags:
        TagValues[targetTag] = NewValue
    else:
        print("Please check whether the target tag %s is in the information record %s." % (targetTag, Infor))
        sys.exit(1)

    Record = []
    Record.append(Infors[0])
    for t in Tags:
        v = TagValues[t]
        record = t + "=" + str(v)
        Record.append(record)
    newInfor = ";".join(Record)
    return newInfor


def genome_seq(fa_file):
    ChrSeq = {}
    for record in FastaParser(fa_file):
        desc = str(record.description)
        desc = desc.strip(">")
        seq = str(record.sequence)
        ChrSeq[desc] = seq
    return ChrSeq


def revComp(sequence):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'} 

    reverse_complement = ""

    for base in sequence:
        reverse_complement = complement[base] + reverse_complement
    return reverse_complement


def vcf_to_vg(vcf_file, fa_file, out_file):
    SVList = ["INS", "DEL", "INV"]
    ChrSeq = genome_seq(fa_file)

    vcf_h = open(vcf_file, "r")
    out_h = open(out_file, "w")
    for line in vcf_h:
        line = line.strip()
        if line.startswith("#"):
            out_h.write("%s\n" % line)
        else:
            lines = line.split("\t")
            Chr, Pos = lines[:2]
            Pos = int(Pos)
            Ref, Alt = lines[3:5]
            Infor = lines[7]
            SVType = Infor_target_values(Infor, "SVTYPE")
            SVLength = Infor_target_values(Infor, "SVLEN")
            SVLength = int(SVLength.strip("-"))
            End = int(Infor_target_values(Infor, "END"))
            if SVType in SVList and SVLength < 1000000:
                if SVType == "DEL":
                    Alt = ChrSeq[Chr][Pos-1]
                    if Ref == "N":
                        Ref = ChrSeq[Chr][Pos-1:End]
                    out_h.write("%s\t%s\t%s\t%s\n" % ("\t".join(lines[:3]), Ref, Alt, "\t".join(lines[5:])))
                elif SVType == "INS":
                    Ref = ChrSeq[Chr][Pos-1]
                    out_h.write("%s\t%s\t%s\t%s\n" % ("\t".join(lines[:3]), Ref, Alt, "\t".join(lines[5:])))
                elif SVType == "INV":
                    Infor = Infor_substution_tag(Infor, "SVLEN", "0")
                    Ref = ChrSeq[Chr][Pos-1:End]
                    Alt = revComp(Ref)
                    out_h.write("%s\t%s\t%s\t%s\t%s\t%s\n" % ("\t".join(lines[:3]), Ref, Alt, "\t".join(lines[5:7]), Infor, "\t".join(lines[8:])))


    vcf_h.close()
    out_h.close()

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-v", "--vcf", help="The input vcf file.")
    parser.add_argument("-f", "--fasta", help="The input fasta file.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    vcf_to_vg(args.vcf, args.fasta, args.out)

if __name__ == "__main__":
    main()


        

