import argparse
import os
import collections

#usage: python ~/github/PanSV/script/buscoSummary.py --file /home/wuzhikun/Project/PanVigna/Evaluation/BUSCO/ragtag/CPG13/run_embryophyta_odb10/short_summary.txt --out temp.txt

def busco_summary(in_file, out_file):
    """
    in_file:
    	***** Results: *****

	C:98.9%[S:96.7%,D:2.2%],F:0.4%,M:0.7%,n:1614	   
	1597	Complete BUSCOs (C)			   
	1561	Complete and single-copy BUSCOs (S)	   
	36	Complete and duplicated BUSCOs (D)	   
	7	Fragmented BUSCOs (F)			   
	10	Missing BUSCOs (M)			   
	1614	Total BUSCO groups searched	

    out_file:
    Sample	Complete (%)	Fragment (%)	Missing (%)	Total
    CPG13	98.9	0.4	0.7	1614
    """
    Samples = []
    BuscoSummary = collections.defaultdict(dict)
    files = in_file.split(',')
    for f in files:
        name = f.strip().split('/')[-3]
        Samples.append(name)
        in_h = open(f, "r")
        for line in in_h:
            line = line.strip()
            if line.startswith("C"):
                lines = line.split(",")
                c, d, f, m, n = lines
                BuscoSummary[name]["Complete"] = c.split(":")[1].split("[")[0].strip("%")
                BuscoSummary[name]["Fragment"] = f.split(":")[1].strip("%")
                BuscoSummary[name]["Missing"] = m.split(":")[1].strip("%")
                BuscoSummary[name]["Total"] = n.split(":")[1]
        in_h.close()

    out_h = open(out_file, "w")
    out_h.write("Sample\tComplete (%)\tFragment (%)\tMissing (%)\tTotal\n")
    for s in sorted(Samples):
        Complete = BuscoSummary[s]["Complete"]
        Fragment = BuscoSummary[s]["Fragment"]
        Missing = BuscoSummary[s]["Missing"]
        Total = BuscoSummary[s]["Total"]
        out_h.write("%s\t%s\n" % (s,"\t".join([Complete, Fragment, Missing, Total])))
    out_h.close()



def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-f", "--file", help="The file of the genome stats.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    busco_summary(args.file, args.out)

if __name__ == "__main__":
    main()

