import argparse
import collections
import operator


#usage: python ~/github/Centromere/script/CenRegionRepeat.py --region  Zm-Mo17_srf_aln_cluster_target.txt --record Zm-Mo17_srf_aln.bed --out temp --satellite sat_temp

#usage: python /home/wuzhikun/github/Centromere/script/CenRegionRepeat.py --record /home/wuzhikun/Project/Centromere/Centromere/SRF/Chromosome/Paspalum_notatum_srf_aln.bed --satellite /home/wuzhikun/Project/Centromere/Centromere/SRF/Chromosome/Paspalum_notatum_srf_aln.satellite.txt


def target_region(bed_file):
    """
    bed_file:
    CM039150.1	136947826	138864444	0.495
    CM039151.1	98033838	98775948	0.329
    CM039152.1	87918991	88620414	0.494
    """
    ChrRegion = {}
    in_h = open(bed_file, "r")
    for line in in_h:
        lines = line.strip().split("\t")
        Chr, Start, End = lines[:3]
        Start = int(Start)
        End = int(End)
        ChrRegion[Chr] = [Start, End]
    in_h.close()
    return ChrRegion


def target_region_record(bed_file, record_file, out_file, sat_file):
    ChrRegion = target_region(bed_file)

    ChrSatelliteNum = collections.defaultdict(lambda: collections.Counter())
    AllSatelliteNum = collections.Counter()
    in_h = open(record_file, "r")
    out_h = open(out_file, "w")
    for line in in_h:
        line = line.strip()
        lines = line.split("\t")
        Chr, Start, End = lines[:3]
        if Chr in ChrRegion:
            Region = ChrRegion[Chr]
            Start = int(Start)
            End = int(End)
            if Start >= Region[0] and End <= Region[1]:
                out_h.write("%s\n" % line)

                ### calculate satellite
                satellite = lines[3]
                sat = int(satellite.split("-")[-1])
                ChrSatelliteNum[Chr][sat] += 1
                AllSatelliteNum[sat] += 1
    in_h.close()
    out_h.close()

    sat_h = open(sat_file,  "w")
    ### sum all sataellite
    SortSatelliteAll =  sorted(AllSatelliteNum.items(), key=operator.itemgetter(1), reverse=True)
    SortSatelliteAll1 = ["%d|%d" % i for i in SortSatelliteAll]
    sat_h.write("%s\t%s\n" % ("Sum" ,"\t".join(SortSatelliteAll1)))

    ### sum chromosomes
    Chrs = sorted(list(ChrSatelliteNum.keys()))
    for c in Chrs:
        SatelliteNum = ChrSatelliteNum[c]
        SortSatellite =  sorted(SatelliteNum.items(), key=operator.itemgetter(1), reverse=True)
        SortSatellite1 = ["%d|%d" % i for i in SortSatellite]
        sat_h.write("%s\t%s\n" % (c, "\t".join(SortSatellite1)))
    sat_h.close()


def all_region_record(record_file, sat_file):

    ChrSatelliteNum = collections.defaultdict(lambda: collections.Counter())
    AllSatelliteNum = collections.Counter()
    in_h = open(record_file, "r")
    for line in in_h:
        line = line.strip()
        lines = line.split("\t")
        Chr, Start, End = lines[:3]
        ### calculate satellite
        satellite = lines[3]
        sat = int(satellite.split("-")[-1])
        ChrSatelliteNum[Chr][sat] += 1
        AllSatelliteNum[sat] += 1
    in_h.close()    

    sat_h = open(sat_file,  "w")
    ### sum all sataellite
    SortSatelliteAll =  sorted(AllSatelliteNum.items(), key=operator.itemgetter(1), reverse=True)
    SortSatelliteAll1 = ["%d|%d" % i for i in SortSatelliteAll]
    sat_h.write("%s\t%s\n" % ("Sum" ,"\t".join(SortSatelliteAll1)))

    ### sum chromosomes
    Chrs = sorted(list(ChrSatelliteNum.keys()))
    for c in Chrs:
        SatelliteNum = ChrSatelliteNum[c]
        SortSatellite =  sorted(SatelliteNum.items(), key=operator.itemgetter(1), reverse=True)
        SortSatellite1 = ["%d|%d" % i for i in SortSatellite]
        sat_h.write("%s\t%s\n" % (c, "\t".join(SortSatellite1)))
    sat_h.close()

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-b', '--region', help='The input bed file containing centromere region.')
    parser.add_argument('-o', '--out', help='The output file.')
    parser.add_argument("-r", "--record", help="The record file.")
    parser.add_argument("-s", "--satellite", help="The satellite number statistics.")
    args = parser.parse_args()
    if args.region:
        target_region_record(args.region, args.record, args.out, args.satellite)
    else:
        all_region_record(args.record, args.satellite)        


if __name__ == '__main__':
    main()


