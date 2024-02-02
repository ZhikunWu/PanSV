import argparse
import json
import pandas as pd
# from prettytable import PrettyTable


#usage: python UltraJson2Table.py --json /home/wuzhikun/Project/Centromere/genome/Models/Vitis_vinifera.ultra.json --out temp.txt

def ultra_json_table(json_file, out_file):
    """
    json_file:
    "Repeats": [{"PassID": 0,
    "SequenceName": "Chr1",
    "Start": 0,
    "Length": 3718,
    "Period": 7,
    "Score": 5276.733709,
    "Log2 Pval": -1979.551511,
    "Substitutions": 41,
    "Insertions": 16,
    "Deletions": 34,
    "Consensus": "CCAAAAC",
    "Sequence": "TAAAACCTAAAACCTAAAACCTAAAACCTAAAACCTAAACCCTAAACCCTAAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTA
    """
    f = open(json_file)
    data = json.load(f)
    passes = data["Passes"]
    repeats = data["Repeats"]
    ### output the repeats
    out_h = open(out_file, "w")
    out_h.write("Chr\tStart\tLength\tPeriod\tConsensus\tSequence\n")
    for repeat in repeats:
        out_h.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (repeat["SequenceName"], repeat["Start"], repeat["Length"], repeat["Period"], repeat["Consensus"], repeat["Sequence"]))
    out_h.close()



def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-j', '--json', help='The input json file.')
    parser.add_argument('-o', '--out', help='The output file.')
    args = parser.parse_args()
    ultra_json_table( args.json, args.out)


if __name__ == '__main__':
    main()

