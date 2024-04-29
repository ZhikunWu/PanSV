import argparse

#usage: python CheckInvalidBase.py --fasta /home/wuzhikun/Project/PanSV/Scaffold/CPG15.genome.fasta > temp.txt

def check_invalid_bases(filename):
    valid_bases = set(['A', 'T', 'C', 'G', 'N'])
    invalid_bases = set()
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            # 如果该行是注释行，跳过
            if line.startswith('>'):
                continue
            
            # 检查非法的碱基
            for base in line:
                if base not in valid_bases:
                    invalid_bases.add(base)
                    

    if invalid_bases:
        print(f"Invalid bases found: {' '.join(invalid_bases)}")
    else:
        print("No invalid bases found.")



def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-f", "--fasta", help="The fasta file.")
    args = parser.parse_args()
    check_invalid_bases(args.fasta)

if __name__ == "__main__":
    main()



