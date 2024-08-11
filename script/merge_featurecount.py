#!/usr/bin/env python
from optparse import OptionParser
import collections
import time
import sys
import os 
#ueage:python /home/wzk/test_data/merge_featurecount.py -i '/home/wzk/test_data/supp/expression/MCF-7_Reads_count.txt,/home/wzk/test_data/supp/expression/MCF-7-2835_Reads_count.txt,/home/wzk/test_data/supp/expression/MCF-7-3440_Reads_count.txt,/home/wzk/test_data/supp/expression/MCF-7-3492_Reads_count.txt' -o /home/wzk/test_data/supp/DEG/All_reads_counts.xls


def this_time():
	timeformat = "%Y-%m-%d %H:%M:%S"
	return time.strftime(timeformat, time.localtime())

def read_featurecount(file):
	'''Read each featurecount file and then build the hash of gene and the record, then return the sample name and this hash'''
	gene_feature = collections.Counter()
	in_handle = open(file,'r')
	for line in in_handle:
		if line.strip().startswith("#"):
			continue
		elif line.strip().startswith("Geneid"):
			#sample = os.path.dirname(line.strip().split("\t")[-1]).split("/")[-1] #for the result of STAR sligned
			sample = line.strip().split("/")[-1].strip("Aligned.sortedByCoord.out.bam")
			print(sample)
		elif line.strip():
			lines = line.strip().split("\t")
			gene = lines[0]
			gene_feature[gene] = line.strip()
	in_handle.close()
	try:
		return sample, gene_feature
	except:
		print("Check whether the sample exists or the input file is empty!")
		sys.exit()
		
			
def merge_files(files,out):
	'''Merge all featuecount result files'''
	sample_featurecount = collections.defaultdict(dict)
	all_files = files.strip().split(",")
	#print(all_files)
	for f in all_files:
		sample, feature_hash = read_featurecount(f)
		sample_featurecount[sample] = feature_hash
	samples = sample_featurecount.keys()
	#the gene number is equal in all samples, so merge the data by regarding the gene name as keys
	sorted_samples = sorted(samples)
	#print(sorted_samples)
	out_handle = open(out,'w')
	out_handle.write("Geneid\tChr\tStart\tEnd\tStrand\tLength\t%s\n" % '\t'.join(sorted_samples))
	sorted_gene = sorted(sample_featurecount[sorted_samples[0]].keys())
	for k in sorted_gene:
		#the information except the read count
		info = sample_featurecount[sorted_samples[0]][k].split("\t")[:-1]
		values = [sample_featurecount[sorted_samples[i]][k].split("\t")[-1] for i in range(len(sorted_samples))]
		#print(values)
		out_handle.write("%s\t%s\n" % ('\t'.join(info),'\t'.join(values)))
	out_handle.close()
	
def main():
	parser = OptionParser(description="Merge the count reads derived from featurecount for all smples.")
	parser.add_option('-i', '--input',dest='input',help='The featureCounts reurlt files which are seperate with comma.')
	parser.add_option('-o','--output',dest='output',help='The output file.')
	(options, args) = parser.parse_args()
	merge_files(options.input, options.output)
	
if __name__ == "__main__":
	main()

	
	
