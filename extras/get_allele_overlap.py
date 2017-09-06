#!/usr/bin/env python

import random

from optparse import OptionParser
p = OptionParser(usage = "python get_allele_overlap.py -i <infile> -o <outfile> -a <allele_file>")
p.add_option("-i", action="store", dest="infile", help="Path to input file")
p.add_option("-o", action="store", dest="outfile", help="Path to output file")
p.add_option("-a", action="store", dest="alleles", help="Path to allele file")
opts, args = p.parse_args()
infile = opts.infile
outfile = opts.outfile
allele_file = opts.alleles

allele_fh = open(allele_file, "r")
alleles = allele_fh.readline().strip("\n").split(",")
allele_fh.close()

a_alleles = []
for allele in alleles:
	if "HLA-A" in allele:
		a_alleles.append(allele)
b_alleles = []
for allele in alleles:
	if "HLA-B" in allele:
		b_alleles.append(allele)
c_alleles = []
for allele in alleles:
	if "HLA-C" in allele:
		c_alleles.append(allele)

allele_sets = []
for i in range(0,1000):
	aset = random.sample(a_alleles, 2)
	bset = random.sample(b_alleles, 2)
	cset = random.sample(c_alleles, 2)
	set = aset + bset + cset
	allele_sets.append(set)

in_fh = open(infile, "r")
out_fh = open(outfile, "w")

epitope_dict = {}
for line in in_fh:
	if "Disease" in line:
		continue
	line = line.strip("\n").split("\t")
	allele = line[1].strip('"')
	adj_allele = "HLA-" + allele[0] + "*" + allele[1:3] + ":" + allele[3:]
	epitope = line[2].strip('"')
	tumor_aff = float(line[3].strip('"'))
	normal_aff = float(line[5].strip('"'))
	if tumor_aff < 500 and normal_aff > 500 and normal_aff >= (5*tumor_aff) and epitope not in epitope_dict:
		epitope_dict[epitope] = [adj_allele]
	elif tumor_aff < 500 and normal_aff > 500 and normal_aff >= (5*tumor_aff) and epitope in epitope_dict:
		if adj_allele not in epitope_dict[epitope]:
			epitope_dict[epitope].append(adj_allele)

for group in allele_sets:
	share_data = {}
	share_counts = {1 : 0, 2 : 0, 3 : 0, 4 : 0, 5 : 0, 6 : 0}
	for epitope in epitope_dict:
		alleles = [HLA for HLA in epitope_dict[epitope] if (HLA in group)]
		count = len(alleles)
		if count != 0:
			share_counts[count] += 1
	
	this_set = ",".join(group)

	outline = this_set + "\t" + str(share_counts[1]) + "\t" + str(share_counts[2]) + "\t" + str(share_counts[3]) + "\t" + str(share_counts[4]) + "\t" + str(share_counts[5]) + "\t" + str(share_counts[6]) + "\n"
	out_fh.write(outline)
	
	
in_fh.close()
out_fh.close()
