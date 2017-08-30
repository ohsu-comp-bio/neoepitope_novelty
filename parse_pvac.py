#!/usr/bin/env python


import argparse
import datetime


def parse_pvac(pvac, output):
	''' Parses output of pVAC-seq to obtain relevant data for annotating neoepitopes
		
		pvac: path to results file from pVAC-seq
		output: path to output file in which to write parsed results
		name: sample name to distinguish file (string)
		
		Return value: none
	'''
	with open(output, "w") as out:
		with open(pvac, "r") as fh:
			for line in fh:
				line = line.strip("\n").split("\t")
				if line[0] != "Chromosome":
					transcript = line[5]
					gene = line[6]
					allele = line[11]
					neoepitope = line[15]
					normal_epitope = line[16]
					normal_bind = line[32]
					tumor_bind = line[33]
					outline = "\t".join([allele, neoepitope, tumor_bind, normal_epitope, normal_bind, transcript, gene])
					out.write(outline + "\n")
					

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', type=str, required=True,
            help='path to input file'
        )
    parser.add_argument('-o', '--outfile', type=str, required=True,
            help='path to output file'
        )
    args = parser.parse_args()
    
    
    print '{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()) + " Beginning parse for " + args.infile
    
    # Parse pVAC-Seq results
    parse_pvac(args.infile, args.outfile)
    
    print '{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()) + " Complete! Output is " + args.outfile