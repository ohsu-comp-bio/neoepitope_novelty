#!/usr/bin/env python

import argparse
import subprocess
import os
from Bio.SubsMat import MatrixInfo


def make_epitope_fasta(epitope_file, output_dir, sample):
	epitope_list = []
	
	with open(epitope_file, "r") as fh:
		for line in fh:
			line = line.split("\t")
    		epitope = line[2]
    		if epitope not in epitope_list:
    			epitope_list.append(epitope)
    
    fasta = output_dir + "/" + sample + ".epitopes.fasta"
    with open(fasta, "w") as fh:
    	for epitope in epitope_list:
    		fh.write("> seq="+ epitope + "\n")
    		fh.write(epitope + "\n")
    
    return fasta


def score_match(pair, matrix):
    if pair not in matrix:
        return matrix[(tuple(reversed(pair)))]
    else:
        return matrix[pair]                                                                                                            


def score_pairwise(seq1, seq2, matrix):
    score = 0
    last_ep = len(seq1) - 1
    for i in range(len(seq1)):
        if i != 1 and i != last_ep:
            pair = (seq1[i], seq2[i])
            score += score_match(pair, matrix)
    return score
blosum = MatrixInfo.blosum62


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required=True,
            help='path input epitope data tsv file'
        )
    parser.add_argument('-o', '--outdir', type=str, required=True,
            help='path to output directory'
        )
    parser.add_argument('-s', '--sample', type=str, required=True,
            help='sample name'
        )
    parser.add_argument('-p', '--blastp', type=str, required=True,
            help='path to blastp executable'
        )
    parser.add_argument('-n', '--netMHCpan', type=str, required=True,
            help='path to netMHCpan executable'
        )
    args = parser.parse_args()
    
    make_epitope_fasta(args.input, args.outdir, args.sample)
