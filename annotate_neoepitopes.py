#!/usr/bin/env python

import argparse
import subprocess
import os
from Bio.SubsMat import MatrixInfo


def make_epitope_fasta(epitope_file, outputdir, name, fasta):
	''' Produces fasta file containing all neoepitope sequences for a sample
		
		epitope_file: parsed output from pVAC-Seq (or alternative program)
		output_dir: directory in which to write fasta
		name: sample name to distinguish file
		fasta: path to output fasta file
		
		No return value
	'''
	epitope_list = []
	
	with open(epitope_file, "r") as fh:
		for line in fh:
			line = line.split("\t")
    		epitope = line[2]
    		if epitope not in epitope_list:
    			epitope_list.append(epitope)
    
    with open(fasta, "w") as fh:
    	for epitope in epitope_list:
    		fh.write("> seq="+ epitope + "\n")
    		fh.write(epitope + "\n")

    
def run_blast(fasta, db, blastp, outputdir, name, type):
	''' Runs blastp
	
		fasta: fasta containing sequences to blast against database
		db: database to blast against
		outputdir: directory in which to write blast output
		name: sample name to distinguish file
		type: blast type - human, bacterial, or viral
		
		No return value
	'''
	outfile = outputdir + "/" + name + "." + type + ".blast.out"
	subprocess.call([blastp, "-outfmt", "6 qseqid sseqid length qstart qend sseq evalue", "-db", db, "-query", fasta, "-matrix", "BLOSUM62", "-evalue", "200000", "-ungapped", "-comp_based_stats", "F", "-out", outfile])


def score_match(pair, matrix):
	''' Gives a score from a matrix for a given pair
	
		pair: pair to score
		matrix: scoring matrix
		
		Return value: score
	''' 
    if pair not in matrix:
        return matrix[(tuple(reversed(pair)))]
    else:
        return matrix[pair]                                                                                                            


def score_pairwise(seq1, seq2, matrix):
	''' For two sequences of equal length, scores them at non-anchor positions given a matrix
	
		seq1: first sequence
		seq2: second sequence, of same length as first sequence
		matrix: scoring matrix
		
		Return value: score
	''' 
    score = 0
    last_ep = len(seq1) - 1
    for i in range(len(seq1)):
        if i != 1 and i != last_ep:
            pair = (seq1[i], seq2[i])
            score += score_match(pair, matrix)
    return score

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required=True,
            help='path input epitope data tsv file'
        )
    parser.add_argument('-o', '--outdir', type=str, required=True,
            help='path to output directory'
        )
    parser.add_argument('-d', '--dbdir', type=str, required=True,
            help='path to directory containing blastp databases'
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
    
    # Set BLOSUM62 as blosum
    blosum = MatrixInfo.blosum62
    
    # Sets paths to blast databases
    humanDB = args.dbdir + "/humanPepDB"
    bacterialDB = args.dbdir + "/bacterialPepDB"
    viralDB = args.dbdir + "/viralPepDB"
    
    # Produce fasta file containing neoepitopes
    fasta_path = outputdir + "/" + sample + ".epitopes.fasta"
    make_epitope_fasta(args.input, args.outdir, args.sample, fasta_path)
    
	
	# Run blast comparing neoepitopes to human, bacterial, and viral peptides
	run_blast(fasta_path, humanDB, args.blastp, args.outdir, args.sample, "human")
	run_blast(fasta_path, bacterialDB, args.blastp, args.outdir, args.sample, "bacterial")
	run_blast(fasta_path, viralDB, args.blastp, args.outdir, args.sample, "viral")