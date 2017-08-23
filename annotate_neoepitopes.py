#!/usr/bin/env python

import argparse
import subprocess
import os
import lib
from Bio.SubsMat import MatrixInfo


def make_epitope_fasta(epitope_file, outputdir, name, fasta):
	''' Produces fasta file containing all neoepitope sequences for a sample
		
		epitope_file: path to parsed output file from pVAC-Seq (or alternative program)
		output_dir: path to directory in which to write fasta
		name: sample name to distinguish file (string)
		fasta: path to output fasta file
		
		No return value
	'''
	# Obtain all unique epitopes
	epitope_list = []
	with open(epitope_file, "r") as fh:
		for line in fh:
			line = line.split("\t")
    		epitope = line[2]
    		if epitope not in epitope_list:
    			epitope_list.append(epitope)
    
    # Write unique epitopes to fasta file
    with open(fasta, "w") as fh:
    	for epitope in epitope_list:
    		fh.write("> seq="+ epitope + "\n")
    		fh.write(epitope + "\n")

    
def run_blast(fasta, db, blastp, outputdir, name, type):
	''' Runs blastp
	
		fasta: path fasta containing sequences to blast against database
		db: path to database to blast against
		outputdir: path to directory in which to write blast output
		name: sample name to distinguish file (string)
		type: blast type - human, bacterial, or viral (string)
		
		No return value
	'''
	outfile = outputdir + "/" + name + "." + type + ".blast.out"
	subprocess.call([blastp, "-outfmt", "6 qseqid sseqid length qstart qend sseq evalue", "-db", db, "-query", fasta, "-matrix", "BLOSUM62", "-evalue", "200000", "-ungapped", "-comp_based_stats", "F", "-out", outfile])


def score_match(pair, matrix):
	''' Gives a score from a matrix for a given pair
	
		pair: pair to score (tuple)
		matrix: scoring matrix (matrix)
		
		Return value: score
	''' 
    if pair not in matrix:
        return matrix[(tuple(reversed(pair)))]
    else:
        return matrix[pair]                                                                                                            


def score_pairwise(seq1, seq2, matrix):
	''' For two sequences of equal length, scores them at non-anchor positions given a matrix
	
		seq1: first sequence (string)
		seq2: second sequence, of same length as first sequence (string)
		matrix: scoring matrix (matrix)
		
		Return value: score
	''' 
    score = 0
    last_ep = len(seq1) - 1
    for i in range(len(seq1)):
        if i != 1 and i != last_ep:
            pair = (seq1[i], seq2[i])
            score += score_match(pair, matrix)
    return score
    
def process_blast(blast_results, type, matrix):
	''' Processes results of blastp to obtain dictionary of best results
		Keys are neoepitope sequences
		Values are lists of associated blast data: E value, sequence of match, raw protein similarity score
		For human blast, transcript and gene of match peptide are also stored
		For bacterial or viral blast, species of match peptide is also stored
	
		blast_results: path to file containing results of blastp
		type: blast type - human, bacterial, or viral (string)
		matrix: scoring matrix to use for comparing sequences (matrix)
	
		Return value: dictionary
	'''
	# Set peptide dictionary
	if type == "human":
		dict = lib.humanPepDB
	elif type == "bacterial":
		dict = lib.bacterialPepDB
	elif type == "viral":
		dict = lib.viralPepDB
	
	# Parse blast output line-by-line
	blast_dict = {}
	with open(blast_results, "r") as fh:
		for line in fh:
		    
		    # Obtain relevant data
		    line = line.strip("\n").split("\t")
    		epitope = line[0].split("=")[1]
    		length = int(line[2])
    		eval = float(line[6])
    		match_pep = line[1]
    		if type == "human":
    			match_transcript = dict[match_pep][1]
    			match_gene = dict[match_pep][0]
    		else:
    			match_species = dict[match_pep]
    		match_seq = line[5]
    		
    		# Check for presence of invalid characters in match seq
    		invalids = ["B", "J", "O", "U", "X", "Z", "*"]
    		invalid_matches = []
    		for char in invalids:
        		if char in match_seq:
            		invalid_matches.append(char)
			
			# If epitope is not already in dictionary, add it
			if epitope not in blast_dict and length == len(epitope) and invalid == []:
				match_ps = score_pairwise(epitope, match_seq, matrix)
				if type == "human":
					blast_dict[epitope] = [eval, match_transcript, match_gene, match_seq, match_ps]
				else:
					blast_dict[epitope] = [eval, match_species, match_seq, match_ps]
			
			# If epitope is in dictionary, but E value for this entry is better, replace data
			elif epitope in blast_dict and length == len(epitope) and eval < blast_dict[epitope][0] and invalid_matches == []:		
				match_ps = score_pairwise(epitope, match_seq, matrix)
				if type == "human":
					blast_dict[epitope] = [eval, match_transcript, match_gene, match_seq, match_ps]
				else:
					blast_dict[epitope] = [eval, match_species, match_seq, match_ps]
			
			# If epitope is in dictionary and E value for this entry is equivalent, compare further
			# If the match sequence is the same as previous entry, store this entry's data too
			# If match sequence is a better match to neoepitope, replace data		
			elif epitope in blast_dict and length == len(epitope) and eval == blast_dict[epitope][0] and invalid_matches == []:
				if type == "human":
					if match_seq == blast_dict[epitope][3] and match_transcript not in blast_dict[epitope][1] and match_gene not in blast_dict[epitope][2]:
            			blast_dict[epitope][1] = blast_dict[epitope][1] + "," + match_transcript
            			blast_dict[epitope][2] = blast_dict[epitope][2] + "," + match_gene
            		else:
            			match_ps = score_pairwise(epitope, match_seq, matrix)
            			if match_seq == epitope or match_ps > blast_dict[epitope][4]:
                			blast_dict[epitope] = [eval, match_transcript, match_gene, match_seq, match_ps]
				else:
					if match_seq == blast_dict[epitope][2] and match_species not in blast_dict[epitope][1]:
            			blast_dict[epitope][1] = blast_dict[epitope][1] + "," + match_species
        			else:
           				match_ps = score_pairwise(epitope, match_seq, blosum)
            			if match_seq == epitope or match_ps > blast_dict[epitope][3]:
                		blast_dict[epitope] = [eval, match_species, match_seq, match_ps]
	
	return blast_dict
	

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