#!/usr/bin/env python

import argparse
import subprocess
import os
import datetime
import pickle
from Bio.SubsMat import MatrixInfo


def make_epitope_fasta(epitope_file, outputdir, name, fasta):
	''' Produces fasta file containing all neoepitope sequences for a sample
		
		epitope_file: path to parsed output file from pVAC-Seq (or alternative program)
		outputdir: path to directory in which to write fasta
		name: sample name to distinguish file (string)
		fasta: path to output fasta file
		
		No return value
	'''
	# Obtain all unique epitopes
	epitope_list = []
	with open(epitope_file, "r") as fh:
		for line in fh:
			line = line.split("\t")
			epitope = line[1]
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
	if os.path.isfile(outfile) == False:
		subprocess.call([blastp, "-outfmt", "6 qseqid sseqid length qstart qend sseq evalue", "-db", db, "-query", fasta, "-matrix", "BLOSUM62", "-evalue", "200000", "-ungapped", "-comp_based_stats", "F", "-out", outfile])
	else:
		print "Blast file for " + type + " peptides already exists - skipping"


def score_match(pair, matrix):
	''' Gives a score from a matrix for a given pair of sequences
	
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
    
def process_blast(blast_results, type, matrix, dict_dir):
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
		dict = pickle.load(open(dict_dir+"humanDict.pickle", "rb"))
	elif type == "bacterial":
		dict = pickle.load(open(dict_dir+"bacterialDict.pickle", "rb"))
	elif type == "viral":
		dict = pickle.load(open(dict_dir+"viralDict.pickle", "rb"))
	
	# Parse blast output line-by-line
	blast_dict = {}
	with open(blast_results, "r") as fh:
		for line in fh:
		    
		    # Obtain relevant data
			line = line.strip("\n").split("\t")
			epitope = line[0].split("=")[1]
			length = int(line[2])
			eval = float(line[6])
			match_pep = line[1].replace("ref", "").replace("|", "")
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
			if epitope not in blast_dict and length == len(epitope) and invalid_matches == []:
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
	
	
def add_affinities(dict, netMHCpan, allele, outputdir, name):
	''' Adds MHC binding affinities for top human blast result peptides to human blast dictionary
	
		dict: blast results dictionary
		netMHCpan: path to netMHCpan executable
		allele: HLA allele
		outputdir: path to output directory in which to write netMHCpan-associate files
		name: sample name to distinguish file (string)
		
		Return value: dictionary
	'''
	# Obtain peptide sequences
	mhc_peps = outputdir + "/" + name + ".mhc.peps"
	with open(mhc_peps, "w") as fh:
		for key in dict:
			seq = dict[key][3]
			fh.write(seq + "\n")
	
	# Run netMHCpan
	mhc_out = outputdir + "/" + name + ".mhc.out"
	subprocess.call([netMHCpan, "-a", allele, "-inptype", "1", "-p", "-xls", "-xlsfile", mhc_out, mhc_peps])
	
	# Process netMHCpan results
	affinity_dict = {}
	with open(mhc_out, "r") as fh:
		for line in fh:
			if line[0] == "0":
				line = line.strip("\n").split("\t")
				pep = line[1]
				nM = line[5]
				affinity_dict[pep] = nM
	
	# Add affinities to dictionary
	for key in dict:
		seq = dict[key][3]
		affinity = affinity_dict[seq]
		dict[key].append(affinity)
	
	return dict


def produce_annotations(epitope_file, human_dict, bacterial_dict, viral_dict, outputdir, name, allele):
	''' Produces an annotation file for predicted neoepitopes
		Contains data about relationship of neoepitope to paired normal peptide, closest human peptide from blast,
			closest bacterial peptide from blast, and closest viral peptide from blast
		
		epitope_file: path to parsed output file from pVAC-Seq (or alternative program)
		human_dict: dictionary containing human blast data
		bacterial_dict: dictionary containing bacterial blast data
		viral_dict: dictionary containing viral blast data
		outputdir: path to directory in which to write output annotation file
		name: sample name to distinguish file (string)
		allele: HLA allele used for analysis
		
		Return value: none
	'''
	outfile = outputdir + "/" + name + ".epitopes.annotated.tsv"
	with open(outfile, "w") as out:
		# Write header to outfile
		out.write("Sample\tAllele\tNeoepitope\tNeoepitope_affinity\tPaired_normal_epitope\tPaired_normal_affinity\tTranscript\tGene\tPaired_BD\tPaired_PS\tBinding_stat\tMatch_transcript\tMatch_gene\tMatch_stat\tMatch_seq\tMatch_exact\tMatch_affinity\tMatch_BD\tMatch_PS\tBac_match\tBac_seq\tBac_exact\tBac_PS\tVir_match\tVir_seq\tVir_exact\tVir_PS\n")
		# Loop through epitope file to obtain info
		with open(epitope_file, "r") as fh:
			for line in fh:
				# Extract data from epitope file
				line = line.strip("\n").split("\t")
				peptide = line[1]
				tum_bind = line[2]
				norm_pep = line[3]
				norm_bind = line[4]
				transcript = line[5]
				gene = line[6]
    			
				tum_ps = float(score_pairwise(peptide, peptide, blosum))
				if norm_pep != "NA":
					# Obtain data re: paired normal epitope
					binding_difference = float(norm_bind) - float(tum_bind)
					if float(norm_bind) > 500 and float(tum_bind) < 500 and float(norm_bind) >= 5*float(tum_bind):
						stat = "novel"
					else:
						stat = "nonnovel"
					peptide_similarity = float(score_pairwise(peptide, norm_pep, blosum))/tum_ps
				else:
					binding_difference = "NA"
					stat = "NA"
					peptide_similarity = "NA"
    			
    			# Obtain data re: closest human peptide from blast
    			if peptide in human_dict:
					blast_match_trans = human_dict[peptide][1]
					blast_match_gene = human_dict[peptide][2]
					if transcript in blast_match_trans:
						match_stat = "transcript_match"
					elif gene in blast_match_gene:
						match_stat = "gene_match"
					else:
						match_stat = "nonmatching"
						match_seq = human_dict[peptide][3]
					if match_seq == peptide:
						match_exact = "exact"
					else:
						match_exact = "inexact"
					match_ps = human_dict[peptide][4]/tum_ps
					match_affinity = human_dict[peptide][5]
					match_bd = float(match_affinity) - float(tum_bind)
				else:
					blast_match_trans = "NA"
					blast_match_gene = "NA"
					match_stat = "NA"
					match_seq = "NA"
					match_exact = "NA"
					match_ps = "NA"
					match_affinity = "NA"
					match_bd = "NA"
        		
        		# Obtain data re: closest bacterial peptide from blast
				if peptide in bacterial_dict:
					bac_match = bacterial_dict[peptide][1]
					bac_seq = bacterial_dict[peptide][2]
					if bac_seq == peptide:
						bac_exact = "exact"
					else:
						bac_exact = "inexact"
					bac_ps = bacterial_dict[peptide][3]/tum_ps
				else:
					bac_match = "NA"
					bac_seq = "NA"
					bac_exact = "NA"
					bac_ps = "NA"
				
				# Obtain data re: closest viral peptide from blast
				if peptide in viral_dict:
					vir_match = viral_dict[peptide][1]
					vir_seq = viral_dict[peptide][2]
					if vir_seq == peptide:
						vir_exact = "exact"
					else:
						vir_exact = "inexact"
					vir_ps = viral_dict[peptide][3]/tum_ps
				else:
					vir_match = "NA"
					vir_seq = "NA"
					vir_exact = "NA"
					vir_ps = "NA"
        		
				# Write data
				outline = "\t".join([allele, peptide, str(tum_bind), norm_pep, str(norm_bind), transcript, gene, str(binding_difference), str(peptide_similarity), stat, blast_match_trans, blast_match_gene, match_stat, match_seq, match_exact, str(match_affinity), str(match_bd), str(match_ps), bac_match, bac_seq, bac_exact, str(bac_ps), vir_match, vir_seq, vir_exact, str(vir_ps)])
				out.write(outline + "\n")
	

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
	parser.add_argument('-a', '--allele', type=str, required=True,
            help='HLA allele (format HLA-A02:01)'
        )
	parser.add_argument('-p', '--blastp', type=str, required=True,
            help='path to blastp executable'
        )
	parser.add_argument('-n', '--netMHCpan', type=str, required=True,
            help='path to netMHCpan executable'
        )
	args = parser.parse_args()
    
	print 'Timestamp: {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()) + " Running annotate_neoepitopes.py for " + args.input + " with allele " + args.allele
    
	# Set BLOSUM62 as blosum
	blosum = MatrixInfo.blosum62
    
	# Locate dictionary and blast database directories
	pickle_dir = os.path.dirname(__file__) + "/dictionaries/"
	blastdb_dir = os.path.dirname(__file__) + "/blast_dbs/"
    
	# Set paths to blast databases
	humanDB = blastdb_dir + "/hg38_peptide_db"
	bacterialDB = blastdb_dir + "/bacterial_peptide_db"
	viralDB = blastdb_dir + "/new_viral_pep_db"
    
	print '{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()) + " Producing fasta file with neoepitope sequences for blast..."
    
	# Produce fasta file containing neoepitopes
	fasta_path = args.outdir + "/" + args.sample + ".epitopes.fasta"
	make_epitope_fasta(args.input, args.outdir, args.sample, fasta_path)
    
	print '{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()) + " Running blast now..."
	
	# Run blast comparing neoepitopes to human, bacterial, and viral peptides
	run_blast(fasta_path, humanDB, args.blastp, args.outdir, args.sample, "human")
	run_blast(fasta_path, bacterialDB, args.blastp, args.outdir, args.sample, "bacterial")
	run_blast(fasta_path, viralDB, args.blastp, args.outdir, args.sample, "viral")
	
	print '{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()) + " Processing blast results..."
	
	# Process blast data and save to dictionaries
	hum_dict = process_blast(args.outdir+"/"+args.sample+".human.blast.out", "human", blosum, pickle_dir)
	bac_dict = process_blast(args.outdir+"/"+args.sample+".bacterial.blast.out", "bacterial", blosum, pickle_dir)
	vir_dict = process_blast(args.outdir+"/"+args.sample+".viral.blast.out", "viral", blosum, pickle_dir)
	
	print '{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()) + " Obtaining binding affinities for human blast peptides..."
	
	# Add NetMHCpan binding affinity data to blast dictionaries
	hum_dict = add_affinities(hum_dict, args.netMHCpan, args.allele, args.outdir, args.sample)
	
	print '{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()) + " Writing annotation file to " + args.outdir
	
	# Produce annotation file
	produce_annotations(args.input, hum_dict, bac_dict, vir_dict, args.outdir, args.sample, args.allele)
	
	print '{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()) + " Complete! Annotation file is " + args.outdir + "/" + args.sample + ".epitopes.annotated.tsv"
