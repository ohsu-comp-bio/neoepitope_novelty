#!/usr/bin/env python


import argparse
import subprocess
import os


def makeHumanDict(fasta):
	''' Prepares a peptide dictionary from an Ensemble human peptide fasta file
		Keys in the resulting dictionary are peptide IDs
		Values are lists: 
			position 0 is the ID for the gene of origin for the peptide
			position 1 is the ID for the transcript of origin for the peptide
					
		fasta: fasta file to make dictionary from
		
		Return value: dictionary
	'''
	dict = {}
	with open(fasta, "r") as fh:
		for line in fh:
    		if line[0] == ">":
        		line = line.strip("\n").split()
        		peptide = line[0].strip(">")
        		gene = line[3].strip("gene:").split(".")[0]
        		transcript = line[4].strip("transcript:").split(".")[0]
      			dict[peptide] = [gene, transcript]
	return dict

	
def makeBacterialDict(fasta)
	''' Prepares a peptide dictionary from a Refseq bacterial fasta file
		Keys in the resulting dictionary are peptide IDs
		Values are the bacterial species of origin for the peptides
				
		fasta: fasta file to make dictionary from
		
		Return value: dictionary
	'''
	dict = {}
	with open(fasta, "r") as fh:
		for line in fh:
    		if line[0] == ">":
    			peptide = line.split()[1]
    			species = line.split("[")[1]
    			dict[peptide] = species
	return dict
	

def makeViralDict(fasta)
	''' Prepares a peptide dictionary from a Refseq viral fasta file
		Keys in the resulting dictionary are peptide IDs
		Values are the viral species of origin for the peptides
		
		fasta: fasta file to make dictionary from
		
		Return value: dictionary
	'''
	dict = {}
	with open(fasta, "r") as fh:
		for line in fh:
    		if line[0] == ">":
    			peptide = line.split()[0].strip(">")
    			species = line.split("[")[1]
    			dict[peptide] = species
	return dict


def makeBlastpDB(makeblastdb, fasta, outputdir, title):
	''' Produces a blast protein database from a fasta file
		
		makeblastdb: path to makeblastdb executable
		fasta: peptide fasta file
		outputdir: directory to which database will be written
		title: name of database
		
		No return value.
	'''
	subprocess.call([makeblastdb, "-in", fasta, "-input_type", "fasta", "-dbtype", "prot", 
					"-title", outputdir+"/"+title, "-out", outputdir+"/"+title])
					

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--blastp', type=str, required=True,
            help='path to blastp executable'
        )
    parser.add_argument('-d', '--dbdir', type=str, required=True,
            help='path to output directory'
        )
    parser.add_argument('-h', '--humanFasta', type=str, required=True,
            help='path to human peptide fasta file from Ensembl'
        )
    parser.add_argument('-b', '--bacterialFasta', type=str, required=True,
            help='path to combined bacterial peptide fasta file derived from RefSeq release'
        )
    parser.add_argument('-v', '--viralFasta', type=str, required=True,
            help='path to combined viral peptide fasta file derived from RefSeq release'
        )
    args = parser.parse_args()
    
    # Produce human peptide dictionary and write to lib
    human_peptides = makeHumanDict(args.humanFasta)
    with open("lib/humanPeptideDict.py", "w") as fh:
    	fh.write(human_peptides)
    
    # Produce bacterial peptide dictionary and write to lib
    bacterial_peptides = makeBacterialDict(args.bacterialFasta)
    with open("lib/bacterialPeptideDict.py", "w") as fh:
    	fh.write(bacterial_peptides)    
    
    # Produce viral peptide dictionary and write to lib
    viral_peptides = makeViralDict(args.viralFasta)
    with open("lib/viralPeptideDict.py", "w") as fh:
    	fh.write(viral_peptides)
    	
    makeBlastpDB(args.blastp, args.humanFasta, args.dbdir, "humanPepDB")
    makeBlastpDB(args.blastp, args.bacterialFasta, args.dbdir, "bacterialPepDB")
    makeBlastpDB(args.blastp, args.viralFasta, args.dbdir, "viralPepDB")