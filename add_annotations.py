#!/usr/bin/env python

# Imports required modules
import subprocess
import os
from Bio.SubsMat import MatrixInfo

# Imports and initializes optparse module                                                                                                                           
from optparse import OptionParser
p = OptionParser(usage = "python compile_results.py -i <indir> -s <sample>")
p.add_option("-i", action="store", dest="infile", help="Path to input file directory")
p.add_option("-o", action="store", dest="outdir", help="Path to output directory")
p.add_option("-s", action="store", dest="sample", help="Sample name")
opts, args = p.parse_args()
sample = opts.sample
infile = opts.infile
outdir = opts.outdir
if outdir[-1] != "/":
    outdir = outdir + "/"

mhcdir = "/mnt/scratch/mhc/"
subprocess.call(["mkdir", mhcdir])

# Defines functions and matrix for protein similarity score                                                                                                                                                                                                                         
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

# Sets paths to relevant files
blastdir = outdir + "blast/"
hum_db_path = "/home/exacloud/lustre1/CompBio/projs/RThompson/immunotherapy/hg38/hg38_peptide_db"
human_peptide_fasta = "/mnt/lustre1/CompBio/projs/RThompson/immunotherapy/hg38/Homo_sapiens.GRCh38.pep.all.fa"
bac_db_path = "/home/exacloud/lustre1/CompBio/projs/RThompson/immunotherapy/jobs/bacterial_peptides/bacterial_pep_fastas/bacterial_peptide_db"
bacterial_peptide_fasta = "/mnt/lustre1/CompBio/projs/RThompson/immunotherapy/jobs/bacterial_peptides/bacterial_pep_fastas/combined_bacteria.faa"
vir_db_path = "/home/exacloud/lustre1/CompBio/projs/RThompson/immunotherapy/jobs/bacterial_peptides/viral_pep_fastas/new_viral_pep_db"
viral_peptide_fasta = "/mnt/lustre1/CompBio/projs/RThompson/immunotherapy/jobs/bacterial_peptides/viral_pep_fastas/combined.viral.protein.faa"
outfile = outdir + sample + ".complete.tsv"
out_fh = open(outfile, "w")
out_fh.write("Sample\tAllele\tTumor_peptide\tTumor_affinity\tNormal_peptide\tNormal_affinity\tBinding_difference\tBLOSUM\tStat\tMatch_transcript\tMatch_gene\tBlast_score\tMatch_seq\tMatch_affinity\tMatch_BD\tMatch_PS\tBac_score\tBac_match\tBac_seq\tBac_PS\tVir_score\tVir_match\tVir_seq\tVir_PS\n")

print "Making blast directory...\n"

# Creates working directory
subprocess.call(["mkdir", blastdir])
blast_out_hum = blastdir + sample + ".hum.blast.out"
blast_out_bac = blastdir + sample + ".bac.blast.out"
blast_out_vir = blastdir + sample + ".vir.blast.out"
fasta = blastdir + sample + ".fasta"

print "\nDone! Collecting unique epitopes now...\n"

epitope_list = []

fasta_fh = open(fasta, "w")
fh1 = open(infile, "r")
for line in fh1:
    line = line.split("\t")
    epitope = line[2]
    normal_pep = line[4]
    if epitope not in epitope_list:
    	epitope_list.append(epitope)
    	fasta_fh.write("> seq="+ epitope + "\n")
    	fasta_fh.write(epitope + "\n")
fasta_fh.close()

epitope_list = []

print "\nDone! Running blast for unique epitopes now...\n"

if os.path.isfile(blast_out_hum) == False:
    print "\nRunning blast for human peptides first...\n"
    subprocess.call(["blastp", "-num_threads", "12", "-outfmt", "6 qseqid sseqid length qstart qend sseq evalue", "-db", hum_db_path, "-query", fasta, "-matrix", "BLOSUM62", "-evalue", "200000", "-ungapped", "-comp_based_stats", "F", "-out", blast_out_hum])
if os.path.isfile(blast_out_bac) == False:
    print "\nNow running for bacterial peptides...\n" 
    subprocess.call(["blastp", "-num_threads", "12", "-outfmt", "6 qseqid sseqid length qstart qend sseq evalue", "-db", bac_db_path, "-query", fasta, "-matrix", "BLOSUM62", "-evalue", "200000", "-ungapped", "-comp_based_stats", "F", "-out", blast_out_bac])
if os.path.isfile(blast_out_vir) == False:
    print "\nNow running for viral peptides...\n"
    subprocess.call(["blastp", "-num_threads", "12", "-outfmt", "6 qseqid sseqid length qstart qend sseq evalue", "-db", vir_db_path, "-query", fasta, "-matrix", "BLOSUM62", "-evalue", "200000", "-ungapped", "-comp_based_stats", "F", "-out", blast_out_vir])

print "\nDone! Collecting peptide to gene/species databases...\n"

human_peptide_dict = {}
hum_fasta_fh = open(human_peptide_fasta, "r")
for line in hum_fasta_fh:
    if line[0] == ">":
        line = line.strip("\n").split()
        peptide = line[0].strip(">")
        gene = line[3].strip("gene:").split(".")[0]
        transcript = line[4].strip("transcript:").split(".")[0]
        human_peptide_dict[peptide] = [gene, transcript]
hum_fasta_fh.close()

bacterial_peptide_dict = {}
bac_fasta_fh = open(bacterial_peptide_fasta, "r")
for line in bac_fasta_fh:
    if line[0] == ">":
        line = line.strip("\n").split()
        peptide = line[0].split("|")[1]
        for item in line:
            if item[0] == "[":
                genus = item.strip("[").strip("]")
                bacterial_peptide_dict[peptide] = genus
bac_fasta_fh.close()

viral_peptide_dict = {}
vir_fasta_fh = open(viral_peptide_fasta, "r")
for line in vir_fasta_fh:
    if line[0] == ">":
        linev1 = line.strip("\n").split()
        peptide = linev1[0].strip(">")
        line = line.split("[")
        species = line[1].strip("]").split()
        species = ".".join(species)
        viral_peptide_dict[peptide] = species
vir_fasta_fh.close()

print "\nDone! Processing blast data...\n"

hum_dict = {}
hum_fh = open(blast_out_hum, "r")
for line in hum_fh:
    line = line.strip("\n").split("\t")
    epitope = line[0].split("=")[1]
    length = int(line[2])
    eval = float(line[6])
    match_pep = line[1]
    match_transcript = human_peptide_dict[match_pep][1]
    match_gene = human_peptide_dict[match_pep][0]
    match_seq = line[5]
    invalids = ["B", "J", "O", "U", "X", "Z", "*"]
    invalid_matches = []
    for char in invalids:
        if char in match_seq:
            invalid_matches.append(char)
    if epitope not in hum_dict and length == len(epitope) and invalid_matches == []:
        match_ps = score_pairwise(epitope, match_seq, blosum)
        hum_dict[epitope] = [eval, match_transcript, match_gene, match_seq, match_ps]
    elif epitope in hum_dict and length == len(epitope) and eval < hum_dict[epitope][0] and invalid_matches == []:
        match_ps = score_pairwise(epitope, match_seq, blosum)
        hum_dict[epitope] = [eval, match_transcript, match_gene, match_seq, match_ps]
    elif epitope in hum_dict and length == len(epitope) and eval == hum_dict[epitope][0] and invalid_matches == []:
        if match_seq == hum_dict[epitope][3]:
            hum_dict[epitope][1] = hum_dict[epitope][1] + "," + match_transcript
            hum_dict[epitope][2] = hum_dict[epitope][2] + "," + match_gene
        else:
            match_ps = score_pairwise(epitope, match_seq, blosum)
            if match_seq == epitope or match_ps > hum_dict[epitope][4]:
                hum_dict[epitope] = [eval, match_transcript, match_gene, match_seq, match_ps]
        
hum_fh.close()        


bac_dict = {}
bac_fh = open(blast_out_bac, "r")
for line in bac_fh:
    line = line.strip("\n").split("\t")
    epitope = line[0].split("=")[1]
    length = int(line[2])
    eval = float(line[6])
    bac_pep = line[1].split("|")[1]
    bac_match = bacterial_peptide_dict[bac_pep]
    match_seq = line[5]
    invalids = ["B", "J", "O", "U", "X", "Z", "*"]
    invalid_matches = []
    for char in invalids:
        if char in match_seq:
            invalid_matches.append(char)
    if epitope not in bac_dict and length == len(epitope) and invalid_matches == []:
        match_ps = score_pairwise(epitope, match_seq, blosum)
        bac_dict[epitope] = [eval, bac_match, match_seq, match_ps]
    elif epitope in bac_dict and length == len(epitope) and eval < bac_dict[epitope][0] and invalid_matches == []:
        match_ps = score_pairwise(epitope, match_seq, blosum)
        bac_dict[epitope] = [eval, bac_match, match_seq, match_ps]
    elif epitope in bac_dict and length == len(epitope) and eval == bac_dict[epitope][0] and invalid_matches == []:
        if match_seq == bac_dict[epitope][2]:
            bac_dict[epitope][1] = bac_dict[epitope][1] + "," + bac_match
        else:
            match_ps = score_pairwise(epitope, match_seq, blosum)
            if match_seq == epitope or match_ps > bac_dict[epitope][3]:
                bac_dict[epitope] = [eval, bac_match, match_seq, match_ps]
bac_fh.close()

vir_dict = {}
vir_fh = open(blast_out_vir, "r")
for line in vir_fh:
    line = line.strip("\n").split("\t")
    epitope = line[0].split("=")[1]
    length = int(line[2])
    eval = float(line[6])
    vir_pep = line[1]
    vir_match = viral_peptide_dict[vir_pep]
    match_seq = line[5]
    invalids = ["B", "J", "O", "U", "X", "Z", "*"]
    invalid_matches = []
    for char in invalids:
        if char in match_seq:
            invalid_matches.append(char)
    if epitope not in vir_dict and length == len(epitope) and invalid_matches == []:
        match_ps = score_pairwise(epitope, match_seq, blosum)
        vir_dict[epitope] = [eval, vir_match, match_seq, match_ps]
    elif epitope in vir_dict and length == len(epitope) and eval < vir_dict[epitope][0] and invalid_matches == []:
        match_ps = score_pairwise(epitope, match_seq, blosum)
        vir_dict[epitope] = [eval, vir_match, match_seq, match_ps]
    elif epitope in vir_dict and length == len(epitope) and eval == vir_dict[epitope][0] and invalid_matches == []:
        if match_seq == vir_dict[epitope][2]:
            vir_dict[epitope][1] = vir_dict[epitope][1] + "," + vir_match
        else:
            match_ps = score_pairwise(epitope, match_seq, blosum)
            if match_seq == epitope or match_ps > vir_dict[epitope][3]:
                vir_dict[epitope] = [eval, vir_match, match_seq, match_ps]
vir_fh.close()

affinity_dict = {}

# For each line in input file, finds best blast match from human, bacterial and viral peptides
fh1.seek(0)
for line in fh1:

    # Parse infile
    line = line.strip("\n").split("\t")
    peptide = line[2]
    tum_bind = line[3]
    norm_pep = line[4]
    norm_bind = line[5]
    binding_difference = float(norm_bind) - float(tum_bind)
    if float(norm_bind) > 500 and float(tum_bind) < 500 and float(norm_bind) >= 5*float(tum_bind):
        stat = "immunogenic"
    else:
        stat = "nonimmunogenic"
    tum_ps = float(score_pairwise(peptide, peptide, blosum))
    peptide_similarity = score_pairwise(peptide, norm_pep, blosum)/tum_ps
	
    if peptide in hum_dict:
        blast_match_trans = hum_dict[peptide][1]
        blast_match_gene = hum_dict[peptide][2]
        blast_score = hum_dict[peptide][0]
        match_seq = hum_dict[peptide][3]
        match_ps = hum_dict[peptide][4]/tum_ps
	
        if match_seq == peptide:
            match_bd = "0"
            match_bind = tum_bind
        elif match_seq == normal_pep:
            match_bd = binding_difference
            match_bind = norm_bind
        else:
            allele = line[1]
            adj_allele = "HLA-" + allele
            mhc_fasta = mhcdir + match_seq + ".fasta"
            mhc_out = mhcdir + match_seq + "." + allele + ".out"
            if match_seq not in affinity_dict:
                fasta_fh = open(mhc_fasta, "w")
                fasta_fh.write("> " + match_seq + "\n")
                fasta_fh.write(match_seq + "\n")
                fasta_fh.close()
                subprocess.call(["netMHCpan", "-a", adj_allele, "-l", str(len(match_seq)), "-inptype", "0", "-xls", "-xlsfile", mhc_out, "-f", mhc_fasta])
                mhc_fh = open(mhc_out, "r")
                mhc_fh.readline()
                mhc_fh.readline()
                mhc_result = mhc_fh.readline().split("\t")
                mhc_fh.close()
                mhc_score = mhc_result[4]
                affinity_dict[match_seq] = {allele : mhc_score}
            else:
                entry = affinity_dict[match_seq]
                if allele not in entry:
                    subprocess.call(["netMHCpan", "-a", adj_allele, "-l", str(len(match_seq)), "-inptype", "0", "-xls", "-xlsfile", mhc_out, mhc_fasta])
                    mhc_fh = open(mhc_out, "r")
                    mhc_fh.readline()
                    mhc_fh.readline()
                    mhc_result = mhc_fh.readline().split("\t")
                    mhc_fh.close()
                    mhc_score = mhc_result[4]
                    affinity_dict[match_seq][allele] = mhc_score	
                else:
                    mhc_score = affinity_dict[match_seq][allele]
            match_bind = mhc_score
            match_bd = float(mhc_score) - float(tum_bind)
    else:
        blast_match = "NA"
        blast_score = "NA"
        match_seq = "NA"
        match_ps = "NA"
        match_bind = "NA"
        match_bd = "NA"
        match_stat = "NA"
    
    if peptide in bac_dict:
        bac_score = bac_dict[peptide][0]
        bac_match = bac_dict[peptide][1]
        bac_seq = bac_dict[peptide][2]
        bac_ps = bac_dict[peptide][3]/tum_ps
    else:
        bac_score = "NA"
        bac_match = "NA"
        bac_seq = "NA"
        bac_ps = "NA"

    if peptide in vir_dict:
        vir_score = vir_dict[peptide][0]
        vir_match = vir_dict[peptide][1]
        vir_seq = vir_dict[peptide][2]
        vir_ps = vir_dict[peptide][3]/tum_ps
    else:
        vir_score = "NA"
        vir_match = "NA"
        vir_seq = "NA"
        vir_ps = "NA"

    line = "\t".join(line)
    line_list = [line, str(binding_difference), str(peptide_similarity), stat, blast_match_trans, blast_match_gene, str(blast_score), match_seq, str(match_bind), str(match_bd), str(match_ps), str(bac_score), bac_match, bac_seq, str(bac_ps), str(vir_score), vir_match, vir_seq, str(vir_ps)]
    newline = "\t".join(line_list) + "\n"
    out_fh.write(newline)

fh1.close()
out_fh.close()

subprocess.call(["rm", "-r", mhcdir])
