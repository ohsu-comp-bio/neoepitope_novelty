# Neoepitope novelty

These tools provide annotation of neoepitope predictions to allow for insights into the level of novelty of a neoepitope.
We provide metrics regarding novelty of MHC binding and of peptide sequence relative to human, bacterial, and viral peptides.
Our manuscript in NAME/LINK demonstrates the utility of investigating these metrics for predicting immunogenic peptides.

# Steps for use:

## 1) Download required software:

- BLAST+ (v2.4.0+)
- NetMHCpan (v2.8)
- Python (v2.7+)
- Python libraries: BioPython (v1.70)

You will also need a neoepitope prediction software!
We used and recommend pVAC-Seq (v4.0.8) for this analysis, as this code best supports that functionality.
However, you're welcome to use any software you like, it will just require a bit more processing on your end!


## 2) Obtain and process required files:

a) Download and unzip Refseq's nonredundant bacterial peptide fasta set:
ftp://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/bacteria.nonredundant_protein.*.protein.faa.gz

Note: our code works for release #83

Concatenate all the bacterial fasta files into one combined, bacterial peptide fasta file

b) Run prepare_dbs.py from within the neoepitope_novelty/ directory to produce a blast protein database from the combined bacterial fasta file:
	
`python prepare_dbs.py -p <path to makeblastdb executable> -d blast_dbs/ -b <path to bacterial fasta>`

c) Download the pickled dictionary of bacterial peptides, bacterialDict.pickle, available on figshare:
https://figshare.com/s/c1094f765bf874bfc4ac
Move the dictionary file into the neoepitope_novelty/dictionaries/ directory


## 3) For each sample and allele, run pVAC-Seq (v4.0.8) according to their recommendations:
http://pvac-seq.readthedocs.io/en/v4.0.8/index.html

You may choose any epitope length of your choice, but otherwise use default parameters.

You may opt to use different parameters or a different tool, but you will then need to parse your results independently.
(See step 4 for instructions.)


## 4) For each run of pVAC-Seq, parse the results to produce a reduced results file:

If pVAC-Seq was run according to our recommendations above, you can use parse_pvac.py to produce the formatted file:

`python parse_pvac.py -i <input file> -o <output file> -s <sample name>`

-i input file, the path to the results from running pVAC-Seq with the recommended settings

-o output file, the desired path to write results to

-s sample name to distinguish run


If you chose to use different parameters for running pVAC-Seq, or used a different program, that's fine too!
You'll just need to parse the results independently.
Produce a tab-delimited file from your results with the following columns in this order and no header:

`<allele> <neoepitope sequence> <neoepitope binding affinity> <normal epitope sequence> <normal epitope binding affinity> <transcript> <gene>
 

## 5) For each sample, run annotate_neoepitopes.py:

`python annotate_neoepitopes.py -i <input file> -o <output directory> -s <sample name> -a <allele> -p <blastp path> -n <netMHCpan path>`

-i	input file to program, the path to the result of parsing a pVAC-Seq results file with parse_pvac.py

-o	path to output directory to which output files will be written

-s	sample name to distinguish run

-a	HLA allele to use for analysis (format example: HLA-A02:01)

-p	path to blastp+ executable

-n	path to netMHCpan executable


Results will be in your chosen output directory under the name [SAMPLE NAME].epitopes.annotated.tsv


## 6) If you use these tools in analyses for your publication, please cite our manuscript:
CITATION HERE