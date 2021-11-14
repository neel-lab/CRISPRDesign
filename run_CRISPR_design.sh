#!/bin/bash
#CRISPR Cas9 Library Design Pipeline:

#Inputs
#Input list: list of genes that need guides:
INPUT_LIST=$1
#Top N: Top sgRNAs per gene to include in final library design:
# Can equal "Total" as well if all guides are to be kept:
TOPN=$2

#Scripts Location:
SCRIPTS=~/Documents/Data/BC_CRISPR_Libs/scripts

#Get the input file prefix:
PREFIX=$(echo $INPUT_LIST | cut -d'.' -f1)
echo "Prefix: $PREFIX"
#Create CRISPR_Seqs Filename:
CRISPRSEQS=$(echo $PREFIX"_CRISPR_Seqs.txt")
echo "CRISPRSeqs Name: $CRISPRSEQS"
#Create Potency Filename:
POTENCY=$PREFIX="_potency.txt"

#1. Retrieve fasta file of sequences for gene list:
# outputs file with suffix "coding_seqs.fa"
echo -e  "1. Retrieving gene list sequences from bioMart\n"
Rscript $SCRIPTS/get_seqs.R $INPUT_LIST
bash $SCRIPTS/edit_fasta.sh $PREFIX"_coding_seqs.fa"
echo -e  "Sequence retrieval Complete\n"

#2. Extract CRISPR sequences:
#Name of output should be file prefix, then "CRISPR_Seqs.txt"
echo -e  "2. Extracting CRISPR Sequences\n"
awk -f $SCRIPTS/get_guides_seqs.awk $PREFIX"_coding_seqs.fa" > $CRISPRSEQS
echo -e  "CRISPR sgRNA Sequence Extraction Complete\n"

#3. Run Azimuth to evaluate guides:
# Returns file called PREFIX"_potency.txt"
echo -e  "3. Grading sgRNA potency using Azimuth...\n"
source activate
python $SCRIPTS/run_azimuth.py $PREFIX"_CRISPR_Seqs.txt"
conda deactivate
echo -e  "Azimuth Complete"

#4. Trim the results to include top 10 guides:
echo -e  "4. Filtering potency data\n"
Rscript $SCRIPTS/guide_processing.R $PREFIX"_CRISPR_Seqs_potency.txt" $TOPN
echo -e  "CRISPR Library Design Complete!"
