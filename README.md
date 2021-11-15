#CRISPR Design Pipeline

This pipeline takes a list of HGNC names in a text file and outputs a table of the top `n` most potent guides as determined by the [Azimuth Algorithm][https://github.com/MicrosoftResearch/Azimuth/tree/v2.0]

All intermediate files generated from this pipeline have the same prefix.  For example, if the name of a file is `crispr_library_genes.txt`, then all intermediate files have the prefix `crispr_library_genes` prepended.  Intermediate files are explained in the pipeline steps below.

## Running the pipeline:

The pipeline can be executed by running the following at the command line:

```
bash run_CRISPR_design.sh ${GENELIST_FILE} ${TOP_N_GUIDES}
```

Where `${GENELIST_FILE}` is a file containing HGNC gene names, one name per row, and `${TOP_N_GUIDES}` is an integer specifying the number of guides per gene to keep, based on the potency.  The final designs are kept in a file called `PREFIX__guideTable_topPotency_info.txt`, where PREFIX is the input file name to the pipeline.

## Pipeline Steps:

### Step 1: Gene Sequence Gathering:

The HGNC names are passed to an R script which reaches out to Ensembl biomart and retrieves every spliced coding transcript of a gene.  These sequences are stored in the `..._coding_seqs.fa` file.  Any genes that couldn't be found will be written to a `..._couldnt_find.tsv` file.  The pipeline will also print out the number of genes that coudn't be found.

### Step 2: Finding CRISPR Guide Candidates: 

CRISPR sequences are parsed from each gene's set of transcripts using the `get_guides_seqs.awk` script.  The Azimuth program requires 3 input features to provide guide potencies: 1. 30-mer guide sequence consisting of 5 bases before the protospacer, the 20-base protospacer, the NGG PAM, and 2 bases after the PAM. 2. The peptide index affected by the guide, and 3. a fraction representing how close to the beginning of a peptide the guide cleaves (0 being the N terminus, 1 being the C terminus).  The cut site is defined as the 16th base within the protospacer.  Peptide index is computed by taking the cut site index and dividing by 3.  The percent peptide is computed by dividing the peptide index by the total number of peptides, which is the length of the coding transcript divided by 3.

The awk script outputs a table called `..._CRISPR_seqs.txt`, with the following columns:

1. Transcript ID _ GeneName
2. 30-mer guide sequence
3. Start position of 30-mer on transcript
4. End position of 30-mer on transcript
5. Peptide index
6. Fraction into peptide (values between 0 and 1)
7. Forward or Reverse Strand (F/R)


### Step 3. Calling Azimuth:

The table containing guide data is passed into Azimuth via a python script which computes the potency information based on a gradient boosted regression tree model.  This potency information is added into the input table and is written out to another file called `..._CRISPR_seqs_potency.txt`

### Step 4. Filtering Top N Most Potent Guides:

An R script called `guide_processing.R` reads the output from the Azimuth script and filters the guides for each gene based on the user-provided threshold to the pipeline. It outputs a table with the following columns:  

1. Gene	: Name of the gene targeted 
2. sgRNA : 20-mer protospacer sequence of guide
3. Context_Sequence: 30-mer input to Azimuth
4. potency: Potency computed by Azimuth
5. TxID_List: Transcripts targeted by this guide
6. StartList: Starting index where guide hits each transcript (in order of TxID_List)
7. EndList: Ending index where guide hits each transcript (in order of TxID_List)
8. frac_protein_list: Fraction of peptide targeted for each transcript (in order of TxID_List)
9. Tx_potency_list: Potency of guide against each transcript (in order of TxID_List)


## Dependencies:

R and python dependencies are stored in `Rpacks` and `CRISPRDesignVenv`, respectively.  The pipeline script should run if executed within this directory.

## Tips for Execution:

- Many sgRNA sequences can be generated when parsing from transcripts.  You can uncomment parts of the `get_guides_seqs.awk` which allows for a thresholding on how many guides are returned.  If not restricted, sgRNA candidate files can be a few megabytes, and can cause Azimuth to take a while inferring the potency of guides.




