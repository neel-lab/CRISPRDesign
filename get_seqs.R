library(biomaRt);library(dplyr);library(reshape2)
#Input should be a list containing all library genes:
args=commandArgs(trailing=T)
geneList_file=args[1] #a file containing all genes for library:
f_pref=unlist(strsplit(split='\\.',x=geneList_file))[1]

#Link to mart:
ensembl<-useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

#Read in the transporter gene names:
transporters<-read.table(file=geneList_file,sep=',',header=F) %>% as.matrix() %>% as.character()

#Get a list of transcript versions:
message("Gathering all transcripts for each gene...")
tx_versions<-getBM(attributes=c('ensembl_transcript_id','hgnc_symbol'),filters='hgnc_symbol',values=transporters,mart=ensembl)
#Which genes weren't found in the database?:
nofind<-setdiff(transporters,unique(tx_versions$hgnc_symbol))
write.table(file=paste(f_pref,'couldnt_find.tsv',sep='_'),sep='\t',quote=F,row.names=F,nofind)
message(paste("Couldn't Find",length(nofind),"gene sequences in the query."))

#Get sequences of transcript variants and merge:
message("Gathering sequences...")
seqs = getSequence(id=tx_versions$ensembl_transcript_id,type='ensembl_transcript_id',seqType='coding',mart=ensembl)

#Merge sequence and gene name information together, and reorganize:
tx_seq_dta<-seqs %>% full_join(tx_versions,by='ensembl_transcript_id') %>% select(coding,ensembl_transcript_id,hgnc_symbol)

#Write everything out as a FASTA file:
exportFASTA(tx_seq_dta,paste(f_pref,'coding_seqs.fa',sep='_'))
