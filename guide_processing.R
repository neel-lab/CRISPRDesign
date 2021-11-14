library(dplyr);library(reshape2)
#Input argument processing:
args=commandArgs(trailing=T)
#1. Potency Table: Tab-delimited output from "run_azimuth.py"
pot_table=args[1]
#2. User specified number of top guides, based on potency score
# Could be set to "Total" to just include all guides:
#top_n_val=args[2] %>% {ifelse(.=='Total',as.character(.),as.integer(.))}
top_n_val=ifelse(args[2]=='Total',args[2],as.integer(args[2]))
print(top_n_val)
f_pref=unlist(strsplit(split='CRISPR\\_Seqs\\_potency\\.',x=pot_table))[1]

#Read in the potency table:
tp_guidePotency<-read.table(file=pot_table,sep='\t',col.names=c('TxID','Context_Sequence','Start','End','AATarget','frac_protein','strand','potency')) %>% 
	mutate(Gene=sapply(as.character(TxID),function(x) unlist(strsplit(split='\\_',x=x))[2]))
#Filtering the data:
# Main datatable has transcript guide efficacy as well. 
# Will filter only for distinct entries of Gene-Guide pairs, then will filter based on potency to pick top 10 guides.

mainTab<-tp_guidePotency %>% 
	#Only select the Gene and Guide Sequence information, then find unique entries:
	select(Gene,Context_Sequence,potency) %>% distinct() %>% 
	# For every gene in the table, pick the top 10 guides based on potency score:
	group_by(Gene) %>% 
	top_n(wt=potency,top_n_val) %>%
	# Filter based on the criterion:
	#{ifelse(top_n_val=='Total',.,top_n(.,wt=potency,top_n_val))} %>% 
	# Parse the sgRNA sequence from the Context Sequence column.
	mutate(sgRNA=sapply(as.character(Context_Sequence),function(x) substr(x,5,24))) %>% 
	# Reorder the data in the table:
	select(Gene,sgRNA,Context_Sequence,potency) %>% 
	ungroup(Gene)

#Helper function to create lists of information:
aggDta<-function(gene,ctxt,tbl){
	return(tbl %>% filter(Gene==gene,Context_Sequence==ctxt))
}
# For all top 10 guides, find the following:
#1. All Transcripts that the guide hits (TxID_List)
#2. Start Position on Transcript (StartList)
#3. End Position on Transcript (EndList)
#4. Fraction from start where protein is cleaved (frac_protein_list)
#5. List of potencies for this guide on all targeted transcripts (Tx_potency_list)
#NOTE: This list of potencies will contain the main potency listed for the guide

mainTab_extras<-mainTab %>% mutate(
	#1.
	TxID_List=apply(mainTab,1,function(r){
		paste(aggDta(r['Gene'],r['Context_Sequence'],tp_guidePotency) %>% select(TxID) %>% as.matrix() %>% as.character(),collapse=';')
	}),
	#2.
	StartList=apply(mainTab,1,function(r){
		paste(aggDta(r['Gene'],r['Context_Sequence'],tp_guidePotency) %>% select(Start) %>% as.matrix() %>% as.character(),collapse=';')
	}),
	#3.
	EndList=apply(mainTab,1,function(r){
		paste(aggDta(r['Gene'],r['Context_Sequence'],tp_guidePotency) %>% select(End) %>% as.matrix() %>% as.character(),collapse=';')
	}),
	#4. 
	frac_protein_list=apply(mainTab,1,function(r){
		paste(aggDta(r['Gene'],r['Context_Sequence'],tp_guidePotency) %>% select(frac_protein) %>% as.matrix() %>% as.character(),collapse=';')
	}),
	#5. 
	Tx_potency_list=apply(mainTab,1,function(r){
		paste(aggDta(r['Gene'],r['Context_Sequence'],tp_guidePotency) %>% select(potency) %>% as.matrix() %>% as.character(),collapse=';')
	})
)
#Save the information into a table:
write.table(file=paste(f_pref,'_guideTable_topPotency_info.txt',sep=''),quote=F,row.names=F,sep='\t',mainTab_extras)
