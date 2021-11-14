#!/usr/bin/awk -f

BEGIN{FS="\n";RS=">";OFS="\t"} 
{
	seqlen=length($2)
	#Forward
	split($2,str,"[ATCG]{25}GG[ATCG]{3}",m)
	for (i in m){
		#Find location of the match:
		match($2,m[i])
		cutsite=(RSTART+21)
		pep_num=int(cutsite/3)
		tot_aa=int(seqlen/3)
		frac_peptide=(pep_num/tot_aa)
		print $1,m[i],RSTART,RSTART+RLENGTH,pep_num,frac_peptide,"F"
	}
	#Reverse
	split($2,str,"[ATCG]{3}CC[ATCG]{25}",m)
	for (i in m){
		#Find location of the match:
		match($2,m[i])
		cutsite=(RSTART+9)
		pep_num=int(cutsite/3)
		tot_aa=int(seqlen/3)
		frac_peptide=(pep_num/tot_aa)
		print $1,m[i],RSTART+RLENGTH,RSTART,pep_num,frac_peptide,"R"
	}
}

