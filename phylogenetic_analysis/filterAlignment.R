library(tidyverse)
library(rtracklayer)

args=commandArgs(T)




aln.file=args[1]
aln.obj=seqinr::read.alignment(aln.file,format="fasta")
aln=aln.obj$seq %>% unlist
names(aln)=aln.obj$nam

cutoff=args[2]


outfile=args[3]

seqs.v=aln %>% names

d=vector(mode="list",length(seqs.v))
names(d)=seqs.v

for(s in seqs.v){
	print(paste0("Calculating coverage for ",as.character(s),"... "))
	s.DNAstr=aln[[as.character(s)]]
	#aln length
	len.total=s.DNAstr %>% as.character %>% strsplit(split="") %>% .[[1]] %>% length
	#aln length wo gaps
	len.aligned=s.DNAstr %>% as.character %>% strsplit(split="") %>% .[[1]] %>% .[.!="-"] %>% length
	#print(paste0("Total alignment length: ", len.total))
	#print(paste0("Number of aligned bases: ", len.aligned))
	p.cov=len.aligned*100/len.total
	print(paste0("Coverage over alignment: ", p.cov))
	d[[s]]=p.cov

}

d.filtered=d[d>cutoff]
d.filtered.names=d.filtered %>% names
aln.filtered=aln[d.filtered.names]

write.table(paste0(">",names(aln.filtered),"\n",toupper(aln.filtered),sep=""),file=as.character(outfile), quote=F,row.names=F,col.names=F)
#export(aln.filtered,con=as.character(outfile),format="fasta")

