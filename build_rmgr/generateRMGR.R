library(tidyverse)
library(GenomicRanges)


args=commandArgs(T)


bed=read.table(as.character(args[1]),sep="\t",fill=T) %>% tibble

eq=read.table(as.character(args[2]),sep="\t",fill=T) %>% tibble
colnames(eq)=c("TEfam","TEorder","TEsuperfam")
outprefix=as.character(args[3])

bed=bed %>% left_join(.,eq,by=c("V4"="TEfam"))
bed=bed %>% mutate(ID=paste0(outprefix,"_TE",seq(1,NROW(bed))), TEfam=V4) %>% dplyr::select(-V4) %>% mutate(strand=ifelse(V6!="+" & V6!="-","*",V6))

rmgr=GRanges(seqnames=bed$V1, ranges=IRanges(start=bed$V2, end=bed$V3,names=bed$ID),strand=bed$strand,Kimura=bed$V7, TEfam=bed$TEfam, TEsuperfam=bed$TEsuperfam, TEorder=bed$TEorder)
save(rmgr,file=paste0(outprefix,".rmgr.RData"))
