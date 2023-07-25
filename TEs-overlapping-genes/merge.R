library(tidyverse)

args=commandArgs(T)


w_reduce.tibble=read.table(args[1],sep="\t",header=T,quote = "",col.names=c("ID","count_overlaps.r","perc_overlapping_frags_on_gene.r","TEsuperfam.r","Number_TE_bps.r","gene_length.r","perc_of_gene_length_overlapping_TEs.r","gene_product.r")) %>% tibble
wo_reduce.tibble=read.table(args[2],sep="\t",header=T,quote = "",col.names=c("ID","count_overlaps.non_r","perc_overlapping_frags_on_gene.non_r","TEsuperfam.non_r","Number_TE_bps.non_r","gene_length.non_r","perc_of_gene_length_overlapping_TEs.non_r","gene_product.non_r")) %>% tibble
#ID      count   % overlapping fragments on this gene    TE superfamily  NBPS    gene length     percOfGeneLengthOverlappingTEsuperfam   Products

tandemsTable=read.table(args[3],sep=" ",fill=T,col.names=c("TEsuperfam","TEsuperfam_frag_ID","Gene","start","end","period_size","copy_number","consensus_size","perc_matches","perc_indels","score","A_composition","C_composition","G_composition","T_composition","entropy","seq1","seq2","seq3","seq4")) %>% tibble
tandemsTable.simple=tandemsTable %>% dplyr::select(TEsuperfam,Gene,copy_number,period_size) %>% group_by(TEsuperfam,Gene) %>% summarize(sum_copy_number=sum(copy_number),avg_period_size=mean(period_size)) %>% ungroup
tandemsTable.merged=tandemsTable.simple %>% left_join(w_reduce.tibble,by=c("TEsuperfam"="TEsuperfam.r","Gene"="ID"))
tandemsTable.merged=tandemsTable.merged %>% mutate(EXT_ID=paste0("TEoverGENE_",seq(1:NROW(tandemsTable.merged))))


TEsuperfams=as.character(args[4])
TEsuperfams=strsplit(TEsuperfams,split=",")[[1]]


perc_gene_length_quantile=as.numeric(args[5]) #0.5
sum_copy_number_quantile=as.numeric(args[6]) #0.99

outprefix=as.character(args[7])


merged.tibble=left_join(wo_reduce.tibble,w_reduce.tibble,by=c("ID"="ID", "TEsuperfam.non_r"="TEsuperfam.r")) %>% dplyr::select(ID,TEsuperfam.non_r,count_overlaps.r,count_overlaps.non_r,perc_overlapping_frags_on_gene.non_r,Number_TE_bps.r,gene_length.r,perc_of_gene_length_overlapping_TEs.r,gene_product.r) %>% mutate(reductionStat=(count_overlaps.non_r-count_overlaps.r)*100/count_overlaps.non_r) %>% filter(count_overlaps.r>0,count_overlaps.non_r>0) %>% arrange(desc(reductionStat))

write.table(merged.tibble,file=paste0(outprefix,".reduced_non-reduced_merged.tsv"),sep="\t",row.names=F,col.names=T,quote=F)

library(grid)
library(gridExtra)
library(ggrepel)
options(ggrepel.max.overlaps = Inf)
p.lst=vector(mode="list",length(TEsuperfams))
names(p.lst)=TEsuperfams


xmax=tandemsTable.merged %>% pull(perc_of_gene_length_overlapping_TEs.r) %>% max
ymax=tandemsTable.merged %>% pull(sum_copy_number) %>% max
sig.v=c()
pdf(paste0(outprefix,".merged.pdf"),width=13,height=8)
for(TEsf in TEsuperfams){


	avg_perc_gene_length=tandemsTable.merged %>% filter(TEsuperfam==TEsf) %>% pull(perc_of_gene_length_overlapping_TEs.r) %>% mean(na.rm=T)
	vline=tandemsTable.merged %>% filter(TEsuperfam==TEsf) %>% pull(perc_of_gene_length_overlapping_TEs.r) %>% quantile(probs=c(perc_gene_length_quantile),na.rm=T)
	hline=tandemsTable.merged %>% filter(TEsuperfam==TEsf) %>% pull(sum_copy_number) %>% quantile(probs=c(sum_copy_number_quantile),na.rm=T)
	avg_sum_copy_number=tandemsTable.merged %>% filter(TEsuperfam==TEsf) %>% pull(sum_copy_number) %>% mean(na.rm=T)
	p.lst[[TEsf]]=tandemsTable.merged %>% 
		dplyr::filter(TEsuperfam==TEsf) %>% 
		mutate(label=ifelse(perc_of_gene_length_overlapping_TEs.r>vline & sum_copy_number>hline, gene_product.r,"")) %>%
		mutate(should_i_color=ifelse(label=="",F,T)) %>%
		ggplot(aes(x=perc_of_gene_length_overlapping_TEs.r,y=sum_copy_number,label=label,colour=should_i_color)) + 
		geom_point(shape=1,alpha=1,size=0.75)+
		geom_text_repel(size=3,min.segment.length=0.25,force=1)+ #size3 solia ser
		geom_vline(xintercept=avg_perc_gene_length,colour="#ca0020",linetype="dashed")+
		geom_vline(xintercept=vline,colour="#f4a582",linetype="dashed")+
		geom_hline(yintercept=avg_sum_copy_number,colour="#0571b0",linetype="dashed")+
		geom_hline(yintercept=hline,colour="#67a9cf",linetype="dashed")+ 
		theme_classic()+ggtitle(TEsf)+theme(legend.position="none")+xlim(0,xmax)+ylim(0,ymax)+xlab("% of gene length")+ylab("Total copy-number")+scale_colour_manual(values=c("grey","black"))


		sig.TE.v=	tandemsTable.merged %>%
                		dplyr::filter(TEsuperfam==TEsf) %>%
                		mutate(significant=ifelse(perc_of_gene_length_overlapping_TEs.r>vline & sum_copy_number>hline, T,F)) %>% pull(significant,EXT_ID)
		sig.v=c(sig.v,sig.TE.v)



}
tandemsTable.merged=tandemsTable.merged %>% mutate(Significant=sig.v[EXT_ID]) %>% arrange(desc(sum_copy_number))
write.table(tandemsTable.merged,file=paste0(outprefix,".reduced_trf_merged.tsv"),sep="\t",row.names=F,col.names=T,quote=F)
grid.arrange(grobs = p.lst, ncol = 3)

dev.off()


save.image(paste0(outprefix,".merged.RData"))
