library(tidyverse)
library(GenomicFeatures)
library(rtracklayer)

args=commandArgs(T)

load(args[1]) # load rmgr
rmgr.tibble=rmgr %>% as.data.frame() %>% tibble()
ann=import(args[2])
ann=ann[ann$type=="mRNA"]
genes.gr=ann

regionSize=as.numeric(args[3])

TEsuperfams=args[4] %>% as.character()
TEsuperfams=strsplit(TEsuperfams,split=",")[[1]]

reduceOverlappingTEs=as.logical(args[5])

outprefix=args[6]
library(grid)
library(gridExtra)


for(TEsf in TEsuperfams){
	TEsf.gr=rmgr.tibble %>% filter(TEsuperfam == TEsf) %>% GRanges

	if(any(reduceOverlappingTEs)){
		warning("Note: Since `reduceOverlappingTEs` is set to TRUE overlapping TE fragments will be merged")
		TEsf.gr=GenomicRanges::reduce(TEsf.gr,min.gapwidth=2L)
		TEsf.gr$TE_ID=paste0(TEsf,"|",seq(1:length(TEsf.gr)))
	}
	else {
		warning("Note: Since `reduceOverlappingTEs` is set to FALSE overlapping TE fragments will NOT be merged")
	}

	o=findOverlapPairs(TEsf.gr,genes.gr,type="within")
	
	nGenesOverlapping=o@second$ID %>% unique %>% length
	if(nGenesOverlapping!=0){
		print(paste0("overlapsFound for ",TEsf))	
		TEsf.gr.o=o@first
		TEsf.gr.o$gene_ID=o@second$ID
		coords.tibble=TEsf.gr.o %>% as.data.frame %>% tibble %>% filter(width>1000) %>% mutate(name=paste0(TE_ID,"_",gene_ID),score=".",phase=".")
		coords.bed=coords.tibble %>% dplyr::select(seqnames,start,end,name,strand,score,phase)
		write.table(coords.bed,file=paste0(outprefix,".",TEsf,".bed"),sep="\t",quote=F,row.names=F, col.names=F)

	}
	else{
		warning(paste0("No overlaps found for ",TEsf))
	}
}


print("All BEDs generated!")
#dev.off()
#save.image(paste0(outprefix,".inspectTandems.RData"))

