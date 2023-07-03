library(tidyverse)
library(GenomicFeatures)
library(rtracklayer)

args=commandArgs(T)

load(args[1]) # load rmgr
rmgr.tibble=rmgr %>% as.data.frame() %>% tibble()
#txdb=makeTxDbFromGFF(args[2])
ann=import(args[2])
ann=ann[ann$type=="mRNA"]
genes.gr=ann
#genes.gr=genes(txdb)

TEsuperfams=args[3] %>% as.character()
TEsuperfams=strsplit(TEsuperfams,split=",")[[1]]

outprefix=args[5]
library(grid)
library(gridExtra)
reduceOverlappingTEs=as.logical(args[4])
fileConn<-file(paste0(outprefix,".summary.txt"))
pdf(paste0(outprefix,".pdf"),width=12,height=4)
s=vector(mode="list",length(TEsuperfams))
names(s)=TEsuperfams
for(TEsf in TEsuperfams){
	TEsf.gr=rmgr.tibble %>% filter(TEsuperfam == TEsf) %>% GRanges

	if(any(reduceOverlappingTEs)){
		warning("Note: Since `reduceOverlappingTEs` is set to TRUE overlapping TE fragments will be merged")
		TEsf.gr=GenomicRanges::reduce(TEsf.gr,min.gapwidth=2L)
	}
	else {
		warning("Note: Since `reduceOverlappingTEs` is set to FALSE overlapping TE fragments will NOT be merged")
	}

	# general metrics for the superfamily
	totalNfrags=length(TEsf.gr)
	print(paste0("calculating overlaps of ", TEsf,"..."))
	totalBases=TEsf.gr %>% as.data.frame %>% tibble %>% pull(width) %>% sum

	o=findOverlapPairs(TEsf.gr,genes.gr,type="within")
	overlappingNfrags=o@first %>% length
	percOverlappingFrags=overlappingNfrags*100/totalNfrags
	overlappingBases=o@first %>% as.data.frame %>% tibble %>% pull(width) %>% sum
	percOverlappingBases=overlappingBases*100/totalBases
	
	nGenesOverlapping=o@second$ID %>% unique %>% length

	strperc=paste0("% of ",TEsf," fragments overlapping genes ",percOverlappingFrags," %\n% of ",TEsf," bases overlapping genes ",percOverlappingBases," %\nNumber of genes overlapping ",TEsf," fragments ", nGenesOverlapping,"\n")
	print(strperc)
	s[TEsf]=strperc

	if(nGenesOverlapping!=0){
	# calc number of fragments per gene
	fragsOnGenes=o@second$ID %>% tibble %>% mutate(.,ID=`.`) %>% group_by(ID) %>% summarise(count=n()) %>% mutate(log2Count=log2(count))
	fragsOnGenes=fragsOnGenes %>% mutate(`% overlapping fragments on this gene`=count*100/overlappingNfrags, `TE superfamily`=TEsf)
	
	idAnno.eq=o@second %>% as_tibble() %>% pull(product,ID) %>% unlist(recursive = T,use.names=T) %>% tibble(ID=names(.),Product=.)
	idAnno.eq=idAnno.eq %>% group_by(ID) %>% summarise(description=c(paste0(`Product`,"; "))) %>% distinct
	idAnno.eq.v=idAnno.eq %>% pull(description,ID)


	# calc number of bases per gene
	genes.tibble=o@second %>% data.frame %>% tibble %>% dplyr::select(ID)
	tes.tibble=o@first %>% data.frame %>% tibble
	geneNames=o@second$ID %>% unique 
	geneLengths=genes.gr %>% as.data.frame %>% filter(ID %in% geneNames) %>% pull(width,ID)
	geneOverlappingMetrics.tibble=cbind(genes.tibble,tes.tibble) %>% tibble %>% group_by(ID) %>% summarise(NBPS=sum(width)) %>% mutate(`gene length`=geneLengths[ID]) %>% mutate(percOfGeneLengthOverlappingTEsuperfam=NBPS*100/`gene length`)

	fragsOnGenes=left_join(fragsOnGenes,geneOverlappingMetrics.tibble,by="ID")

	# calc stats
	ngenes=NROW(fragsOnGenes)
	avg=fragsOnGenes %>% pull(count) %>% mean %>% format(digits=3)
	std=fragsOnGenes %>% pull(count) %>% sd %>% format(digits=3)
	med=fragsOnGenes %>% pull(count) %>% median %>% format(digits=3)
	mx=fragsOnGenes %>% pull(count) %>% max %>% format(digits=3)

	statsText=paste0("N = ",ngenes,"\n",
			 "Mean = ",avg,"\n",
			 "SD = ",std,"\n",
			 "Median = ",med,"\n",
			 "Max = ",mx)

	# plot
	p=fragsOnGenes %>% ggplot(aes(x=count)) + geom_density() + 
		ggtitle(TEsf)+theme_classic()
	p.hist=fragsOnGenes %>% ggplot(aes(x=count)) + geom_histogram() + 
		ggtitle(TEsf)+theme_classic()
	plog=fragsOnGenes %>% ggplot(aes(x=log2Count)) + geom_density() + 
		ggtitle(TEsf)+theme_classic()
	p.txt=ggplot()+annotate("text",x=1,y=1,size=4,label=statsText)+theme_void()
	grid.arrange(p,plog,p.txt,ncol=3)
	

	# tables
	fragsOnGenes %>% arrange(desc(percOfGeneLengthOverlappingTEsuperfam)) %>% dplyr::select(-log2Count) %>% mutate(Products=idAnno.eq.v[ID]) %>% write.table(paste0(outprefix,".",TEsf,".tsv"),sep="\t",col.names=T,row.names=F,quote=F)
#	geneOverlappingMetrics.tibble %>% arrange(desc(percOfGeneLengthOverlappingTEsuperfam)) %>%  write.table(paste0(outprefix,".",TEsf,".geneMetrics.tsv"),sep="\t",col.names=F,row.names=F,quote=F)
	}
	else{
		warning(paste0("Not overlaps found for ",TEsf))
	}
}
writeLines(unlist(s), fileConn)
close(fileConn)
dev.off()
save.image(paste0(outprefix,".RData"))

