library(tidyverse)
library(GenomicFeatures)
library(rtracklayer)

args=commandArgs(T)

load(args[1]) # load rmgr
rmgr.tibble=rmgr %>% as.data.frame() %>% tibble()
txdb=makeTxDbFromGFF(args[2])
exons.gr=exons(txdb) %>% reduce
introns.gr=intronsByTranscript(txdb) %>% reduce
flank_size=500
ann=import(args[2])
ann=ann[ann$type=="mRNA"]
genes.gr=ann
fiveUTR.gr=flank(genes.gr, flank_size)
threeUTR.gr=flank(genes.gr, flank_size, start=FALSE)
#genes.gr=genes(txdb)

# total bases for subregions
allgenesbps=genes.gr %>% GenomicRanges::reduce(.) %>% as.data.frame %>% tibble %>% pull(width) %>% sum
allexonsbps=exons.gr %>% GenomicRanges::reduce(.) %>% as.data.frame %>% tibble %>% pull(width) %>% sum
allintronsbps=introns.gr %>% GenomicRanges::reduce(.) %>% as.data.frame %>% tibble %>% pull(width) %>% sum
allfiveUTRbps=fiveUTR.gr %>% GenomicRanges::reduce(.) %>% as.data.frame %>% tibble %>% pull(width) %>% sum
allthreeUTRbps=threeUTR.gr %>% GenomicRanges::reduce(.) %>% as.data.frame %>% tibble %>% pull(width) %>% sum

TEsuperfams=args[4] %>% as.character()
TEsuperfams=strsplit(TEsuperfams,split=",")[[1]]
colourPalette=read.table(args[3], sep="\t", comment.char="@",col.names=c("TEsuperfam","colour")) %>% tibble %>% pull(colour,TEsuperfam)

outprefix=args[7]
library(grid)
library(gridExtra)
reduceOverlappingTEs=as.logical(args[5])
textTag=as.character(args[6])
fileConn<-file(paste0(outprefix,".summary.txt"))
pdf(paste0(outprefix,".pdf"),width=12,height=4)
s=vector(mode="list",length(TEsuperfams))
names(s)=TEsuperfams
pdf.lst=vector(mode="list",length(TEsuperfams)) # probability density function
names(pdf.lst)=TEsuperfams

percOverlappingTEs.v=rep(0,length(TEsuperfams))
names(percOverlappingTEs.v)=TEsuperfams
percOverlappingbps.v=rep(0,length(TEsuperfams))
names(percOverlappingbps.v)=TEsuperfams
nOverlappingGenes.v=rep(0,length(TEsuperfams))
names(nOverlappingGenes.v)=TEsuperfams
exons.bps.v=rep(0,length(TEsuperfams))
names(exons.bps.v)=TEsuperfams
introns.bps.v=rep(0,length(TEsuperfams))
names(introns.bps.v)=TEsuperfams
fiveUTR.bps.v=rep(0,length(TEsuperfams))
names(fiveUTR.bps.v)=TEsuperfams
threeUTR.bps.v=rep(0,length(TEsuperfams))
names(threeUTR.bps.v)=TEsuperfams
genes.bps.v=rep(0,length(TEsuperfams))
names(genes.bps.v)=TEsuperfams


if(textTag=="custom"){

	aliasDesc=	ann %>% as.data.frame %>% tibble %>% dplyr::select(protein_match,alias) %>% unnest %>% distinct %>%
		mutate(size=nchar(alias)) %>% group_by(protein_match) %>% filter(size>1) %>% filter(size==min(size)) %>% ungroup %>% pull(alias,protein_match)
	productDesc=	ann %>% as.data.frame %>% tibble %>% dplyr::select(protein_match,product) %>% unnest %>% mutate(size=nchar(product)) %>% group_by(protein_match) %>% filter(size==min(size)) %>% ungroup %>% pull(product,protein_match)


	description.v=ann %>% as.data.frame %>% tibble %>% dplyr::select(ID,protein_match) %>% mutate(description=ifelse(is.na(aliasDesc[protein_match])==T,ifelse(is.na(productDesc[protein_match])==T," ",productDesc[protein_match]),aliasDesc[protein_match])) %>% pull(description,ID)
}
if(textTag=="ncbi"){
	description.v=ann %>% as.data.frame %>% tibble %>% dplyr::select(gene,ID) %>% mutate(description=ifelse(is.na(gene)==T," ",gene)) %>% pull(description,ID)
}


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
	percOverlappingTEs.v[TEsf]=percOverlappingFrags
	overlappingBases=o@first %>% as.data.frame %>% tibble %>% pull(width) %>% sum
	percOverlappingBases=overlappingBases*100/totalBases
	percOverlappingbps.v[TEsf]=percOverlappingBases
	
	nGenesOverlapping=o@second$ID %>% unique %>% length
	nOverlappingGenes.v[TEsf]=nGenesOverlapping
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
	p=fragsOnGenes %>% ggplot(aes(x=count,fill=`TE superfamily`)) + geom_histogram(binwidth=1) + 
		ggtitle(TEsf)+theme_classic() + xlab("Number of fragments") + ylab("Density") + scale_fill_manual(values=colourPalette)+theme(legend.position="none")+xlim(0,50)
	p.hist=fragsOnGenes %>% ggplot(aes(x=count)) + geom_histogram() + 
		ggtitle(TEsf)+theme_classic() + xlab("Number of fragments") + ylab("Density")
	plog=fragsOnGenes %>% ggplot(aes(x=log2Count)) + geom_density() + 
		ggtitle(TEsf)+theme_classic() + xlab("Number of fragments") + ylab("Density")
	p.txt=ggplot()+annotate("text",x=1,y=1,size=3,label=statsText)+theme_void()
	pdf.lst[[TEsf]]=grid.arrange(p,p.txt,ncol=1)
	

	# tables
	fragsOnGenes %>% arrange(desc(percOfGeneLengthOverlappingTEsuperfam)) %>% dplyr::select(-log2Count) %>% mutate(Products=description.v[ID]) %>% write.table(paste0(outprefix,".",TEsf,".tsv"),sep="\t",col.names=T,row.names=F,quote=F)
#	geneOverlappingMetrics.tibble %>% arrange(desc(percOfGeneLengthOverlappingTEsuperfam)) %>%  write.table(paste0(outprefix,".",TEsf,".geneMetrics.tsv"),sep="\t",col.names=F,row.names=F,quote=F)


	# overlapping with exons and introns
	o.exons=findOverlapPairs(TEsf.gr,exons.gr,type="within")
	exons.bps=o.exons@first %>% as.data.frame %>% tibble %>% pull(width) %>% sum 
	o.introns=findOverlapPairs(TEsf.gr,introns.gr,type="within")
	introns.bps=o.introns@first %>% as.data.frame %>% tibble %>% pull(width) %>% sum 

	o.fiveUTR=findOverlapPairs(TEsf.gr,fiveUTR.gr,type="within")
	fiveUTR.bps=o.fiveUTR@first %>% as.data.frame %>% tibble %>% pull(width) %>% sum 
	o.threeUTR=findOverlapPairs(TEsf.gr,threeUTR.gr,type="within")
	threeUTR.bps=o.threeUTR@first %>% as.data.frame %>% tibble %>% pull(width) %>% sum 
	# overlappingBases

	exons.bps.v[TEsf]=exons.bps*100/allexonsbps #overlappingBases #totalBases
	introns.bps.v[TEsf]=introns.bps*100/allintronsbps #overlappingBases #totalBases
	fiveUTR.bps.v[TEsf]=fiveUTR.bps*100/allfiveUTRbps #overlappingBases #totalBases
	threeUTR.bps.v[TEsf]=threeUTR.bps*100/allthreeUTRbps #overlappingBases #totalBases

	genes.bps.v[TEsf]=overlappingBases*100/allgenesbps
	}
	else{
		warning(paste0("Not overlaps found for ",TEsf))
	}
}
writeLines(unlist(s), fileConn)
close(fileConn)

percs.tibble=tibble(`TE superfamily`=names(percOverlappingTEs.v),
		    `% of TE fragments overlapping genes`=percOverlappingTEs.v,
		    `% of TE bps overlapping genes`=percOverlappingbps.v,
		    `Number of overlapping genes`=nOverlappingGenes.v,
		    `% of bps overlapping genes`=genes.bps.v,
		    `% of bps overlapping exons`=exons.bps.v,
		    `% of bps overlapping introns`=introns.bps.v,
		    `% of bps overlapping 5'-UTR`=fiveUTR.bps.v,
		    `% of bps overlapping 3'-UTR`=threeUTR.bps.v)

percs.tibble$`TE superfamily`=factor(percs.tibble$`TE superfamily`,levels=TEsuperfams)

#p.percTEs=percs.tibble %>% ggplot(aes(x=`TE superfamily`,y=`% of fragments overlapping genes`))+geom_bar(stat="identity")+theme_classic()+theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),axis.title.x=element_blank())
p.nGenes=percs.tibble %>% ggplot(aes(x=`TE superfamily`,y=`Number of overlapping genes`,fill=`TE superfamily`))+geom_bar(stat="identity",width=0.8)+theme_classic()+theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),axis.title.x=element_blank(),legend.position="none")+ggtitle("Genes overlapping TEs")+ylab("Number of genes")+scale_fill_manual(values=colourPalette)
p.genes=percs.tibble %>% ggplot(aes(x=`TE superfamily`,y=`% of bps overlapping genes`,fill=`TE superfamily`))+geom_bar(stat="identity",width=0.8)+theme_classic()+theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),axis.title.x=element_blank(),legend.position="bottom")+ggtitle("Percentage of the genes length")+ylab("% of bp")+scale_fill_manual(values=colourPalette)

g_legend<-function(a.gplot){
        tmp <- ggplot_gtable(ggplot_build(a.gplot))
        leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
        legend <- tmp$grobs[[leg]]
        return(legend)
}

p.legend=g_legend(p.genes)
p.genes=p.genes+theme(legend.position="none")
p.exons=percs.tibble %>% ggplot(aes(x=`TE superfamily`,y=`% of bps overlapping exons`,fill=`TE superfamily`))+geom_bar(stat="identity",width=0.8)+theme_classic()+theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),axis.title.x=element_blank(),legend.position="none")+ggtitle("Exons")+ylab("% of bp")+scale_fill_manual(values=colourPalette)
p.introns=percs.tibble %>% ggplot(aes(x=`TE superfamily`,y=`% of bps overlapping introns`,fill=`TE superfamily`))+geom_bar(stat="identity",width=0.8)+theme_classic()+theme(axis.text.x =element_text(angle = 45, vjust = 0.5, hjust=1),axis.title.x=element_blank(),legend.position="none")+ggtitle("Introns")+ylab("% of bp")+scale_fill_manual(values=colourPalette)
p.fiveUTR=percs.tibble %>% ggplot(aes(x=`TE superfamily`,y=`% of bps overlapping 5'-UTR`,fill=`TE superfamily`))+geom_bar(stat="identity",width=0.8)+theme_classic()+theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),axis.title.x=element_blank(),legend.position="none")+ggtitle(paste0("Upstream (",flank_size," bp)"))+ylab("% of bp")+scale_fill_manual(values=colourPalette)
p.threeUTR=percs.tibble %>% ggplot(aes(x=`TE superfamily`,y=`% of bps overlapping 3'-UTR`,fill=`TE superfamily`))+geom_bar(stat="identity",width=0.8)+theme_classic()+theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),axis.title.x=element_blank(),legend.position="none")+ggtitle(paste0("Downstream (",flank_size," bp)"))+ylab("% of bp")+scale_fill_manual(values=colourPalette)



library(reshape2)
p.void=ggplot() + theme_void()
general.p=grid.arrange(p.nGenes,p.genes,p.void,p.legend,ncol=2)
subregions.p=grid.arrange(p.introns,p.exons,p.fiveUTR,p.threeUTR,ncol=2)
grid.arrange(general.p,subregions.p,ncol=2)
#grid.arrange(p.nGenes,p.genes,p.introns,p.exons,p.fiveUTR,p.threeUTR,ncol=6)

grid.arrange(grobs=pdf.lst,ncol=length(pdf.lst))
#grid.arrange(general.p,subregions.p,ncol=1)
dev.off()
save.image(paste0(outprefix,".RData"))

