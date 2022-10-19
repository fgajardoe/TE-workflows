library(tidyverse)
library(rtracklayer)
library(GenomicFeatures)

library(ape)
library(geiger)
library(ggtree)

args=commandArgs(T)
treePath=args[1]
#centroids.tibble.RData=args[2]
#spsNames=args[3]
outprefix=as.character(args[2])

args.sps=args[3:length(args)]
labels=args.sps[1:(length(args.sps)/2)]
files=args.sps[(1+length(args.sps)/2):(length(args.sps))]

names(files)=labels

c.lst=vector(mode="list", length(labels))
names(c.lst)=labels

for(sp in labels){
	print(paste0("working on ",sp))
	load(as.character(files[sp]))
	centroids.tibble.sp=centroids.tibble %>% mutate(specie=sp)
	c.lst[[sp]]=centroids.tibble.sp
}
c.allSpecies.tibble=do.call("rbind",c.lst)

c.allSpecies.tibble$perc_fragments_overlapping_txs=as.numeric(c.allSpecies.tibble$perc_fragments_overlapping_txs)
c.allSpecies.tibble$perc_bps_overlapping_txs=as.numeric(c.allSpecies.tibble$perc_bps_overlapping_txs)

#tandemRData=as.character(args[5])


# Tree plotting


tree=read.tree(treePath)
tipsToRemove=c("ORF14_Cynopoecilus_melanotaenia_TE860828/1-98","ORF2_Nematolebias_whitei_TE1960123/1-80")
outgroup="ORF23_Oryzias_latipes_TE89720/1-502"

tree.reduced=drop.tip(tree, tipsToRemove)

tree.rooted=root(tree.reduced,outgroup=outgroup)

centroidTips=tree.rooted$tip.label %>% gsub("ORF\\d+_","",x=.,perl=T) %>% gsub("/\\d+-\\d+","",x=.,perl=T)
tree.centroidTips=tree.rooted
tree.centroidTips$tip.label=centroidTips


# plots

c.allSpecies.tibble=c.allSpecies.tibble %>% filter(centroidID %in% tree.centroidTips$tip.label)
#c.allSpecies.tibble %>% mutate(node=centroidID)

p.nFrag= c.allSpecies.tibble %>% 
	 ggplot(aes(x=n_fragments,y=centroidID,fill=specie)) + geom_bar(stat="identity", width=0.5)+ theme_minimal() + theme(axis.text.y=element_blank(), axis.title.y=element_blank())

p.nBPs= c.allSpecies.tibble %>% 
	 ggplot(aes(x=n_bps,y=centroidID,fill=specie)) + geom_bar(stat="identity", width=0.5)+ theme_minimal() + theme(axis.text.y=element_blank(), axis.title.y=element_blank())


p.nRTs= c.allSpecies.tibble %>% 
	 ggplot(aes(x=n_RT_domains,y=centroidID,fill=specie)) + geom_bar(stat="identity", width=0.5)+ theme_minimal() + theme(axis.text.y=element_blank(), axis.title.y=element_blank())

p.nORFs= c.allSpecies.tibble %>% 
	 ggplot(aes(x=n_ORFs,y=centroidID,fill=specie)) + geom_bar(stat="identity", width=0.5)+ theme_minimal() + theme(axis.text.y=element_blank(), axis.title.y=element_blank())

p.perc_fragments_overlapping_txs= c.allSpecies.tibble %>% 
	 ggplot(aes(x=perc_fragments_overlapping_txs,y=centroidID,fill=specie)) + geom_bar(stat="identity", width=0.5)+ theme_minimal() + theme(axis.text.y=element_blank(), axis.title.y=element_blank())

p.perc_bps_overlapping_txs= c.allSpecies.tibble %>% 
	 ggplot(aes(x=perc_bps_overlapping_txs,y=centroidID,fill=specie)) + geom_bar(stat="identity", width=0.5)+ theme_minimal() + theme(axis.text.y=element_blank(), axis.title.y=element_blank())

p.length= c.allSpecies.tibble %>% 
	 ggplot(aes(x=width,y=centroidID,fill=specie)) + geom_bar(stat="identity", width=0.5)+ theme_minimal() + theme(axis.text.y=element_blank(), axis.title.y=element_blank())


library(reshape2)
c.allSpecies.tibble.melt=c.allSpecies.tibble %>% dplyr::select(centroidID,n_ORFs,n_RT_domains) %>% reshape2::melt() 
vars=c.allSpecies.tibble.melt %>% pull(variable) %>% unique
totalForCat=vector(mode="list",length(vars))
names(totalForCat)=vars
for(v in vars){
	maxVal=c.allSpecies.tibble.melt[c.allSpecies.tibble.melt$variable==as.character(v),] %>% pull(value) %>% max
	totalForCat[[v]]=maxVal

}
totalForCat=unlist(totalForCat)
c.allSpecies.tibble.melt=c.allSpecies.tibble.melt %>% mutate(value=value*100/totalForCat[variable])

p.heat <- c.allSpecies.tibble.melt %>% ggplot(aes(x=variable, y=centroidID)) + 
    geom_tile(aes(fill=value)) + scale_fill_viridis_c() + 
    theme_minimal() + xlab(NULL) + ylab(NULL)



node.tibble=tibble(node_num=seq(1,tree.centroidTips$Nnode), node=node_num+Ntip(tree.centroidTips), bootstrap=as.numeric(tree.centroidTips$node.label), is.supported=F) %>% mutate(bootstrap=ifelse(bootstrap>80,bootstrap,NA))


supported.nodes=node.tibble %>% filter(bootstrap>80) %>% pull(node)

tree.centroidTips.grouped=groupClade(tree.centroidTips,supported.nodes)

print("saving test...")
save.image("test.RData")

# algunos parches

# este es porque hay un prob con la columna overlappingTx en species sin transcriptoma
c.allSpecies.tibble=c.allSpecies.tibble %>% mutate(overlappingTx=ifelse(n_fragments_overlapping_txs=="no-data",F,overlappingTx)) 
c.allSpecies.tibble=c.allSpecies.tibble %>% mutate(typeOfOverlap=ifelse(overlappingTx==T, "centroid",ifelse(perc_fragments_overlapping_txs > 0 ,"some other fragments","no overlaps")))

tr=ggtree(tree.centroidTips.grouped, aes(colour=specie)) %<+% c.allSpecies.tibble %<+% node.tibble  +
#	geom_tiplab(align=T)+
	geom_tippoint(aes(shape = typeOfOverlap,fill=overlappingTx),size=3) +
	scale_fill_manual(values=c("#e0e0e0","#d53e4f"))+
	scale_shape_manual(values=c(21,4,1)) +
	geom_text(aes(label=bootstrap),nudge_x=-0.15, nudge_y=0.4) +
	theme(legend.position = "right") 


library(aplot)
pdf(paste0(outprefix,".plotTree.pdf"),height=15,width=15)

p.heat %>% insert_left(tr) %>% 
	insert_right(p.length) %>%
	insert_right(p.nFrag)

	#insert_right(p.perc_fragments_overlapping_txs) %>% insert_right(p.nBPs) %>% insert_right(p.perc_bps_overlapping_txs) %>% insert_right(p.nORFs) %>% insert_right(p.nRTs)
#p.nFrag

dev.off()
pdf(paste0(outprefix,".plotTree.onlyTree.pdf"),height=10,width=20)
tr+geom_tiplab(align=F)+geom_text(aes(label=bootstrap),nudge_x=-0.05, nudge_y=0.05) +
dev.off()


save.image(paste0(outprefix,".plotTree.RData"))
