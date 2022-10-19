library(tidyverse)
library(rtracklayer)
library(GenomicFeatures)
library(reshape2)

args=commandArgs(T)
rmgrPath=args[1]
specialFamiliesAusStr=args[2]

tandemRData=as.character(args[5])
outprefix=as.character(args[6])

specialFamiliesAus=strsplit(specialFamiliesAusStr,",")[[1]]

load(rmgrPath)
#load("../METAMATRIX/aus.rmgr")
#specialFamiliesAus=c("rnd-5_family-204","rnd-6_family-117")



RTEBovB_expanded=rm.gr %>% as_tibble %>% mutate(Kimura=as.numeric(Kimura)) %>% filter(is.na(Kimura)==F) %>% filter(TEfam %in% specialFamiliesAus, width>200) %>% mutate(name=paste0("RTEBOVB_",seq(1,NROW(.))))

# exportamos bed para obtener seqs
##rm.gr %>% as_tibble %>% mutate(Kimura=as.numeric(Kimura)) %>% filter(is.na(Kimura)==F) %>% filter(TEfam %in% specialFamiliesAus, width>200) %>% mutate(name=paste0("RTEBOVB_",seq(1,NROW(.)))) %>% select(seqnames,start,end,name,Kimura,strand) %>% write.table(file="RTEBOVB_TREE_FOR_SUMMARY/RTEBovB.bed",quote=F,row.names=F,sep="\t",col.names=F)
#rm.gr %>% as_tibble %>% mutate(Kimura=as.numeric(Kimura)) %>% filter(is.na(Kimura)==F) %>% filter(TEfam=="rnd-5_family-204", width>200) %>% mutate(name=paste0("rnd-5_family-204_",seq(1,NROW(.)))) %>% dplyr::select(seqnames,start,end,name,Kimura,strand) %>% write.table(file="RTEBovB_rnd-5_family-204.bed",quote=F,row.names=F,sep="\t",col.names=F)
#rm.gr %>% as_tibble %>% mutate(Kimura=as.numeric(Kimura)) %>% filter(is.na(Kimura)==F) %>% filter(TEfam=="rnd-6_family-117", width>200) %>% mutate(name=paste0("rnd-6_family-117_",seq(1,NROW(.)))) %>% dplyr::select(seqnames,start,end,name,Kimura,strand) %>% write.table(file="RTEBovB_rnd-6_family-117.bed",quote=F,row.names=F,sep="\t",col.names=F)


#importamos la anotacion de Austrolebi
gff=import(as.character(args[3]))
#gff=import("aus.genes.gff3")
gff=gff[gff$type=="mRNA"]
names(gff)=gff$ID

# filtramos (con los mismos criterios de antes) el GRange de austrolebias

RTEBovB_expanded.gr=rm.gr[is.na(as.numeric(rm.gr$Kimura))==F & rm.gr$TEfam %in% specialFamiliesAus]
RTEBovB_expanded.gr=RTEBovB_expanded.gr[width(RTEBovB_expanded.gr) >200]






# Tree plotting

library(ape)
library(geiger)
library(ggtree)

treePath=as.character(args[4])
rte=read.tree(treePath)

#set.seed(12345)
#rte.sample=rte$tip.label %>% tibble %>% sample_n(100) %>% c





# arbolito circular, destacando grupos
pdf(paste0(outprefix,".tree.pdf"),height=40, width=40)
rte.grouped=groupClade(rte,c(1557,1706))
p.tree=
	ggtree(rte.grouped,aes(color=group),layout='fan')+
	geom_text2(aes(label=label, subset=!isTip), hjust=-.2)+
	#geom_text(aes(label=node), hjust=-.3)+
	scale_color_manual(values=c("#404040", "#2c7bb6", "#ca0020")) #, "#1a9641"))

p.tree.wNode=
	ggtree(rte.grouped,aes(color=group),layout='fan',branch.length="none")+
	geom_text2(aes(label=label, subset=!isTip), hjust=-.2)+
	geom_text(aes(label=node), hjust=-.3)+
	scale_color_manual(values=c("#404040", "#2c7bb6", "#ca0020")) #, "#1a9641"))

p.tree.collapsed=
	collapse(p.tree,node=1706)+
	geom_point2(aes(subset=(node==1706)), shape=23, size=30, fill='#ca0020')+
	theme(legend.position='none')

p.tree.wNode
p.tree.collapsed

# aislando el arbol individual de el grupo3 (aquel que cambio mas rapido)

p.tree.rec=
        ggtree(rte.grouped,aes(color=group),layout='rectangular',branch.length="none")+
        geom_text2(aes(label=label, subset=!isTip), hjust=-.2)+
        geom_text(aes(label=node), hjust=-.3)+
        scale_color_manual(values=c("#404040", "#2c7bb6", "#ca0020")) #, "#1a9641"))

p.group3=viewClade(p.tree.rec, 1707)


dev.off()


# identifica las inserciones en cada grupo
rteGroup3=tips(rte,1706)
#rteGroup3=read.table("group3.lst") %>% tibble %>% pull


#p.tree=ggtree(rte,layout='circular',branch.length='none')


# Heatmaps v2

# generamos GRanges de regiones upstream y downstream de genes y TEs
flank_size=5000

# genes
g.fKbUpstream=gff %>% flank(flank_size,start=T) 
g.fKbDownstream=gff %>% flank(flank_size,start=F) 

#overlaps
# genes
g.upOverlaps=findOverlapPairs(RTEBovB_expanded.gr,g.fKbUpstream)@second
g.downOverlaps=findOverlapPairs(RTEBovB_expanded.gr,g.fKbDownstream)@second
g.intraOverlaps=findOverlapPairs(RTEBovB_expanded.gr,gff)@second
# tes
t.upOverlaps=findOverlapPairs(RTEBovB_expanded.gr,g.fKbUpstream)@first
t.downOverlaps=findOverlapPairs(RTEBovB_expanded.gr,g.fKbDownstream)@first
t.intraOverlaps=findOverlapPairs(RTEBovB_expanded.gr,gff)@first

# agregamos info de la anotacion de genes a la de TEs
t.upOverlaps$GeneID=g.upOverlaps$ID
t.downOverlaps$GeneID=g.downOverlaps$ID
t.intraOverlaps$GeneID=g.intraOverlaps$ID

#matrix

# overlaps in matrix
t.upOverlaps.tibble=t.upOverlaps %>% mcols %>% as_tibble(rownames="ID_GRanges") %>% mutate(upOverlap=T)
t.downOverlaps.tibble=t.downOverlaps %>% mcols %>% as_tibble(rownames="ID_GRanges") %>% mutate(downOverlap=T)
t.intraOverlaps.tibble=t.intraOverlaps %>% mcols %>% as_tibble(rownames="ID_GRanges") %>% mutate(intraOverlap=T)


RTEBovB_expanded$GeneID=NA

RTEBovB_expanded=RTEBovB_expanded %>% mutate(ID_GRanges=paste0(seqnames,"_",start,"_",end,"_",strand)) %>%
	left_join(.,t.upOverlaps.tibble,by="ID_GRanges",suffix=c("",".up")) %>% dplyr::select(-TEfam.up, -TEorder.up, -TEsuperfam.up, -Kimura.up) %>%
	left_join(.,t.downOverlaps.tibble,by="ID_GRanges", suffix=c("",".down")) %>% dplyr::select(-TEfam.down, -TEorder.down, -TEsuperfam.down, -Kimura.down) %>%
	left_join(.,t.intraOverlaps.tibble,by="ID_GRanges", suffix=c("",".intra")) %>% dplyr::select(-TEfam.intra, -TEorder.intra, -TEsuperfam.intra, -Kimura.intra, -GeneID)



# intronless in matrix
g.intronless.tibble=read.table("aus.genes.gff3.TxIntronlessAndLength.tab",header=T) %>% as_tibble

RTEBovB_expanded=RTEBovB_expanded %>%
	left_join(., g.intronless.tibble, by=c("GeneID.up"="Tx")) %>% dplyr::rename(isIntronless.up=isIntronless) %>% dplyr::select(-Gid,-Len) %>%
	left_join(., g.intronless.tibble, by=c("GeneID.down"="Tx")) %>% dplyr::rename(isIntronless.down=isIntronless) %>% dplyr::select(-Gid,-Len) %>%
	left_join(., g.intronless.tibble, by=c("GeneID.intra"="Tx")) %>% dplyr::rename(isIntronless.intra=isIntronless) %>% dplyr::select(-Gid,-Len)


# creamos una columna con un ID que sirve para los arboles
RTEBovB_expanded=RTEBovB_expanded %>% mutate(treeID=paste("RTEBovB",seqnames,start,end,sep="_"))


# filtramos la tabla para generar heatmaps
RTEBovB_expanded.group3=RTEBovB_expanded %>% filter(treeID %in% rteGroup3) 




# load data of tandem groups
#source("/storage2/Thesis/Package/texplorer-reordered/lib/base-complement.R")
loadRData <- function(fileName){
#loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
}


tiles.tibble=loadRData(tandemRData)


#ploting heatmaps
set.seed(12345)
pdf(paste0(outprefix,".heatmaps.pdf"))

# heatmap up,intra,downOverlap(BOOLEAN)
RTEBovB_expanded.group3 %>% dplyr::select(ID_GRanges,TEorder,TEsuperfam,TEfam,upOverlap,intraOverlap,downOverlap,Kimura) %>% mutate(upOverlap=ifelse(is.na(upOverlap),0,1), intraOverlap=ifelse(is.na(intraOverlap),0,1), downOverlap=ifelse(is.na(downOverlap),0,1)) %>%
	dplyr::select(-Kimura) %>%
	melt() %>% ggplot(aes(x=variable,y=ID_GRanges,fill=value)) +geom_tile(width=0.7, height=0.7)+ # x el momento usamos un sample de 200 TEs
	theme_minimal()+
	scale_fill_gradient2(low = "white", high = "black")



nx <- 1
x.breaks <- seq(-nx / 2, 1 + nx/2, length = 5)
x.label <- seq(0, 1, length = 5)
scale_x_continuous(breaks = x.breaks, labels = x.label)

# heatmap Kimura
heatmap.kimura=RTEBovB_expanded %>% filter(treeID %in% rte$tip.label) %>% dplyr::select(treeID,TEorder,TEsuperfam,TEfam,Kimura) %>%
	melt() %>% ggplot(aes(x=variable,y=treeID,fill=value)) +geom_tile()+
	theme_minimal()+
	scale_fill_gradient2(low = "#e0e0e0", mid = "green", high= "red")+theme(axis.text.y = element_blank())+scale_x_discrete(breaks = x.breaks, labels = x.label)
heatmap.kimura.g3=RTEBovB_expanded.group3 %>% dplyr::select(treeID,TEorder,TEsuperfam,TEfam,Kimura) %>%
	melt() %>% ggplot(aes(x=variable,y=treeID,fill=value)) +geom_tile()+
	theme_minimal()+
	scale_fill_gradient2(low = "#e0e0e0", mid = "green", high= "red")+theme(axis.text.y = element_blank())+scale_x_discrete(breaks = x.breaks, labels = x.label)


#dev.off()
# heatmap tandem Score
heatmap.tandem=tiles.tibble %>% left_join(RTEBovB_expanded, by=c("insertion"="ID_GRanges")) %>%  dplyr::select(treeID,Score) %>%
	#mutate(Score=log10(Score))%>% distinct %>% 
	melt %>% tibble %>% ggplot(aes(x=variable,y=treeID,fill=value)) +geom_tile() + theme_minimal()+scale_fill_gradient2(low = "#e0e0e0",midpoint=0.5,high = "darkblue",mid="blue")+theme(axis.text.y = element_blank())+scale_x_discrete(breaks = x.breaks, labels = x.label)


# estabamos haciendo una sola matrix, para un heatmap con diferentes paletas x columna
RTEBovB_expanded %>% left_join(tiles.tibble,by=c("ID_GRanges"="insertion")) %>% filter(treeID %in% rte$tip.label) %>% select(treeID, Kimura,Score,isIntronless.up,isIntronless.intra,isIntronless.down) %>% melt() %>% head # <- no funca el melt



pdf(paste0(outprefix,".panel.pdf"),height=30,width=30)

#heatmap.kimura.g3 %>% aplot::insert_left(p.group3)
heatmap.kimura.g3 %>% aplot::insert_left(p.group3) %>% aplot::insert_right(heatmap.tandem)
#heatmap.tandem %>% aplot::insert_left(p.tree.rec)
heatmap.kimura %>% aplot::insert_left(p.tree.rec)
heatmap.kimura %>% aplot::insert_left(p.tree.rec) %>% aplot::insert_right(heatmap.tandem)

#p.tree.rec %>% aplot::insert::righ
dev.off()

rm(rm.gr)
gc()
save.image(paste0(outprefix,".FinalSummaryFigure.RData"))
