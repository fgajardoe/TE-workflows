
library(tidyverse)
args=commandArgs(T)

g_legend<-function(a.gplot){
        a.gplot=a.gplot+theme(legend.position="right")
        tmp <- ggplot_gtable(ggplot_build(a.gplot))
        leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
        legend <- tmp$grobs[[leg]]
        return(legend)
}

outprefix=as.character(args[1])


min_num_bps=19012566 # mean ; primer quantile 368995



library(ape)
library(ggtree)
tree=read.tree(args[5])
p.tree=ggtree(tree) + geom_tiplab(align=TRUE) + hexpand(.01)

pairedFilesArgs=args[6:length(args)]




gffFiles=pairedFilesArgs[1:(length(pairedFilesArgs)/3)]
sizeFiles=pairedFilesArgs[(1+(length(pairedFilesArgs)/3)):(2*length(pairedFilesArgs)/3)]
tableFiles=pairedFilesArgs[(1+((length(pairedFilesArgs)/3)*2)):length(pairedFilesArgs)]



# cargamos los datos de cada archivo de tamaños genómicos
gz=c()
i=1
for(f in sizeFiles){
	print(paste0("working on ",f))
	sp=strsplit(f,split="\\.")[[1]][1] %>% as.character()
	gz[[sp]]=read.table(as.character(f), fill=T, quote="", col.names=c("seqnames","width")) %>% tibble %>% pull(width) %>% sum
	print(paste0("Genome size of ",sp,": ",as.character(gz[[sp]])))
}

# diccionario con genome sizes
gz=unlist(gz)

# cargamos los RData que contienen las tablas con los datos
tf.lst=vector(mode="list",length(tableFiles))
names(tf.lst)=tableFiles

for(tf in tableFiles){
	print(paste0("loading ",tf,"..."))
	load(tf)
	tf.lst[[tf]]=TEsuperfam.nBPs

}



# cargamos las anotaciones de genes/CDS/exones/etc...
library(rtracklayer)
speciesOrder=tableFiles %>% gsub(".RData","",.)

gff.lst=vector(mode="list",length(gffFiles))
names(gff.lst)=speciesOrder
i=1

plotCDS=F
if(any(plotCDS)){
	for(ann in gffFiles){
		print(paste0("loading genome annotation ",ann))
		gff.gr=import(ann)
		sp=speciesOrder[i]	

		# calculamos el numero de bps en cdss
		gff.lst[[sp]]=gff.gr %>% as.data.frame %>% as.tibble %>% filter(type=="CDS") %>% pull(width) %>% sum
		i=i+1
	}

	# diccionario con el numero de bps en CDSs por especie.
	gff.lst=unlist(gff.lst)


	# tibble con las bps por CDSs
	d.coding=tibble(TEorder="Coding regions",specie=names(gff.lst),`Total bps`=gff.lst)
}
mergedTable=do.call("rbind",tf.lst)

mergedTable=mergedTable %>% mutate(perc_of_genome_size=nBPs*100/gz[specie])
mergedTable$TEsuperfam=as.character(mergedTable$TEsuperfam)

save(mergedTable,file=paste0(outprefix,".merged.RData"))

# filtramos la tabla luego de guardarla, para que solo tenga efecto en la parte de visualización y no en los datos.

# params:
collapse_lower_than_mean_contrib=as.logical(args[2]) #F
onlyCoreTEsuperfams=as.logical(args[3]) #T

#options(scipen=100000000000)
elementsCollapsedToOther=strsplit(as.character(args[4]),",")
elementsCollapsedToOther=elementsCollapsedToOther[[1]]


numberOfSpecies=mergedTable$specie %>% unique %>% length
coreTEsupefams=mergedTable %>% select(TEsuperfam,specie) %>% group_by(TEsuperfam) %>% summarise(nSpecies=n_distinct(specie)) %>% filter(nSpecies==numberOfSpecies) %>% pull(TEsuperfam)

if(any(collapse_lower_than_mean_contrib)){
	mergedTable=mergedTable %>% mutate(TEsuperfam=ifelse(TEsuperfam %in% coreTEsupefams,ifelse(nBPs<min_num_bps,"Lower than mean",TEsuperfam),TEsuperfam))
}

`%notin%` <<- Negate(`%in%`)
if(any(onlyCoreTEsuperfams)){
	mergedTable=mergedTable %>% mutate(TEsuperfam.complete=TEsuperfam) %>% mutate(TEsuperfam=ifelse(TEsuperfam %notin% coreTEsupefams, "No core",TEsuperfam))
}
mergedTable$TEsuperfam=factor(mergedTable$TEsuperfam)


# Plotting TE superfams
maxVal=mergedTable %>% filter(nBPs==max(.$nBPs)) %>% select(nBPs) %>% unlist
x.lim=c(0,maxVal/1000000)
orders=mergedTable %>% pull(TEorder) %>% unique

p.lst=vector(mode="list",length(orders))
names(p.lst)=orders
p2.lst=vector(mode="list",length(orders))
names(p2.lst)=orders




# barplots TE superfams separado por ordenes en graficos distintos (paleta chicago)

library(ggsci)
for(o in orders){
	p.lst[[o]]=mergedTable %>% filter(TEorder==o) %>% ggplot(aes(x=nBPs/1000000,y=TEsuperfam,fill=specie)) + geom_bar(stat="identity",position="dodge",width=0.8)+theme_classic()+theme(legend.position="none")+xlab("Mbps")+xlim(x.lim) + ggtitle(as.character(o)) +scale_fill_uchicago()
	p2.lst[[o]]=mergedTable %>% filter(TEorder==o) %>% ggplot(aes(x=perc_of_genome_size,y=TEsuperfam,fill=specie)) + geom_bar(stat="identity",position="dodge",width=0.8)+theme_classic()+theme(legend.position="none")+xlab("% of genome size")+xlim(0,25) + ggtitle(as.character(o)) +scale_fill_uchicago()
}
p.legend=p2.lst[[o]] + theme(legend.position="bottom")
p.legend=g_legend(p.legend)
p2.lst[["legend"]]=p.legend
p.lst[["legend"]]=p.legend

library(grid)
library(gridExtra)

plot.ncol=2
title=paste0("Genome size contribution of TE superfamilies (bps)")
title2=paste0("Genome size contribution of TE superfamilies (%)")

title.grob=textGrob(as.character(title),gp=gpar(fontsize=20,font=3))
title2.grob=textGrob(as.character(title2),gp=gpar(fontsize=20,font=3))
grid=do.call("arrangeGrob",c(p.lst, ncol=plot.ncol))
grid2=do.call("arrangeGrob",c(p2.lst, ncol=plot.ncol))

pdf(paste0(outprefix,".pdf"), height=26 # solia ser 13 para LINE,DNA,LTR,SINE
    , width=13)
grid.arrange(grid,top=title.grob)
grid.arrange(grid2,top=title2.grob)
dev.off()


## Genome size contrib TE orders and shared
#speciesOrder=c("Austrolebias_charrua","Nematolebias_whitei","Cynopoecilus_melanotaenia","Austrofundulus_limnaeus","Kryptolebias_marmoratus","Nothobranchius_furzeri","Oryzias_latipes")

mergedTable=mergedTable %>% mutate(TEorder=as.character(TEorder))
if(length(elementsCollapsedToOther)>0){
	mergedTable=mergedTable %>% mutate(TEorder=as.character(TEorder))
	mergedTable=mergedTable %>% mutate(TEorder=ifelse(TEorder %in% elementsCollapsedToOther, "Other TEs",TEorder))

}

# Remove insertions with no order classification
mergedTable=mergedTable %>% filter(is.na(TEorder)==F)



# data for genome size contrib plot (TEorders)
print("Calculating number of bps for TE orders")
d=mergedTable %>% group_by(TEorder,specie) %>% summarise(`Total bps`=sum(nBPs))

if(any(plotCDS)){
	d=rbind(d,d.coding)
}


d.order=d$TEorder %>% unique
d.order=d.order[c(
			    6, #coding
			    1, #dna core
			    2, #line
			    5, #sine
			    3, #ltr
			    4)] #other
d.order=c(d.order,"Other genomic elements")
#d.order=factor(d.order,levels=d.order, ordered=T)
d.order.pal=c(
		   "#ff6f00ff", #coding
		   "#8a4198ff", #dna core
		   "#c71000ff", #line
		   "#008ea0ff", #sine
		   "#5a9599ff", #ltr
		   "#d1cfa1ff", #other TEs
		   "#e2e2e2fb") #other genomic elements

if(!any(plotCDS)){
	d.order.pal=d.order.pal[c(2:length(d.order.pal))]
}


d$specie=factor(d$specie,levels=rev(speciesOrder),ordered=T)

# data for genome size contrib plot (TEorders and conservation)
print("Calculating number of bps for TE orders highlighting core and unique TE families")
d.cons=mergedTable %>% mutate(Conservation=ifelse(TEsuperfam %in% coreTEsupefams,"Core","No core"))  %>% group_by(TEorder,specie,Conservation) %>% summarise(`Total bps`=sum(nBPs))

d.cons$specie=factor(d.cons$specie,levels=rev(speciesOrder),ordered=T)


d.cons.order=d.cons$TEorder %>% unique
d.cons.order=d.cons.order[c(
			    #10, #coding
			    1, #dna core
			    2, #line
			    5, #sine
			    3, #ltr
			    4)] #other


d.cons.order.pal=c(
		   #"#ff6f00ff", #coding
		   "#8a4198ff", #dna core
		   "#c71000ff", #line
		   "#008ea0ff", #sine
		   "#5a9599ff", #ltr
		   "#d1cfa1ff") #other

d.cons$TEorder=factor(d.cons$TEorder,ordered=T,levels=as.character(d.cons.order))


# data for genome size contrib plot (TEsuperfams core)
print("Calculating number of bps for core TE superfams")

d.tesuperfam=mergedTable %>% mutate(TEsuperfam=as.character(TEsuperfam.complete)) %>% mutate(Conservation=ifelse(TEsuperfam %in% coreTEsupefams,"Core","No core"))  %>% group_by(TEsuperfam,specie,Conservation) %>% summarise(`Total bps`=sum(nBPs))
nSpecies.v=d.tesuperfam %>% distinct(TEsuperfam,specie) %>% group_by(TEsuperfam) %>% summarise(.,nSpecies=length(TEsuperfam)) %>% pull(nSpecies, TEsuperfam)

# agregamos info para determinar si son unicos (unique)
d.tesuperfam=d.tesuperfam %>% mutate(nSpecies=nSpecies.v[TEsuperfam]) 
d.tesuperfam=d.tesuperfam %>% mutate(Conservation=ifelse(nSpecies==1,"Unique",Conservation))

d.tesuperfam$specie=factor(d.tesuperfam$specie,levels=rev(speciesOrder),ordered=T)


d.tesuperfam.order=d.tesuperfam$TEsuperfam %>% unique
#d.tesuperfam.order=d.tesuperfam.order[c(
#			    #10, #coding
#			    1, #dna core
#			    2, #line
#			    5, #sine
#			    3, #ltr
#			    4)] #other



d.tesuperfam$TEsuperfam=factor(d.tesuperfam$TEsuperfam,ordered=T,levels=as.character(d.tesuperfam.order))


# plotting

library(ggsci)
library(RColorBrewer)
pdf(paste0(outprefix,".genomeSizePerSpecie.pdf"),width=12,height=5)

# agregamos entradas para llegar al al genome size

for(sp in speciesOrder){
	print(paste0("calculating substraction from genome size for ",sp,"..."))
	sumSp=d %>% filter(specie==sp) %>% pull(`Total bps`) %>% sum
	d.other=tibble(TEorder="Other genomic elements",specie=sp,`Total bps`=gz[[sp]]-sumSp)
	d=rbind(d,d.other)

}
d$TEorder=factor(d$TEorder,levels=d.order,ordered=T)

# barplot species tree
p.barplot=d %>% ggplot(aes(x=`Total bps`/1000000000,y=specie, fill=TEorder))+geom_bar(stat="identity",position=position_stack(reverse = TRUE),width=0.5) +theme_classic()+xlab("Gbps")+scale_fill_manual(values=c(d.order.pal))+theme(axis.title.y=element_blank())
p.barplot.perc=d %>% ggplot(aes(x=`Total bps`/1000000000,y=specie, fill=TEorder))+geom_bar(stat="identity",position=position_fill(reverse = T),width=0.5) +theme_classic()+xlab("Proportion of genome size")+scale_fill_manual(values=c(d.order.pal))+theme(axis.title.y=element_blank())



TipsToKeep=speciesOrder
TipsToRemove=tree$tip.label %>% .[. %notin% TipsToKeep]
#TipsToRemove=tree$tip.label %>% .[. %notin% c("Austrolebias_charrua","Cynopoecilus_melanotaenia","Nematolebias_whitei","Austrofundulus_limnaeus","Kryptolebias_marmoratus","Nothobranchius_furzeri","Oryzias_latipes") ]
tree.panel=drop.tip(tree,TipsToRemove)
#tree.panel=root(tree.panel,outgroup="Anableps_anableps")
tree.panel=ape::rotate(phy=tree.panel,c(1,2))

p.tree.panel=ggtree(tree.panel) + hexpand(.01)
#p.tree.panel=viewClade(p.tree, MRCA(p.tree, speciesOrder))

library(aplot)
p.barplot %>% insert_left(p.tree.panel) %>% insert_right(p.barplot.perc)



# barplot core and no core
d.cons %>% ggplot(aes(x=`Total bps`/1000000000,y=specie, fill=TEorder,group=Conservation))+geom_bar(stat="identity",position=position_stack(reverse = TRUE),width=0.5)+facet_grid(~ Conservation) +theme_classic()+xlab("Gbps")+scale_fill_manual(values=c(d.cons.order.pal))




# All core, no core and unique filled barplot

n_superfams=d.tesuperfam %>% pull(TEsuperfam) %>% unique %>% length
n_superfams.core=d.tesuperfam %>% filter(Conservation=="Core") %>% pull(TEsuperfam) %>% unique %>% length
pal.tesuperfams=colorRampPalette(brewer.pal(9,name="Set1"))(n_superfams)
pal.tesuperfams2=colorRampPalette(brewer.pal(9,name="Set1"))(n_superfams.core)
d.tesuperfam %>% ggplot(aes(x=`Total bps`/1000000000,y=specie, fill=TEsuperfam,group=Conservation))+geom_bar(stat="identity",position="fill")+
													     facet_grid(~ Conservation)+
													     theme_classic()+
													     theme(legend.position="bottom")+
													     xlab("Proportion of base pairs")+
													     scale_fill_manual(values = pal.tesuperfams)
# Core barplot filled

# dict TEorder - TEsuperfam												
TEorder.v=mergedTable %>% pull(TEorder,TEsuperfam.complete)

# generate TEorder specific palette
d.tesuperfam.core=d.tesuperfam %>% ungroup %>% filter(Conservation=="Core") %>% mutate(TEsuperfam=as.character(TEsuperfam)) %>% mutate(TEorder=TEorder.v[TEsuperfam])

# remove unknown TE superfamilies
hideUnknown=T
if(any(hideUnknown)){
	d.tesuperfam.core=d.tesuperfam.core %>% filter(TEsuperfam!="Unknown")
}


nTEsuperfamsPerTEorder.v=d.tesuperfam.core %>% ungroup %>% select(TEsuperfam,TEorder) %>% distinct(TEsuperfam,TEorder) %>% group_by(TEorder) %>% summarise(nTEsuperfams=n()) %>% pull(nTEsuperfams,TEorder)



#pal.TEorders=pal_futurama(palette="planetexpress")(nTEsuperfamsPerTEorder.v %>% length)
#pal.TEorders=pal.TEorders[c()]
pal.TEorders=c(
		   DNA="#8a4198ff", #dna core
		   LINE="#c71000ff", #line
		   SINE="#008ea0ff", #sine
		   LTR="#5a9599ff", #ltr
		   `Other TEs`="#d1cfa1ff") #other
#names(pal.TEorders)=names(nTEsuperfamsPerTEorder.v)


pal.TEorderSpecific=c()
for(TEo in (d.tesuperfam.core$TEorder %>% unique)){
	d.tesuperfam.core

	TEorderColor=pal.TEorders[TEo]
	n_superfams_in_TEorder=as.numeric(nTEsuperfamsPerTEorder.v[TEo])
	pal.TEorderSpecific.TEo=colorRampPalette(c(TEorderColor,"black"))(3)
	pal.TEorderSpecific.TEo=colorRampPalette(c(pal.TEorderSpecific.TEo[1],pal.TEorderSpecific.TEo[2]))(n_superfams_in_TEorder)
	TEsuperfamsInTEorder=d.tesuperfam.core %>% ungroup %>% select(TEsuperfam,TEorder) %>% distinct %>% filter(TEorder==TEo) %>% pull(TEsuperfam)
	names(pal.TEorderSpecific.TEo)=TEsuperfamsInTEorder
	pal.TEorderSpecific=c(pal.TEorderSpecific,pal.TEorderSpecific.TEo)

}
library(scales)
#pdf("test.pdf")
#show_col(pal.TEorderSpecific)
#dev.off()

# plot all core TE superfams with all species in the same plot
d.tesuperfam.core %>% ggplot(aes(x=`Total bps`/1000000000,y=specie, fill=TEsuperfam,group=TEorder))+geom_bar(stat="identity",position="fill")+
													     theme_classic()+
													     theme(legend.position="bottom")+
													     xlab("Proportion")+
													     ggtitle("Core TE superfamilies")+
													     scale_fill_manual(values = pal.TEorderSpecific)

dev.off()



library(ggrepel)

pdf(paste0(outprefix,".genomeSizePerSpecie.TEsuperfam.pdf"),width=12,height=4)
flag=1
scaleSize=50000000
for(sp in speciesOrder){
	print(paste0("plotting TE superfams barplot with TE order colours for ",sp))
	p.sp=d.tesuperfam.core %>% filter(specie==sp) %>% ggplot(aes(x=`Total bps`,y=specie, fill=TEsuperfam,group=TEorder))+
		geom_bar(stat="identity",position="stack",width=1)+
		theme_void()+
		#theme_classic()+
		theme(legend.position="bottom", axis.title.y=element_blank(),axis.text.y=element_blank(), plot.title = element_text(hjust = 0.05,size=18,face="italic"))+
		xlab("Gbps")+
		ggtitle(gsub(pattern="_",replacement=" ",x=sp))+
#		labs(caption=sp)+
		scale_fill_manual(values = pal.TEorderSpecific)+
		geom_text(aes(y="TEsuperfamily",label=TEsuperfam,
			      size=`Total bps`,
			      alpha=`Total bps`
			      ),angle=45,hjust = 0,
			  position=position_stack(vjust=0.5))+
		geom_segment(aes(x=0,y="scale",xend=scaleSize,yend="scale"))+annotate("text",label=paste0("\n",scaleSize/1000000," Mbps"),x=scaleSize/2,y="scale")

	if(flag==1){												   #0.5
		p.sp.legend=g_legend(p.sp)
		flag=0
	}
	p.sp=p.sp + theme(legend.position="none")
	print(p.sp)
}
grid.arrange(p.sp.legend)
dev.off()














save.image(paste0(outprefix,".RData"))
