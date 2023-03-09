
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
tree=read.tree(args[8])
refSpecie=as.character(args[5])
colorFile=as.character(args[6])
ylimit=as.character(args[7])
p.tree=ggtree(tree) + geom_tiplab(align=TRUE) + hexpand(.01)

pairedFilesArgs=args[9:length(args)]




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
#speciesOrder=speciesOrder[c(1,2,3,4,6,5)]

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
if(args[4] != "none"){
elementsCollapsedToOther=strsplit(as.character(args[4]),",")
elementsCollapsedToOther=elementsCollapsedToOther[[1]]
} else {
elementsCollapsedToOther="none"
}

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
#if(elementsCollapsedToOther!="none"){
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
d.order=c(d.order,"No-repetitive")
#d.order=c(d.order,"Other genomic elements")
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
	d.other=tibble(TEorder="No-repetitive",specie=sp,`Total bps`=gz[[sp]]-sumSp)
	d=rbind(d,d.other)

}
d$TEorder=factor(d$TEorder,levels=d.order,ordered=T)

# barplot species tree
barplots.xlim=max(gz)/1000000000
barplots.xlim=barplots.xlim %>% round(1)
p.barplot=d %>% ggplot(aes(x=`Total bps`/1000000000,y=specie, fill=TEorder))+geom_bar(stat="identity",position=position_stack(reverse = TRUE),width=0.5) +theme_classic()+xlab("Gbps")+scale_fill_manual(values=c(d.order.pal))+theme(axis.title.y=element_blank())+xlim(values=c(0,barplots.xlim))
p.barplot.perc=d %>% ggplot(aes(x=`Total bps`/1000000000,y=specie, fill=TEorder))+geom_bar(stat="identity",position=position_fill(reverse = T),width=0.5) +theme_classic()+xlab("Proportion of genome size")+scale_fill_manual(values=c(d.order.pal))+theme(axis.title.y=element_blank())



TipsToKeep=speciesOrder
TipsToRemove=tree$tip.label %>% .[. %notin% TipsToKeep]
#TipsToRemove=tree$tip.label %>% .[. %notin% c("Austrolebias_charrua","Cynopoecilus_melanotaenia","Nematolebias_whitei","Austrofundulus_limnaeus","Kryptolebias_marmoratus","Nothobranchius_furzeri","Oryzias_latipes") ]
tree.panel=drop.tip(tree,TipsToRemove)
#tree.panel=root(tree.panel,outgroup="Anableps_anableps")
tree.panel=ape::rotate(phy=tree.panel,c(1,2))

p.tree.panel=ggtree(tree.panel) + hexpand(.01)
#p.tree.panel=viewClade(p.tree, MRCA(p.tree, speciesOrder))

# barplot core and no core

# preparing text data
d.cons.text=d.cons %>% ungroup() %>% dplyr::select(specie,Conservation,`Total bps`) %>% group_by(specie,Conservation)%>% summarise(`Total Mbps`=sum(`Total bps`)/1000000) %>% mutate(label=paste0(format(`Total Mbps`,digit=2)," Mbps (",Conservation,")")) %>% mutate(TEorder="Other TEs")

# ploting
p.coreNoCore=d.cons %>% ggplot(aes(x=`Total bps`/1000000000,y=specie, fill=TEorder,group=Conservation))+geom_bar(stat="identity",position="dodge",width=0.5)+
		xlim(c(0,barplots.xlim))+
	#position_stack(reverse = TRUE),width=0.5)+facet_grid(~ Conservation) +
		xlab("Gbps")+
		theme_void()+
#		scale_fill_manual(values=c("#4d4d4d","#878787"))+
		scale_fill_manual(values=c(d.cons.order.pal))+
		geom_text(data=d.cons.text,aes(label=label,x=0.7*barplots.xlim,group=`Conservation`,colour=`TEorder`),size=3, position=position_dodge(width = .9))+
		scale_colour_manual(values=c("#4d4d4d","#878787"))+
		theme(axis.text.y=element_blank(),legend.position="none",axis.title.y=element_blank())

# composite panel
library(aplot)
p.barplot %>% insert_left(p.tree.panel) %>% insert_right(p.coreNoCore)
#p.barplot %>% insert_left(p.tree.panel) %>% insert_right(p.barplot.perc)







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

dev.off() # cierra .genomeSizePerSpecie.pdf



library(ggrepel)

pdf(paste0(outprefix,".genomeSizePerSpecie.TEsuperfam.pdf"),width=18, height=4*length(speciesOrder))
#,height=4)
flag=1
scaleSize=50000000

p.TEsuperfamBarplot.lst=vector(mode="list",length(speciesOrder))
names(p.TEsuperfamBarplot.lst)=speciesOrder


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
	#print(p.sp)
	p.TEsuperfamBarplot.lst[[sp]]=p.sp
}
grid.arrange(p.sp.legend)

do.call("grid.arrange", c(p.TEsuperfamBarplot.lst, ncol=1))



dev.off() # cierra .genomeSizePerSpecie.TEsuperfam.pdf




# log2fold-change calculation



library(ggsci)
mergedTable$specie=factor(mergedTable$specie,levels=speciesOrder,ordered=T)


if(colorFile=="none"){
	pal.species=rep("black",length(speciesOrder))
}
if(colorFile=="default"){
	pal.species=pal_uchicago()(length(speciesOrder))
} 
if(colorFile!="none" && colorFile!="default"){

	pal.species=read.table(colorFile, sep="\t",comment.char="@",col.names=c("Specie","Colour")) %>% tibble %>% pull(Colour, Specie)

	#pal.species=c(Austrolebias_charrua="#767676FF",
	#	      Nematolebias_whitei="#155F83FF",
	#	      Cynopoecilus_melanotaenia="#FFA319FF",
	#	      Austrofundulus_limnaeus="#800000FF",
	#	      Kryptolebias_marmoratus="#8A9045FF",
	#	      Nothobranchius_furzeri="#C16622FF",
	#	      Oryzias_latipes="#8F3931FF")
}




coreTEsupefams.noUnk=coreTEsupefams[coreTEsupefams!="Unknown"]
coreTEsupefams.noUnk=c(coreTEsupefams.noUnk,"No-repetitive")

TEorders=TEorder.v %>% unique

p.logfc.lst=vector(mode="list",length(c(coreTEsupefams.noUnk,TEorders)))
names(p.logfc.lst)=c(coreTEsupefams.noUnk,TEorders)
#nBPs

mergedTable.NonRepetitive=tibble(TEsuperfam="No-repetitive",
				 #nBPs= (d %>% filter(TEorder=="Other genomic elements") %>% pull(`Total bps`)),
				 nBPs= (d %>% filter(TEorder=="No-repetitive") %>% pull(`Total bps`)),
				 TEorder="No-repetitive",
				 #TEorder="None",
				 specie=(d %>% filter(TEorder=="No-repetitive")%>% pull(specie)),
				 #specie=(d %>% filter(TEorder=="Other genomic elements")%>% pull(specie)),
				 perc_of_genome_size=(nBPs/gz[specie])*100,
				 TEsuperfam.complete="No-repetitive")

mergedTable=rbind(mergedTable,mergedTable.NonRepetitive)

# define ymax and min
#ymaxNumber=12
ymaxNumber=0


# Log2FC analysis for TEorders
mergedTable.TEorders=mergedTable %>% select(-TEsuperfam) %>% group_by(TEorder,specie) %>% summarise(sum(nBPs)) %>% mutate(nBPs=`sum(nBPs)`) %>% select(-`sum(nBPs)`)
print("ploting Log2FC for TEorders...")
for(TEo in TEorders){
	print(TEo)
	nBPsRef=mergedTable.TEorders %>% filter(TEorder==TEo, specie==refSpecie) %>% pull(nBPs) %>% unlist

	p.logfc.lst[[TEo]]=mergedTable.TEorders %>% filter(TEorder==TEo, specie!=refSpecie) %>% mutate(log2FC=log2(nBPsRef/nBPs)) %>% ggplot(aes(x=specie,y=log2FC,fill=specie))+geom_bar(stat="identity",width=0.7)+theme_minimal()+ggtitle(TEo,paste0(refSpecie,": ",round(nBPsRef/1000000)," Mbps"))+ theme(title=element_text(size=14),axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=10),axis.text.y=element_text(size=12), legend.position="none")+scale_fill_manual(values=pal.species)+geom_hline(yintercept=c(0), linetype="dotted",colour="#767676FF")
	#+annotation_logticks(base=2,sides="l")


	#+scale_y_continuous(trans = log2_trans(),
        #             breaks = trans_breaks("log2", function(x) 2^x),
         #            labels = trans_format("log2", math_format(2^.x))) 

	ymaxNumberFam=mergedTable.TEorders %>% filter(specie!=refSpecie) %>% mutate(log2FC=log2(nBPsRef/nBPs)) %>% pull(log2FC) %>% abs %>% max %>% ceiling
	if(ymaxNumberFam>ymaxNumber){
		ymaxNumber=ymaxNumberFam
	}
}


#ymaxNumber=0
for(TEsf in coreTEsupefams.noUnk){
	nBPsRef=mergedTable %>% filter(TEsuperfam==TEsf, specie==refSpecie) %>% pull(nBPs) %>% unlist
	p.logfc.lst[[TEsf]]=mergedTable %>% filter(TEsuperfam==TEsf, specie!=refSpecie) %>% mutate(log2FC=log2(nBPsRef/nBPs)) %>% ggplot(aes(x=specie,y=log2FC,fill=specie))+geom_bar(stat="identity",width=0.7)+theme_minimal()+ggtitle(TEsf,paste0(refSpecie,": ",round(nBPsRef/1000000)," Mbps"))+ theme(title=element_text(size=14),axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=10),axis.text.y=element_text(size=12), legend.position="none")+scale_fill_manual(values=pal.species)+geom_hline(yintercept=c(0), linetype="dotted",colour="#767676FF") 

		#ggplot(aes(x=specie,y=log2FC,fill=specie))+geom_bar(stat="identity",width=0.7)+theme_minimal()+ggtitle(TEsf)+ theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1), legend.position="none")+scale_fill_manual(values=pal.species)+geom_hline(yintercept=c(0), linetype="dotted",colour="#767676FF")
	ymaxNumberFam=mergedTable %>% filter(specie!=refSpecie) %>% mutate(log2FC=log2(nBPsRef/nBPs)) %>% pull(log2FC) %>% abs %>% max %>% ceiling
	if(ymaxNumberFam>ymaxNumber){
		ymaxNumber=ymaxNumberFam
	}
}



if(ylimit=="default"){
	forceThisYlimit=NULL
} else {
	forceThisYlimit=as.numeric(ylimit)
}
if(!is.null(forceThisYlimit)){
	ymaxNumber=as.numeric(forceThisYlimit)
}

print(paste0("ylim is: +/-", ymaxNumber))

#for(TEsf in coreTEsupefams.noUnk){
if(ylimit!="disable"){
	for(TEsf in names(p.logfc.lst)){
		p.logfc.lst[[TEsf]]=p.logfc.lst[[TEsf]]+ylim(c(-ymaxNumber,ymaxNumber))
	}
}
#
p.logfc.perc.lst=vector(mode="list",length(coreTEsupefams.noUnk))
names(p.logfc.perc.lst)=coreTEsupefams.noUnk
#%

# define ymax and min
#ymaxNumber=mergedTable %>% filter(specie!=refSpecie) %>% mutate(log2FC=log2(perc_of_genome_sizeRef/perc_of_genome_size)) %>% pull(log2FC) %>% abs %>% max %>% ceiling


for(TEsf in coreTEsupefams.noUnk){
	perc_of_genome_sizeRef=mergedTable %>% filter(TEsuperfam==TEsf, specie==refSpecie) %>% pull(perc_of_genome_size) %>% unlist
	p.logfc.perc.lst[[TEsf]]=mergedTable %>% filter(TEsuperfam==TEsf, specie!=refSpecie) %>% mutate(log2FC=log2(perc_of_genome_sizeRef/perc_of_genome_size)) %>% ggplot(aes(x=specie,y=log2FC,fill=specie))+geom_bar(stat="identity",width=0.7)+theme_minimal()+ggtitle(TEsf)+ theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1), legend.position="none")+scale_fill_manual(values=pal.species)+ylim(c(-ymaxNumber,ymaxNumber))+geom_hline(yintercept=c(0), linetype="dotted",colour="#767676FF")

}



# preparing PDF page size
nplots=length(p.logfc.lst) %>% as.numeric()
maxNumberOfPlots=1000
multiples_of_5=(1:maxNumberOfPlots)%%5
names(multiples_of_5)=seq(1:maxNumberOfPlots)
multiples_of_5=multiples_of_5 %>% .[.==0]
page_height.v=rep(0,maxNumberOfPlots)
names(page_height.v)=seq(1,maxNumberOfPlots)
z=1
for(i in names(multiples_of_5)){ for(j in seq(1:5)){ page_height.v[z]=i; z=z+1 }  }
page_height=page_height.v[nplots] %>% as.numeric()

p.logfc.legend=p.logfc.lst[[1]]+theme(legend.position="right")
p.logfc.legend=g_legend(p.logfc.legend)
p.logfc.lst[["Lengend"]]=p.logfc.legend
p.logfc.perc.lst[["Lengend"]]=p.logfc.legend
grid=do.call("grid.arrange", c(p.logfc.lst, ncol=5))
grid.perc=do.call("grid.arrange", c(p.logfc.perc.lst, ncol=5))

save(p.logfc.lst,file=paste0(outprefix,".logfc.RData"))

title.grob=textGrob(as.character(paste0("Log2 fold change of the number of bases in ",refSpecie," vs all other species")),gp=gpar(fontsize=20,font=3))
title.grob.perc=textGrob(as.character(paste0("Log2 fold change of the proportion of bases in ",refSpecie," vs all other species")),gp=gpar(fontsize=20,font=3))

# plotting PDF
pdf(paste(outprefix,".log2FC.pdf",sep=""),height=page_height,width=20)
grid.arrange(grid,top=title.grob)
grid.arrange(grid.perc,top=title.grob.perc)
dev.off()

print("I'll start calculate contribution of TEs to delta genome size...")
tableTibble=tibble(`title`=NA,`value`=NA,`specie`=NA,`refSpecie`=NA)
#tableTibble=tibble(`title`="Reference specie",`value`=NA,`specie`=refSpecie,`rSpecie`=NA)

print("tibble created. now I'll screen each specie with the reference")

for(refSpecieForAllvsAll in speciesOrder){
	refSpecieForAllvsAllGenomeSize=as.numeric(gz[[refSpecieForAllvsAll]])
	print("genome size obtained for reference")

	# calculating the contribution of each TE order to the differences in the genome size between the reference specie and each other in the panel
	for(sp in speciesOrder){
		print(paste0("==> ",sp," and ",refSpecieForAllvsAll," (ref)"))
		if(sp != refSpecieForAllvsAll){

			spGenomeSize=as.numeric(gz[[sp]])
			print(paste0("genome size obtained for ",sp))

			deltaGS=abs(refSpecieForAllvsAllGenomeSize-spGenomeSize)

			print("delta genome size calculated")



			totalDeltaTEo=0
			tableTibble=rbind(tableTibble,c("Delta genome size", deltaGS, sp,refSpecieForAllvsAll))
			for(TEo in TEorders){

				nBPsRef=mergedTable.TEorders %>% filter(TEorder==TEo, specie==refSpecieForAllvsAll) %>% pull(nBPs) %>% unlist
				nBPsSp=mergedTable.TEorders %>% filter(TEorder==TEo, specie==sp) %>% pull(nBPs) %>% unlist
				deltaTEo=nBPsRef-nBPsSp
				percDeltaGSdueTEo=deltaTEo*100/deltaGS
				totalDeltaTEo=totalDeltaTEo+deltaTEo
				# report
				print(paste0("==> For ",refSpecieForAllvsAll," vs ",sp," <=="))
				print(paste0("Delta genome size = ",deltaGS," bps"))
				print(paste0("Delta TE order(",TEo,") = ",deltaTEo))
				print(paste0("Weight of Delta TE order(",TEo,") on Delta genome size",percDeltaGSdueTEo))

				# table
				tableTibble=rbind(tableTibble,c(paste0("Delta TE order(",TEo,")"),deltaTEo,sp,refSpecieForAllvsAll))
				tableTibble=rbind(tableTibble,c(paste0("Weight of Delta TE order(",TEo,") on Delta genome size"),percDeltaGSdueTEo,sp,refSpecieForAllvsAll))

			}
			print(paste0("Sum of TEorders differences: ",totalDeltaTEo))
			for(TEo in TEorders){
				nBPsRef=mergedTable.TEorders %>% filter(TEorder==TEo, specie==refSpecieForAllvsAll) %>% pull(nBPs) %>% unlist
				nBPsSp=mergedTable.TEorders %>% filter(TEorder==TEo, specie==sp) %>% pull(nBPs) %>% unlist
				deltaTEo=nBPsRef-nBPsSp
				percSumDeltaTEodueTEo=deltaTEo*100/totalDeltaTEo
				print(paste0("Weight of Delta TE order(",TEo,") on sum of TEorders differences: ",percSumDeltaTEodueTEo,sp))
				tableTibble=rbind(tableTibble,c(paste0("Weight of Delta TE order(",TEo,") on sum of TEorders differences"),percSumDeltaTEodueTEo,sp,refSpecieForAllvsAll))


			}

		}
	}
}

tableTibble=tableTibble %>% filter(is.na(title)!=T)
#complete diagonal with 0s
for(sp in speciesOrder){
	tableTibble=rbind(tableTibble,c("Delta genome size",0,sp,sp))
	for(TEo in TEorders){
		tableTibble=rbind(tableTibble,c(paste0("Weight of Delta TE order(",TEo,") on sum of TEorders differences"),0,sp,sp))
		tableTibble=rbind(tableTibble,c(paste0("Weight of Delta TE order(",TEo,") on Delta genome size"),0,sp,sp))
		tableTibble=rbind(tableTibble,c(paste0("Delta TE order(",TEo,")"),0,sp,sp))
	}

}

write.table(tableTibble,file=paste0(outprefix,".TE-contribution-to-delta-genome-size.tab"),quote=F,col.names=F,row.names=F,sep="\t")


# now we plot all these numbers.
pdf(paste0(outprefix,".heatmaps.pdf"),width=6,height=6)
tableTibble$specie=factor(tableTibble$specie,levels=speciesOrder)
tableTibble$refSpecie=factor(tableTibble$refSpecie,levels=speciesOrder)


# eliminamos la mitad del heatmap
#cross=crossing(speciesOrder,speciesOrder) %>% mutate(comparisonID=paste0(`speciesOrder...1`,"_and_",`speciesOrder...2`))
#tableTibble=tableTibble %>% mutate(comparisonID=paste0(refSpecie,"_and_",specie))
library(reshape2)

tableTibble.gz.m=	tableTibble %>%
			filter(title=="Delta genome size") %>%
			mutate(value=as.numeric(value)) %>%
			select(refSpecie,specie,value) %>% 
			acast(specie ~ refSpecie)
# dejamos el heatmap triangular
tableTibble.gz.m[lower.tri(tableTibble.gz.m)]=NA
tableTibble.gz=tableTibble.gz.m %>% melt %>% mutate(`specie`=`Var1`,`refSpecie`=`Var2`) %>% select(-Var1,-Var2) %>% dplyr::filter(is.na(`value`)==F)
			
#tableTibble %>% filter(title=="Delta genome size") %>% mutate(value=as.numeric(value)) %>% 
tableTibble.gz %>% 
	ggplot(aes(x=`specie`,y=`refSpecie`,fill=`value`))+geom_tile()+geom_text(aes(label=format(`value`/1000000,digit=2), angle=45),size=3)+theme_minimal()+scale_fill_gradient2(mid = "#ffeda0", high = "#f03b20",low="#91cf60")+
	theme(axis.text.x=element_text(angle=45,hjust=0,vjust=1),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
	      plot.caption = element_text(hjust=0.5, size=rel(1.2)),
	      axis.title.x=element_blank(),
	      legend.position="none")+
ggtitle(expression(Delta*" Genome size"))+
#ggtitle(paste0("\u0394 Genome size\n"))+
labs(caption="Number of Mbps")+
xlab("Specie")+ylab("Reference specie")+scale_x_discrete(position = "top")
dev.off()

pdf(paste0(outprefix,".heatmaps-TEo.pdf"),width=12,height=6)

for(TEo in TEorders){
	tableTibble.TEo.m=tableTibble %>%
                filter(title==paste0("Delta TE order(",TEo,")")) %>%
                mutate(value=as.numeric(value)) %>% 
		select(refSpecie,specie,value) %>% 
		acast(specie ~ refSpecie)

	# dejamos el heatmap triangular
	tableTibble.TEo.m[lower.tri(tableTibble.TEo.m)]=NA
	tableTibble.TEo=tableTibble.TEo.m %>% melt %>% mutate(`specie`=`Var1`,`refSpecie`=`Var2`) %>% select(-Var1,-Var2) %>% dplyr::filter(is.na(`value`)==F)

	# absolute heatmap
	p.teo.1=tableTibble.TEo %>% 
		ggplot(aes(x=`specie`,y=`refSpecie`,fill=`value`))+geom_tile()+geom_text(aes(label=format(`value`/1000000,digit=2),angle=45),size=3)+theme_minimal()+scale_fill_gradient2(mid = "#ffeda0",low="#91cf60", high = "#f03b20")+
		theme(axis.text.x=element_text(angle=45,hjust=0,vjust=1),
		      panel.grid.major = element_blank(),
		      panel.grid.minor = element_blank(),
		      plot.caption = element_text(hjust=0.5, size=rel(1.2)),
		      axis.title.x=element_blank(),
		      legend.position="none")+
		labs(caption=expression(Delta*" TE order content (Mbps)"))+
		#labs(caption="Difference in TE order content\n(Mbps)")+
		xlab("Specie")+ylab("Reference specie")+scale_x_discrete(position = "top")


	# weighted heatmap
	tableTibble.TEo.w.m=tableTibble %>%
                filter(title==paste0("Weight of Delta TE order(",TEo,") on Delta genome size")) %>%
                mutate(value=as.numeric(value)) %>% 
		select(refSpecie,specie,value) %>% 
		acast(specie ~ refSpecie)

	# dejamos el heatmap triangular
	tableTibble.TEo.w.m[lower.tri(tableTibble.TEo.w.m)]=NA
	tableTibble.TEo.w=tableTibble.TEo.w.m %>% melt %>% mutate(`specie`=`Var1`,`refSpecie`=`Var2`) %>% select(-Var1,-Var2) %>% dplyr::filter(is.na(`value`)==F)

	#p.teo.2=tableTibble %>% filter(title==paste0("Weight of Delta TE order(",TEo,") on Delta genome size")) %>% mutate(value=as.numeric(value)) %>% 
	p.teo.2=tableTibble.TEo.w %>%
		ggplot(aes(x=`specie`,y=`refSpecie`,fill=`value`))+geom_tile()+geom_text(aes(label=format(`value`,digit=2),angle=45),size=3)+theme_minimal()+scale_fill_gradient2(mid = "#f5f5f5",high="#d8b365", low = "#5ab4ac")+
		theme(axis.text.x=element_text(angle=45,hjust=0,vjust=1), 
		      panel.grid.major = element_blank(), 
		      panel.grid.minor = element_blank(),
		      plot.caption = element_text(hjust=0.5, size=rel(1.2)),
		      axis.title.x=element_blank(),
		      legend.position="none")+
#		ggtitle("% over the difference in genome sizes")+
		#labs(caption=paste0("Weight of TE order over the\ndifference in genome sizes"))+
		labs(caption=expression("Weight of "*Delta*" TE order over "*Delta*" genome size (%)"))+
		xlab("Specie")+ylab("Reference specie")+scale_x_discrete(position = "top")
	#show(p.teo.1)
	#show(p.teo.2)
	grid.arrange(	p.teo.1,p.teo.2,ncol=2,
			top=textGrob(TEo,just = "centre",gp=gpar(fontsize=25,font=3)))
}


dev.off()



save.image(paste0(outprefix,".RData"))
