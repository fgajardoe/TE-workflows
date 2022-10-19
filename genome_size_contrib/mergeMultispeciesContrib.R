
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



pairedFilesArgs=args[5:length(args)]



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
d=rbind(d,d.coding)


d.order=d$TEorder %>% unique
d.order=d.order[c(
			    6, #coding
			    1, #dna core
			    2, #line
			    5, #sine
			    3, #ltr
			    4)] #other

d.order.pal=c(
		   "#ff6f00ff", #coding
		   "#8a4198ff", #dna core
		   "#c71000ff", #line
		   "#008ea0ff", #sine
		   "#5a9599ff", #ltr
		   "#d1cfa1ff") #other

d$specie=factor(d$specie,levels=rev(speciesOrder),ordered=T)
d$TEorder=factor(d$TEorder,ordered=T,levels=as.character(d.order))


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


#d.tesuperfam.order.pal=c(
#		   #"#ff6f00ff", #coding
#		   "#8a4198ff", #dna core
#		   "#c71000ff", #line
#		   "#008ea0ff", #sine
#		   "#5a9599ff", #ltr
#		   "#d1cfa1ff") #other

d.tesuperfam$TEsuperfam=factor(d.tesuperfam$TEsuperfam,ordered=T,levels=as.character(d.tesuperfam.order))


# plotting

library(ggsci)
library(RColorBrewer)
pdf(paste0(outprefix,".genomeSizePerSpecie.pdf"),width=12,height=7)
d %>% ggplot(aes(x=`Total bps`/1000000000,y=specie, fill=TEorder))+geom_bar(stat="identity",position=position_stack(reverse = TRUE)) +theme_classic()+xlab("Gbps")+scale_fill_manual(values=c(d.order.pal))
d.cons %>% ggplot(aes(x=`Total bps`/1000000000,y=specie, fill=TEorder,group=Conservation))+geom_bar(stat="identity",position=position_stack(reverse = TRUE))+facet_grid(~ Conservation) +theme_classic()+xlab("Gbps")+scale_fill_manual(values=c(d.cons.order.pal))



n_superfams=d.tesuperfam %>% pull(TEsuperfam) %>% unique %>% length
pal.tesuperfams=colorRampPalette(brewer.pal(9,name="Set1"))(n_superfams)
pal.tesuperfams2=colorRampPalette(brewer.pal(8,name="Dark2"))(n_superfams)
d.tesuperfam %>% ggplot(aes(x=`Total bps`/1000000000,y=specie, fill=TEsuperfam,group=Conservation))+geom_bar(stat="identity",position="fill")+
													     facet_grid(~ Conservation)+
													     theme_classic()+
													     theme(legend.position="bottom")+
													     xlab("Proportion of base pairs")+
													     scale_fill_manual(values = pal.tesuperfams)

#d.tesuperfam %>% ggplot(aes(x=`Total bps`/1000000000,y=specie, fill=TEsuperfam,group=Conservation))+geom_bar(stat="identity",position="fill")+
#													     facet_grid(~ Conservation)+
#													     theme_classic()+
#													     theme(legend.position="bottom")+
#													     xlab("Gbps")+
#													     scale_fill_manual(values = pal.tesuperfams2)
#+scale_fill_rickandmorty() # +scale_fill_manual(values=c(d.cons.order.pal))
													     

												     #position_stack(reverse = TRUE))
dev.off()


















save.image(paste0(outprefix,".RData"))
