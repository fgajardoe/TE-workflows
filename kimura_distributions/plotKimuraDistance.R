## Init
library(tidyverse)
library(RColorBrewer)
library(gridExtra)
source("plotKimuraDistance.functions.R")


system("pwd")

# functions

addOther=function(d,genomeSizes){
	species.nr=d %>% pull(Specie) %>% unique
	for(sp in species.nr){
		gSize=as.numeric(genomeSizes[sp])
		rSize=d %>% filter(Specie==sp) %>% pull(`Genome Size Contribution`) %>% as.numeric %>% sum
		f=c("Non-repetitive",as.character(sp),gSize-rSize)
		d=rbind(d,f)
	}
	return(d)

}


plotGenomeSizeContribution=function(d.bed.lst, TEclassification=NULL, genomeSizes, TEclassificationLevel="TEorder"){

        for(sp in names(d.bed.lst)){
		sp.bed=d.bed.lst[[sp]]
		#read.table("Nematolebias_whitei.bed",sep="\t",comment.char="@",col.names=c("seqname","start","end","name","score","strand","kimura")) 
		d.bed.lst[[sp]]$fragmentSize=d.bed.lst[[sp]]$end-d.bed.lst[[sp]]$start
                d.bed.lst[[sp]]$Specie=as.character(sp)

        }

	d=do.call("rbind",d.bed.lst)

	if(!is.null(TEclassification)){
		TEorder.v=TEclassification %>% pull(TEorder,family_id)
		TEsuperfam.v=TEclassification %>% pull(TEsuperfam,family_id)
		d$TEorder=TEorder.v[d$name]
		d$TEsuperfam=TEsuperfam.v[d$name]
	}

	d=d%>% mutate(TEorder=ifelse(is.na(TEorder),"Unknown",TEorder), TEsuperfam=ifelse(is.na(TEsuperfam),"Unknown",TEsuperfam)) %>% aggregate(fragmentSize ~ TEorder + Specie, ., sum) %>% mutate(`Genome Size Contribution`=fragmentSize) %>% select(-fragmentSize)
	d=addOther(d,genomeSizes)
	d$`Genome Size Contribution`=as.numeric(d$`Genome Size Contribution`)
	return(d)

}

plotKimuraCompareBetweenSpecies=function(d.bed.lst, TEclassification=NULL){
        for(sp in names(d.bed.lst)){
                d.bed.lst[[sp]]$Specie=as.character(sp)
        }
        d=do.call("rbind",d.bed.lst)

        if(!is.null(TEclassification)){
                TEorder.v=TEclassification %>% pull(TEorder,family_id)
                TEsuperfam.v=TEclassification %>% pull(TEsuperfam,family_id)
                d$TEorder=TEorder.v[d$name]
                d$TEsuperfam=TEsuperfam.v[d$name]
        }
        return(d)
}

g_legend<-function(a.gplot){
        tmp <- ggplot_gtable(ggplot_build(a.gplot))
        leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
        legend <- tmp$grobs[[leg]]
        return(legend)
}

# parametros de entrada
args=commandArgs(trailingOnly=TRUE)

TElevel=as.character(args[9])
normalizeByTotal=as.character(args[10])

outprefix=args[11]
inTable=args[1]
eqFile=args[2]

inTable=read.table(as.character(inTable),sep="\t",col.names=c("specie","bedFile","genomeSizeFile"))
bedFiles=inTable %>% pull(bedFile,specie)


# read and sum fragments(contigs, scaffolds, etc) sizes
genomeSizes=c()
genomeSizesFiles=inTable %>% pull(genomeSizeFile,specie)
for(sp in names(genomeSizesFiles)){
	genomeSizeFile=genomeSizesFiles[sp]
	gz.tibble=read.table(genomeSizeFile, sep="\t",col.names=c("seqname","length")) %>% tibble
	genomeSizes[sp]=gz.tibble %>% pull(length) %>% sum 
}


# Parameters Kimura plot
title=as.character(args[5])
ymax=as.numeric(args[4]) #150000
xmax=as.numeric(args[3]) #50
#TEorderSelection=T
TEorderSelection=args[6]
#pal="Dark2"
pal.nlimit=8
plot.ncol=5


# Parameters Genome size contribution plot
speciesOrder=args[7]
speciesOrder=strsplit(speciesOrder,",")[[1]]
#speciesOrder=c("Austrolebias_charrua","Nematolebias_whitei","Cynopoecilus_melanotaenia","Austrofundulus_limnaeus","Kryptolebias_marmoratus","Nothobranchius_furzeri","Oryzias_latipes")
speciesOrder=rev(speciesOrder)

TEorderOrder=args[8]
TEorderOrder=strsplit(TEorderOrder,",")[[1]]
TEorderOrder=c(TEorderOrder,"Non-repetitive")
#TEorderOrder=c("LINE","SINE","LTR","Retroposon","DNA","RC","Satellite","rRNA","tRNA","Unknown","Non-repetitive")
TEorderOrder=rev(TEorderOrder)
pal.genomeSizeContrib="Spectral"
pal.genomeSizeContrib.nlimit=9

# MAIN

# for que lee los beds
d.lst=vector(mode="list",length(inTable$specie))
names(d.lst)=inTable$specie

# alias for species names
sp.alias=rep(0,length(inTable$specie))
names(sp.alias)=inTable$specie


for(sp in inTable$specie){
	bed=read.table(bedFiles[sp],sep="\t",comment.char="@",col.names=c("seqname","start","end","name","score","strand","kimura")) %>% as_tibble %>% filter(kimura!="ND") %>% mutate(kimura=as.numeric(kimura))
	d.lst[[sp]]=bed


	sp.alias[sp]=gsub("\\w+_",paste0(substr(sp,0,1),". "),sp)

}


# load TE equivalence table
eq=read.table(as.character(eqFile),col.names=c("family_id","TEsuperfam","TEorder")) %>% as_tibble


# Genome size contribution

d.genomeSizeContrib=plotGenomeSizeContribution(d.lst,eq,genomeSizes)
pal.genomeSizeContrib=colorRampPalette(brewer.pal(pal.genomeSizeContrib.nlimit, pal.genomeSizeContrib))(length(TEorderOrder))

d.genomeSizeContrib$Specie=factor(d.genomeSizeContrib$Specie,levels=speciesOrder)
d.genomeSizeContrib$TEorder=factor(d.genomeSizeContrib$TEorder,levels=TEorderOrder)


p.genomeSize=d.genomeSizeContrib %>% ggplot(aes(x=`Genome Size Contribution`/1000000000, y=Specie,fill=TEorder))+geom_bar(stat="identity",position="stack")+ggtitle(paste0(title," - Genome size contribution"))+theme_classic()+xlab("Gbp")+scale_fill_manual(values=rev(pal.genomeSizeContrib))


# plot Kimura distributions

d=plotKimuraCompareBetweenSpecies(d.lst,eq) # where `d.bed.lst` is a list of "BED" dataframes `cmel.bed` and `nwhi.bed`
						# and `eq` is a tibble with `family_id`, `TEorder` and `TEsuperfam` columns.
d$kimura_int=as.integer(d$kimura)





# Calculamos el numero de fragmentos por kimura
if(TElevel=="TEorder"){
	d.counts=d %>% count(kimura_int,TEorder,Specie)
	d.text=d %>% group_by(TEorder,Specie) %>% summarise(bps=sum(abs(end-start)), Specie=unique(Specie)) %>% mutate(perc=100*bps/genomeSizes[Specie])
}
if(TElevel=="TEsuperfamily"){
	d.counts=d %>% count(kimura_int,TEsuperfam,Specie)
	d.text=d %>% group_by(TEsuperfam,Specie) %>% summarise(bps=sum(abs(end-start)), Specie=unique(Specie)) %>% mutate(perc=100*bps/genomeSizes[Specie])
}

#colnames(d.counts)=c("kimura_int","TEclass","Specie")


if(any(normalizeByTotal)){
	maxValue.v=d.counts %>% group_by(Specie) %>% summarise(maxValue=sum(n)) %>% pull(maxValue,Specie)
	ymax=(ymax/max(maxValue.v))*100
	d.counts$n=(d.counts$n/maxValue.v[d.counts$Specie])*100
}


#TEsuperfam.nBps.v=d %>% group_by(TEsuperfam) %>% summarise(nBpsTotalSuperfam=sum(abs(end-start))) %>% pull(nBpsTotalSuperfam,TEsuperfam)



TEorders=strsplit(TEorderSelection,",")[[1]]


library(ggsci)

#pal=colorRampPalette(brewer.pal(pal.nlimit, pal))(length(d.lst))
pal=pal_uchicago(palette="default",alpha=0.7)(length(d.lst))
rm(d.lst)

p.lst=vector(mode="list",length(TEorders))
names(p.lst)=TEorders



# flag para hacer la leyenda una pura vez
renderLegendPlot=T


alias.sp=c(	"Austrolebias_charrua"="A. charrua",
		"Cynopoecilus"		
)


rm(d)
gc()




# Por cada TEorder o TEsuperfam
for(TEo in TEorders){
	print(paste0("Plotting ",TEo,"..."))
	if(TElevel=="TEorder"){
		d.text.TEo=d.text %>% filter(TEorder==TEo)
		d.text.TEo=d.text.TEo %>% ungroup %>% select(Specie,bps,perc)
		d.text.TEo=d.text.TEo %>% mutate(x=18,y=ymax-(seq(ymax/10,ymax*10,by=ymax/10)[1:NROW(.)]))

		nTEoInData=d.counts %>% filter(TEorder==TEo) %>% NROW %>% as.numeric
		if(nTEoInData==0){
			print(paste0("TE order ",TEo," not found in data. Skipping"))
			p.lst[[TEo]]=ggplot()+geom_text(aes(x=1,y=1,label="Order\nnot found in the data"),size=6)+theme_void()+ggtitle(as.character(TEo))
			next
		}

		p.lst[[TEo]]= d.counts %>%
				filter(TEorder==TEo) %>% 
				ggplot(.,aes(x=kimura_int,y=n,colour=Specie))+
				geom_line(size=1.1)+
				xlim(values=c(0,xmax))+ylim(values=c(0,ymax))+
				theme_classic()+
				scale_color_manual(values=pal)+
				ggtitle(as.character(TEo))+
				theme(legend.position="none")+
				ylab("Count")+xlab("Kimura distance")+
			#	geom_text(data=d.text.TEo,aes(x=x,y=y,label=paste0(bps," bps (",perc,"%)"), colour=Specie), position=position_stack(vjust = 0.5), inherit.aes = T)
				geom_text(data=d.text.TEo,aes(x=x,y=y,label=paste0(format(bps/1000000,digit=2)," Mbps (",format(perc,digit=2),"%) "), colour=Specie), 
			  	  inherit.aes = F, hjust=0,size=5)

		if(any(renderLegendPlot)){
			p.legend=d.counts %>% filter(TEorder==TEo) %>% ggplot(.,aes(x=kimura_int,y=n,colour=Specie))+geom_point()+geom_line()+xlim(values=c(0,xmax))+ylim(values=c(0,ymax))+ theme_classic()+scale_color_manual(values=pal)+ggtitle(as.character(TEo))+ylab("Count")+xlab("Kimura distance")
			renderLegendPlot=F
		}
	}
	if(TElevel=="TEsuperfamily"){
		d.text.TEo=d.text %>% filter(TEsuperfam==TEo)
		d.text.TEo=d.text.TEo %>% ungroup %>% select(Specie,bps,perc)
		#d.text.TEo=d.text.TEo %>% mutate(x=18,y=200000-(seq(20000,2000000,by=20000)[1:NROW(.)]))
		d.text.TEo=d.text.TEo %>% mutate(x=18,y=ymax-(seq(ymax/10,ymax*10,by=ymax/10)[1:NROW(.)]))

		nTEoInData=d.counts %>% filter(TEsuperfam==TEo) %>% NROW %>% as.numeric
		if(nTEoInData==0){
			print(paste0("TE superfamily ",TEo," not found in data. Skipping"))
			p.lst[[TEo]]=ggplot()+geom_text(aes(x=1,y=1,label="Superfamily\nnot found in the data"),size=6)+theme_void()+ggtitle(as.character(TEo))
			next
		}

		p.lst[[TEo]]= d.counts %>% 
			filter(TEsuperfam==TEo) %>% 
			ggplot(.,aes(x=kimura_int,y=n,colour=Specie))+
			geom_line(size=1.1)+
			xlim(values=c(0,xmax))+ylim(values=c(0,ymax))+ 
			theme_classic()+
			scale_color_manual(values=pal)+
			ggtitle(as.character(TEo))+
			theme(legend.position="none")+
			ylab("Count")+xlab("Kimura distance")+
			geom_text(data=d.text.TEo,aes(x=x,y=y,label=paste0(format(bps/1000000,digit=2)," Mbps (",format(perc,digit=2),"%) "), colour=Specie), 
				  #position=position_stack(vjust = 0.5),
			  	  inherit.aes = F, hjust=0,size=5)
	if(any(renderLegendPlot)){	
			p.legend=d.counts %>% filter(TEsuperfam==TEo) %>% ggplot(.,aes(x=kimura_int,y=n,colour=Specie))+geom_point()+geom_line()+xlim(values=c(0,xmax))+ylim(values=c(0,ymax))+ theme_classic()+scale_color_manual(values=pal)+ggtitle(as.character(TEo))+ylab("Count")+xlab("Kimura distance")
			renderLegendPlot=F
		}
	}
}

print("building plot legend ...")
p.legend=g_legend(p.legend)
p.lst[["legend"]]=p.legend

save.image("test.RData")

library(grid)

if(TElevel=="TEorder"){
TElevelStr="TE orders"
} else {
TElevelStr="TE superfamilies"
}

title=paste0(title," - ",TElevelStr," - Kimura distribution ","\n")
title.grob=textGrob(as.character(title),gp=gpar(fontsize=20,font=3))
grid=do.call("arrangeGrob",c(p.lst, ncol=plot.ncol))

#pequenio ajuste para el pdf
#if(TElevel=="TEorder"){
#	page_height=5
#}
#if(TElevel=="TEsuperfamily"){
#	page_height=10
#}
#segundo pequenio ajuste para el pdf
nplots=length(p.lst) %>% as.numeric()
maxNumberOfPlots=1000
multiples_of_5=(1:maxNumberOfPlots)%%5
names(multiples_of_5)=seq(1:maxNumberOfPlots)
multiples_of_5=multiples_of_5 %>% .[.==0]
page_height.v=rep(0,maxNumberOfPlots)
names(page_height.v)=seq(1,maxNumberOfPlots)
z=1
for(i in names(multiples_of_5)){ for(j in seq(1:5)){ page_height.v[z]=i; z=z+1 }  }
page_height=page_height.v[nplots] %>% as.numeric()






# plotting PDF
pdf(paste(outprefix,".kimuraDistrib.pdf",sep=""),height=page_height,width=20)
grid.arrange(grid,top=title.grob)
dev.off()

pdf(paste0(outprefix,".genomeSizeContrib.pdf"),height=5,width=7)
p.genomeSize
dev.off()

save.env=T
if(any(save.env)){
	save.image(paste0(outprefix,".RData"))
}
