library(tidyverse)
library(grid)
library(gridExtra)

args=commandArgs(T)

loadRData=function(fileName){
	load(fileName)
	get(ls()[ls() != "fileName"])
}


logfcPath=as.character(args[1])
kimuraPath=as.character(args[2])
logfc.plots=loadRData(logfcPath)
kimura.plots=loadRData(kimuraPath)
outprefix=as.character(args[3])

catNames=c(names(logfc.plots),names(kimura.plots)) %>% tibble %>% mutate(v=1)
colnames(catNames)=c("id","v")
matchingNames=catNames %>% group_by(id) %>% summarise(t=sum(v)) %>% filter(t==2) %>% pull(id)

arrangeOp="h"
#arrangeOp="v"

#ncols=5
#if(arrangeOp=="h"){
	ncols=length(matchingNames)+1
	nrows=2
#}
#if(arrangeOp=="v"){
#	nrows=length(matchingNames)+1
#	ncols=2
#}

unitSize=5

if(arrangeOp=="h"){
	pdf(paste0(outprefix,".LogFC_kimura_plots.pdf"),width=unitSize*ncols,
	    						height=unitSize*nrows)
}
if(arrangeOp=="v"){
	pdf(paste0(outprefix,".LogFC_kimura_plots.pdf"),height=unitSize*ncols,
	    						width=unitSize*nrows)
}
if(arrangeOp != "h" && arrangeOp != "v"){
	error("the value given to `arrangeOp` is not valid.")
}

p.l=logfc.plots[matchingNames]
p.k=kimura.plots[matchingNames]

p.par=vector(mode="list",length(matchingNames))
names(p.par)=matchingNames
for(TEo in matchingNames){
	p.l[[TEo]]=p.l[[TEo]]+theme(title=element_text(size=20),
				    axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=15),
				    axis.text.y=element_text(size=15),
				    axis.title.x=element_blank())
	p.k[[TEo]]=p.k[[TEo]]+theme(axis.title.y=element_text(size=18),
				    title=element_text(size=20),
				    #axis.title.y=element_text(size=18),
				    axis.title.x=element_text(size=18),
				    axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=15),
				    axis.text.y=element_text(size=15))
	
	if(arrangeOp=="h"){
	p.par[[TEo]]=grid.arrange(p.l[[TEo]],p.k[[TEo]],ncol=1)
	}
	if(arrangeOp=="v"){
	p.l[[TEo]]=p.l[[TEo]]+theme(axis.text.x=element_blank())
	p.par[[TEo]]=grid.arrange(p.l[[TEo]],p.k[[TEo]],nrow=1)
	}
}

p.norep=logfc.plots[["No-repetitive"]]+theme(title=element_text(size=20),
				    axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=15),
				    axis.text.y=element_text(size=15),
				    axis.title.x=element_blank())

p.legend=logfc.plots[["Lengend"]]

if(arrangeOp=="h"){
	p.special=grid.arrange(p.norep,p.legend,ncol=1)
}
if(arrangeOp=="v"){
	p.special=grid.arrange(p.norep,p.legend,nrow=1)
}
p.par[["Special"]]=p.special

if(arrangeOp=="h"){
do.call("grid.arrange", c(p.par, ncol=ncols))
}
if(arrangeOp=="v"){
do.call("grid.arrange", c(p.par, nrow=ncols))
}
dev.off()

save.image(paste0(outprefix,".LogFC_kimura_plots.RData"))
