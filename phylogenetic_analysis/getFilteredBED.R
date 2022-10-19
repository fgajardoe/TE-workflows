library(tidyverse)
library(rtracklayer)
library(GenomicFeatures)

loadRData <- function(fileName){
#loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
}


args=commandArgs(T)
rmgrPath=as.character(args[1])
rm.gr=loadRData(rmgrPath)
specialFamiliesAusStr=as.character(args[2])
specialFamiliesAus=strsplit(specialFamiliesAusStr,",")[[1]]


kimura.cutoff=args[3]
if(kimura.cutoff=="max"){
	maxk=rm.gr %>% as.data.frame %>% as_tibble(rownames="fragmentID") %>% filter(is.na(Kimura)==F) %>% pull(Kimura) %>% max %>% as.numeric()
	print(paste0("kimura.cutoff set to 'max'... the max value is ", maxk))
	kimura.cutoff=maxk
}else{
	kimura.cutoff=as.numeric(kimura.cutoff)
}


width.limit=as.numeric(args[4])
#width.limit=read.table(args[4]) %>% as.numeric()

prefix=as.character(args[5])


# exportamos bed para obtener seqs
rm.gr %>% as.data.frame %>% as_tibble(rownames="fragmentID") %>%
	mutate(Kimura=as.numeric(Kimura)) %>% 
	filter(is.na(Kimura)==F) %>% 
	filter(TEsuperfam %in% specialFamiliesAus, width>width.limit) %>% 
	mutate(name=fragmentID) %>% 
#	mutate(name=paste0(prefix,"_",seqnames,"_",start,"_",end)) %>% 
	dplyr::select(seqnames,start,end,name,Kimura,strand) %>% 
	filter(Kimura<kimura.cutoff) %>% 
	write.table(file=paste0(prefix,".filtered-by-length-and-kimura.bed"),quote=F,row.names=F,sep="\t",col.names=F)

save.image(paste0(prefix,".getFilteredBED.RData"))
