library(tidyverse)
library(GenomicFeatures)
library(rtracklayer)

loadRData <- function(fileName){
#loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
}


args=commandArgs(T)


rmgr=loadRData(args[1])
cons=import(args[2])

centroids=import(args[3])

TEfamStr=as.character(args[4])

#centroidStr="Nematolebias_whitei_TE31915"


outprefix=as.character(args[5])

cons.sizes.tibble=tibble(name=cons %>% names, size=cons %>% width)
cons.len=cons.sizes.tibble %>% filter(name==TEfamStr) %>% pull

centroids.sizes=tibble(name=centroids %>% names, size=centroids %>% width) %>% pull(size,name)

#d=rmgr %>% as.data.frame %>% tibble %>% filter(TEfam==TEfamStr 
d=rmgr %>% as.data.frame %>% as_tibble(rownames="TE_ID") %>% filter(TEfam==TEfamStr) %>% mutate(label=ifelse(TE_ID %in% names(centroids.sizes), TE_ID, ""))

p=d %>% ggplot(aes(x=width))+geom_density()+geom_vline(xintercept=cons.len,colour="red") +theme_classic()

h=0.0001
for(c in names(centroids.sizes)){
	c.len=centroids.sizes[c]
	p=p+geom_vline(xintercept=c.len,colour="blue",alpha=0.5)+geom_text(data=d,aes(x=width,y=h,label=label),angle=90,size=5,check_overlap=F,nudge_x=-10,colour="darkblue")
	#h=h+0.0001
}
pdf(paste0(outprefix,".TEfragmentLengthDistribution.pdf"))
p+ggtitle(outprefix)
dev.off()

#centroid.len=rmgr[centroidStr] %>% data.frame %>% tibble %>% pull(width)





