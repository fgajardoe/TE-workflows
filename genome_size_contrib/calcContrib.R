args=commandArgs(T)


rmgrPath=as.character(args[1])
load(rmgrPath)


outprefix=as.character(args[3])


library(tidyverse)

rmgr.tibble=rmgr %>% as.data.frame %>% as_tibble(rownames="TE_ID")


ordersToConsider=strsplit(as.character(args[2]),split=",")[[1]]   #c("LINE","SINE","DNA","LTR")

TEsuperfam.nBPs.lst=vector(mode="list",length(ordersToConsider))
names(TEsuperfam.nBPs.lst)=ordersToConsider

p.lst=vector(mode="list",length(ordersToConsider))
names(p.lst)=ordersToConsider


for(o in ordersToConsider){
	print(paste0("working on order ",o))
	TEsuperfam.nBPs.lst[[o]]=rmgr.tibble %>% filter(TEorder==as.character(o)) %>% group_by(TEsuperfam) %>% summarise(nBPs = sum(width)) %>% arrange(desc(nBPs)) %>% mutate(TEorder=as.character(o))
	p.lst[[o]]=TEsuperfam.nBPs.lst[[o]] %>% ggplot(aes(x=nBPs,y=TEsuperfam)) + geom_bar(stat="identity",width=0.8)+theme_classic()+theme(legend.position="none")+xlab("Mbps")	

}
TEsuperfam.nBPs=do.call("rbind",TEsuperfam.nBPs.lst)
TEsuperfam.nBPs$TEorder=factor(TEsuperfam.nBPs$TEorder)
TEsuperfam.nBPs$TEsuperfam=factor(TEsuperfam.nBPs$TEsuperfam)
TEsuperfam.nBPs=TEsuperfam.nBPs %>%  mutate(specie=outprefix)


TEsuperfam.nBPs %>% select(TEsuperfam) %>% mutate(TEsuperfam=as.character(TEsuperfam)) %>% unique %>% write.table(.,file=paste0("TEsuperfams_in_",outprefix,".lst"), row.names=F, col.names=F, quote=F,sep="\t")
save(TEsuperfam.nBPs,file=paste0(outprefix,".RData"))


maxVal=TEsuperfam.nBPs %>% filter(nBPs==max(.$nBPs)) %>% select(nBPs) %>% unlist
x.lim=c(0,maxVal)
for(o in ordersToConsider){
	p.lst[[o]]=p.lst[[o]] + xlim(x.lim) + ggtitle(as.character(o))

}




library(grid)
library(gridExtra)

plot.ncol=2
title=paste0("Genome size contribution of TE superfamilies\nin the ",outprefix," genome")

title.grob=textGrob(as.character(title),gp=gpar(fontsize=20,font=3))
grid=do.call("arrangeGrob",c(p.lst, ncol=plot.ncol))




pdf(paste0(outprefix,".pdf"), height=13, width=13)
grid.arrange(grid,top=title.grob)
dev.off()

