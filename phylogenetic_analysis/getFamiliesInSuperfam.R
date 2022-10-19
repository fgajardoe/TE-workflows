library(tidyverse)

args=commandArgs(T)

load(args[1])
sf=as.character(args[2])
outprefix=args[3]

lst=rmgr %>% as.data.frame %>% as.tibble %>% filter(TEsuperfam==sf) %>% pull(TEfam) %>% unique
tibble(ids=lst) %>% write.table(file=paste0(outprefix,".TEfamsInSuperfam.lst"),quote=F,col.names=F,row.names=F) 
