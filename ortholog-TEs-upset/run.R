library(ggupset)
library(tidyverse)

args=commandArgs(T)
outprefix=as.character(args[1])
TEclass=read.table(as.character(args[2]), sep="\t", col.names=c("TEfam","TEsuperfam","TEorder")) %>% tibble
TEclassSimple=TEclass %>% distinct(TEsuperfam,TEorder) 
TEordersPal=read.table(args[3],sep="\t",comment.char="@",col.names=c("TEorder","colour")) %>% tibble %>% pull(colour,TEorder)
argsLsts=args[4:length(args)]

d.lst=vector(mode="list", length(argsLsts))

species=argsLsts %>% gsub(".TEsuperfams.lst","",x=.)

names(d.lst)=species
for(sp in species){
	d.lst[[sp]]=read.table(paste0(sp,".TEsuperfams.lst"),sep="\t",col.names=c("TEsuperfam")) %>%
		tibble %>%
		mutate(specie=sp)

}
d=do.call("rbind",d.lst)
d.2=d %>% select(`TEsuperfam`,specie) %>% distinct(`TEsuperfam`,specie) %>% group_by(`TEsuperfam`) %>% summarise(speciesArr=list(`specie`))


intersections=as.list(species)
intersections[[length(intersections)+1]]=species
pal.v=read.table("colourFileAustrolebiasPanel.tab",comment.char="@") %>% tibble %>% pull(V2,V1)
pdf(file=paste0(outprefix,".upset.pdf"),width=14,height=4)
d.2 %>% distinct(`TEsuperfam`,speciesArr) %>%
	ggplot(aes(x=speciesArr))+
	geom_bar()+
	#geom_bar(data=TEclassSimple,aes())+
	geom_text(stat='count', aes(label=after_stat(count),angle=90),size=3, hjust=-0.5) +
	scale_x_upset(order_by="degree",reverse=T,sets=species,intersections=intersections)+
	theme_minimal()+
	theme_combmatrix(combmatrix.panel.line.size = 0, combmatrix.label.make_space = FALSE)+
	xlab("Intersections")+
	theme(axis.text.y=element_text(size=12))


library(ComplexUpset)
library(reshape2)
#set_size(8, 3)
d.matrix=d %>% mutate(value=1) %>% dcast(TEsuperfam~specie)
d.matrix[is.na(d.matrix)==T]=0
allowedOrders=c("DNA","LINE","SINE","LTR")
d.matrix=d.matrix %>% tibble %>% left_join(.,TEclassSimple) %>% mutate(TEorder=ifelse(TEorder %in% allowedOrders,TEorder,"Other-TEs")) %>% as.data.frame
upset(d.matrix,
      rev(species),
      name='TEsuperfam', 
      base_annotations=list(
			    `Intersection size`=intersection_size(counts=T,
								  mapping=aes(fill=TEorder))+
			    scale_fill_manual(values=TEordersPal)),
      width_ratio=0.1,
      set_sizes=FALSE,
      sort_intersections_by='degree',
      sort_sets=FALSE)



dev.off()

save.image(paste0(outprefix,".RData"))
