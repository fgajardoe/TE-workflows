




library(tidyverse)



args=commandArgs(T)
mergedTable=read.table(args[1],sep="\t",quote="",fill=T,header=T) %>% tibble
pantherSp=read.table(args[2],sep="\t",quote="",fill=T,col.names=c("Gene","PantherFam","Description","Evalue","Score","Coords")) %>% tibble
pantherdb=as.character(args[3])
pvalCutoff=as.numeric(args[4])
outprefix=as.character(args[5])




# panther db
if(pantherdb=="download"){
	print("Generating object with PANTHER FAMILY - GO SLIM ID equivalence ...")
	library(PANTHER.db)
	GO_IDS=keys(PANTHER.db,keytype="GOSLIM_ID")
	panther.go.v=select(PANTHER.db, keys=GO_IDS, columns=c("FAMILY_ID"), keytype="GOSLIM_ID") %>% tibble %>% pull(GOSLIM_ID,FAMILY_ID)
	save(panther.go.v,file="panther_GO_FAM_equivalence.RData")
}else{
	print("Object with PANTHER FAMILY - GO SLIM ID equivalence provided :) ")
	load(pantherdb)
}
pantherSp.v=pantherSp %>% pull(PantherFam, Gene)

library(clusterProfiler)

# universes
allGenome=pantherSp %>% dplyr::select(PantherFam, Gene) %>% mutate(term=panther.go.v[PantherFam],gene=Gene) %>% filter(is.na(term)==F) %>% dplyr::select(term,gene)
allTandem=mergedTable %>% dplyr::select(Gene,Significant) %>% mutate(PantherFam=pantherSp.v[Gene]) %>% mutate(GO_ID=panther.go.v[PantherFam]) %>% filter(is.na(GO_ID)==F) %>% mutate(term=GO_ID, gene=Gene) %>% dplyr::select(term,gene)

# sets
setHighlyTandem=mergedTable %>% dplyr::select(Gene,Significant) %>% mutate(PantherFam=pantherSp.v[Gene]) %>% mutate(GO_ID=panther.go.v[PantherFam]) %>% filter(is.na(GO_ID)==F) %>% filter(Significant==T) %>% mutate(term=GO_ID, gene=Gene) %>% pull(gene)
setAllTandem=mergedTable %>% dplyr::select(Gene,Significant) %>% mutate(PantherFam=pantherSp.v[Gene]) %>% mutate(GO_ID=panther.go.v[PantherFam]) %>% filter(is.na(GO_ID)==F) %>% mutate(term=GO_ID, gene=Gene) %>% pull(gene)

setGSEA_HighlyTandem=mergedTable %>% dplyr::select(Gene,sum_copy_number,Significant) %>% mutate(PantherFam=pantherSp.v[Gene]) %>% mutate(GO_ID=panther.go.v[PantherFam]) %>% filter(is.na(GO_ID)==F) %>% filter(Significant==T) %>% dplyr::select(Gene,sum_copy_number) %>% mutate(rank=rank(sum_copy_number,ties.method = "random")) %>% arrange(desc(rank)) %>% pull(rank,Gene)
setGSEA_AllTandem=mergedTable %>% dplyr::select(Gene,sum_copy_number,Significant) %>% mutate(PantherFam=pantherSp.v[Gene]) %>% mutate(GO_ID=panther.go.v[PantherFam]) %>% filter(is.na(GO_ID)==F) %>% dplyr::select(Gene,sum_copy_number) %>% mutate(rank=rank(sum_copy_number,ties.method = "random")) %>% arrange(desc(rank)) %>% pull(rank,Gene)

setGSEA_percGeneLen_HighlyTandem=mergedTable %>% dplyr::select(Gene,perc_of_gene_length_overlapping_TEs.r,Significant) %>% mutate(PantherFam=pantherSp.v[Gene]) %>% mutate(GO_ID=panther.go.v[PantherFam]) %>% filter(is.na(GO_ID)==F) %>% filter(Significant==T) %>% dplyr::select(Gene,perc_of_gene_length_overlapping_TEs.r) %>% mutate(rank=rank(perc_of_gene_length_overlapping_TEs.r,ties.method = "random")) %>% arrange(desc(rank)) %>% pull(rank,Gene)
setGSEA_percGeneLen_AllTandem=mergedTable %>% dplyr::select(Gene,perc_of_gene_length_overlapping_TEs.r,Significant) %>% mutate(PantherFam=pantherSp.v[Gene]) %>% mutate(GO_ID=panther.go.v[PantherFam]) %>% filter(is.na(GO_ID)==F) %>% dplyr::select(Gene,perc_of_gene_length_overlapping_TEs.r) %>% mutate(rank=rank(perc_of_gene_length_overlapping_TEs.r,ties.method = "random")) %>% arrange(desc(rank)) %>% pull(rank,Gene)

# enrichments
enrichment_AllTandem=enricher(setAllTandem,TERM2GENE=allGenome,pvalueCutoff=pvalCutoff)
enrichment_HighlyTandemVsGenome=enricher(setHighlyTandem,TERM2GENE=allGenome,pvalueCutoff=pvalCutoff)
enrichment_HighlyTandemVsAllTandem=enricher(setHighlyTandem,TERM2GENE=allTandem,pvalueCutoff=pvalCutoff)

enrichmentGSEA_AllTandem=GSEA(setGSEA_AllTandem,TERM2GENE=allGenome,pvalueCutoff=pvalCutoff)
enrichmentGSEA_HighlyTandemVsGenome=GSEA(setGSEA_HighlyTandem,TERM2GENE=allGenome,pvalueCutoff=pvalCutoff)
enrichmentGSEA_HighlyTandemVsAllTandem=GSEA(setGSEA_HighlyTandem,TERM2GENE=allTandem,pvalueCutoff=pvalCutoff)

enrichmentGSEA_percGeneLen_AllTandem=GSEA(setGSEA_percGeneLen_AllTandem,TERM2GENE=allGenome,pvalueCutoff=pvalCutoff)
enrichmentGSEA_percGeneLen_HighlyTandemVsGenome=GSEA(setGSEA_percGeneLen_HighlyTandem,TERM2GENE=allGenome,pvalueCutoff=pvalCutoff)
enrichmentGSEA_percGeneLen_HighlyTandemVsAllTandem=GSEA(setGSEA_percGeneLen_HighlyTandem,TERM2GENE=allTandem,pvalueCutoff=pvalCutoff)

print("--------------------------------------------")
print(paste0("Summary of significant enrichments (Adj. p-value < ",pvalCutoff,")"))
print("--------------------------------------------")
n1=enrichment_AllTandem@result$qvalue %>% .[.<pvalCutoff] %>% length
n2=enrichment_HighlyTandemVsGenome@result$qvalue %>% .[.<pvalCutoff] %>% length
n3=enrichment_HighlyTandemVsAllTandem@result$qvalue %>% .[.<pvalCutoff] %>% length
n4=enrichmentGSEA_AllTandem@result$qvalue %>% .[.<pvalCutoff] %>% length
n5=enrichmentGSEA_HighlyTandemVsGenome@result$qvalue %>% .[.<pvalCutoff] %>% length
n6=enrichmentGSEA_HighlyTandemVsAllTandem@result$qvalue %>% .[.<pvalCutoff] %>% length
n7=enrichmentGSEA_percGeneLen_AllTandem@result$qvalue %>% .[.<pvalCutoff] %>% length
n8=enrichmentGSEA_percGeneLen_HighlyTandemVsGenome@result$qvalue %>% .[.<pvalCutoff] %>% length
n9=enrichmentGSEA_percGeneLen_HighlyTandemVsAllTandem@result$qvalue %>% .[.<pvalCutoff] %>% length
print("Gene set enrichment method (enrich)")
print(paste0("All tandem vs genome: ",n1))
print(paste0("Highly tandem vs genome: ",n2))
print(paste0("Highly tandem vs All tandem: ",n3))
print("Ranked gene set enrichment method (GSEA ranked by copy number)")
print(paste0("All tandem vs genome: ",n4))
print(paste0("Highly tandem vs genome: ",n5))
print(paste0("Highlt tandem vs All tandem: ",n6))
print("Ranked gene set enrichment method (GSEA ranked by % of gene length)")
print(paste0("All tandem vs genome: ",n7))
print(paste0("Highly tandem vs genome: ",n8))
print(paste0("Highly tandem vs All tandem: ",n9))

save.image(paste0(outprefix,".enrichment.RData"))
