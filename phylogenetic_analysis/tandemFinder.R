library(tidyverse)
library(rtracklayer)
library(GenomicRanges)

# function
loadRData <- function(fileName){
#loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
}


plotDensity=function(gr,region,wsize=30,wstep=4){


#	tiles=tile(region,wsize)
	tiles=slidingWindows(region,wsize,wstep)[[1]]
	
	tileIDs=seq(1,length(tiles),1) %>% as_tibble %>% mutate(ID=paste0("TILE_",value)) %>% pull(ID)
	names(tiles)=tileIDs
	tiles.all=tibble(ID=names(tiles))
	tiles.cov=findOverlapPairs(gr,tiles)@second %>% names %>% tibble(ID=.) %>% dplyr::count(ID)
	covRow=left_join(tiles.all,tiles.cov,by="ID") %>% mutate(n=ifelse(is.na(n)==T, 0,n)) %>% pull(n)  %>% as.matrix %>% t 
	graphics::image(1:NCOL(covRow),1:NROW(covRow),t(covRow))
#	return(covRow)
#	graphics::image(1:NCOL(covRow),1,covRow)	



}


plotDotplotWithDensity=function(seq.v,wsize,wstep,nmatch,title,gr,region){

	nf=layout(matrix(c(1,2), ncol=1),heights=c(3,2))
	seqinr::dotPlot(seq.v,seq.v,wsize=wsize,wstep=wstep,nmatch = nmatch,main=as.character(title),ylab="",xlab="")
	plotDensity(ins.gr,tile,wsize=wsize,wstep=wstep)

}


# main
args=commandArgs(T)
#rmgr.path="AUS.rmgr.Rdata"
rmgr.path=as.character(args[1])
rmgr=loadRData(rmgr.path)

#gff="aus.genes.gff3"
gff=as.character(args[2])
gff=import(as.character(gff))


#TEfam="ACH_rnd-1_family-146"
TEfam=as.character(args[4])

#genomeSeq="austrolebiask3g3.l10000.fasta"
genomeSeq=as.character(args[3])
genomeSeq=import(genomeSeq)



#genomeSeq.len="austrolebiask3g3.l10000.len"
#genomeSeq.len=read.table(genomeSeq.len) %>% as_tibble %>% pull(V2,V1)

# genera un vector con los largos de los scaffolds a partir de la info del ensamble
genomeSeq.len=genomeSeq@ranges %>% as.data.frame %>% tibble %>% pull(width,names)

#tileWidth=1000
tileWidth=as.numeric(args[5])
#tileWidthMinLength=500
tileWidthMinLength=as.numeric(args[6])


#minNumFragments=2
minNumFragments=as.numeric(args[7]) #for dotplot


# identity limit to consider as part of a cluster
#identityLimit=0.9
identityLimit=as.numeric(args[8])

#outprefix="RTEBovB_Tandems"
outprefix=as.character(args[10])
wrkFolder=paste0("tandemFinderRun_",outprefix)


if(!dir.exists(wrkFolder)){
	print("Creating a folder for this run...")
	system(paste0("mkdir ",wrkFolder))
}else{
	print("Previous run folder detected...")
}


tiles=tileGenome(genomeSeq.len, tilewidth=tileWidth, cut.last.tile.in.chrom=F)
tileIDs=seq(1,length(tiles),1) %>% as_tibble %>% mutate(ID=paste0("TILE_",value)) %>% pull(ID)
names(tiles)=tileIDs
tiles=tiles %>% width %>% unlist %>% .[.>=as.numeric(tileWidthMinLength)] %>% names %>% tiles[.]

rmgr.filtered=rmgr[rmgr$TEfam %in% TEfam]

overlaps=findOverlapPairs(tiles,rmgr.filtered)
names(overlaps)=names(overlaps@first)


tiles.tibble=tibble(tiles=overlaps@first %>% names,insertion=overlaps@second %>% names)
tileNfragments=tiles.tibble %>% group_by(tiles) %>% tally() %>% arrange(desc(n))

# obtenemos un subset de tiles para inspeccionar
tiles.subset=tileNfragments %>% filter(n>1) %>% pull(tiles) %>% tiles[.]



pdf(paste0(outprefix,".dotplots.pdf"),height=14,width=14)


tileScore=vector(mode="list",length(tiles.subset))
names(tileScore)=tiles.subset %>% names
tileNFrag=vector(mode="list",length(tiles.subset))
names(tileNFrag)=tiles.subset %>% names
tileNClust=vector(mode="list",length(tiles.subset))
names(tileNClust)=tiles.subset %>% names

runCDHIT=T
#for each Tile
cdhitExec=as.character(args[9])
#cdhitExec="/opt/cdhit/cd-hit-est"
for(tileName in tiles.subset %>% names){
	print(paste0("working on ",tileName))

	tile=tiles.subset[[tileName]]

	# create insertions FASTA

	ins=tiles.tibble %>% filter(tiles==tileName) %>% pull(insertion)
	ins.gr=rmgr[ins]

	# modificamos limites muy extremos (como start=0 y end=largo de scaffold)
	goodCoords=tibble(name=names(ins.gr),start=start(ins.gr),end=end(ins.gr),seqnames=as.character(seqnames(ins.gr))) %>% mutate(start.good=ifelse(start==0,1,start), end.good=ifelse(end==genomeSeq.len[seqnames],genomeSeq.len[seqnames]-1,end))
	start(ins.gr)=goodCoords %>% pull(start.good)
	end(ins.gr)=goodCoords %>% pull(end.good)



	#si el archivo no existe
	seqPath=paste0(wrkFolder,"/",tileName,".insertions.fasta")
	if(!any(file.exists(seqPath))){
		print("Generating FASTA with insertion sequences.")
		ins.seq=BSgenome::getSeq(genomeSeq,ins.gr)
		export(ins.seq,seqPath)


	}else{
		print("There are [insertions] FASTA files of a previous run. I will use these instead.")
		ins.seq=import(seqPath)
	}
	if(any(runCDHIT)){
		if(!any(file.exists(paste0(wrkFolder,"/",tileName,".cdhit")))){
			system(paste0(cdhitExec," -i ",seqPath," -o ",wrkFolder,"/",tileName,".cdhit"," -c ",identityLimit))
		}else{
			print("There are files of a previous run of CD-HIT. I will use these instead...\nYes, I'm lazy.")
		}
		nClusters=system(paste0("less ",wrkFolder,"/",tileName,".cdhit |grep -c \">\""), intern = TRUE) %>% as.numeric
		nFragments=ins.seq %>% length
		scoreTandem=1-(nClusters/nFragments)

		tileScore[[tileName]]=scoreTandem
		tileNFrag[[tileName]]=nFragments
		tileNClust[[tileName]]=nClusters
	}
	if(nFragments>=minNumFragments){
		#dotplot
		tileSeqPath=paste0(wrkFolder,"/",tileName,".fasta")
		if(!any(file.exists(tileSeqPath))){
			print("Generating FASTA with tile sequence.")
			wSeq=BSgenome::getSeq(genomeSeq, tile)
			export(wSeq,tileSeqPath)
		}else{
			print("There are [tile] FASTA files of a previous run. I will use these instead.")
			wSeq=import(tileSeqPath)
		}
		wSeq.str= wSeq %>% as.character()
		wSeq.v=wSeq.str %>% strsplit("") %>% .[[1]]
		tryCatch(
			 error=function(cnd){
				 print("Some tiles have been skipped, check with warnings()")
				 warning(paste0(tileName," skipped."))
			 },
			 plotDotplotWithDensity(wSeq.v,wsize=30,wstep=4,nmatch=30,tileName,ins.gr,tile))
	}
}

print("CD-HIT finished")

if(any(runCDHIT)){

	tileScore=unlist(tileScore)
	tileNFrag=unlist(tileNFrag)
	tileNClust=unlist(tileNClust)
	tiles.subset.score=tibble(tiles=names(tileScore),nFragments=tileNFrag,nClusters=tileNClust,Score=tileScore)
	tiles.tibble=tiles.tibble %>% left_join(tiles.subset.score,by="tiles") %>% mutate(Score=ifelse(is.na(Score),0,Score), nFragments=ifelse(is.na(nFragments),0,nFragments), nClusters=ifelse(is.na(nClusters),0,nClusters))
}

dev.off()

save(tiles.tibble,file=paste0(outprefix,".tandemFinder.RData"))

#save.image(paste0(outprefix,".RData"))
