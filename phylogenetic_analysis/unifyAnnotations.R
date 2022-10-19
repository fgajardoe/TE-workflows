library(tidyverse)
library(rtracklayer)
library(Gviz)
library(ggsci)

options(ucscChromosomeNames=FALSE)
args=commandArgs(T)

#CDD database
rps=args[1]
rps=read.table(as.character(rps),fill=T,quote="",col.names=c("ORF_ID","session_ordinal","Type","ID_PSSM","from","to","E_Value","bitscore","accession","short_name","incomplete","superfamily_PSSM_ID")) %>% tibble
rps$from=as.numeric(rps$from)
rps$to=as.numeric(rps$to)


#RT only
hmm=args[2]
hmm=read.table(as.character(hmm),fill=T,quote="",col.names=c("ORF_ID","hitIndex","score","bias","c_Evalue","i_Evalue","hmmFrom","hmmTo","hmmIntegrity","aliFrom","aliTo","aliIntegrity","envFrom","envTo","envIntegrity","acc")) %>% tibble 

if(NROW(hmm)!=0){
	hmm=hmm %>% mutate(RT_ID=paste0("RT_",seq(1,NROW(hmm))))
}

#else{
#
#}


bls=args[3]
bls=read.table(as.character(bls),fill=T,quote="",col.names=c("ORF_ID","hitAccession","identity","length","mismatch","gapopen","qstart","qend","start","send","evalue","bitscore","qlen","slen","qcovs","scovs")) %>% tibble


centroidsBed=args[4]

orfsBed=args[5]
orfsBed=read.table(orfsBed,fill=T,quote="",col.names=c("centroidID","start","end","name","strand","score","ORF_ID")) %>% tibble

famCentroidEq=read.table(as.character(args[6]),sep="\t", col.names=c("TEinsertion","Cluster","centroidID")) %>% tibble

load(args[7])
rmgr.tibble=rmgr %>% as.data.frame %>% as_tibble(rownames="TEinsertion")
rm(rmgr)
gc()

# generamos un GenomicRanges con las coordenadas de inserciones consideradas en clusters (centroides)
TEinsertionsCentroidEq.tibble=rmgr.tibble %>% filter(TEinsertion %in% famCentroidEq$TEinsertion)
TEinsertionsCentroidEq.v=famCentroidEq %>% pull(centroidID,TEinsertion)
TEinsertionsCentroidEq.tibble=TEinsertionsCentroidEq.tibble %>% mutate(centroidID=TEinsertionsCentroidEq.v[TEinsertion])

TEinsertionsCentroidEq.gr=makeGRangesFromDataFrame(TEinsertionsCentroidEq.tibble,keep.extra.columns=T)
names(TEinsertionsCentroidEq.gr)=TEinsertionsCentroidEq.gr$TEinsertion

seqnamesSizes=read.table(as.character(args[8]), fill=T, quote="", col.names=c("seqnames","width")) %>% tibble

transcriptsAnnotationPath=as.character(args[9])
if(transcriptsAnnotationPath!="none"){
	transcriptsAnnotation=import(transcriptsAnnotationPath)
}else{
	print("transcript annotation BED not provided. Ignoring")
}

TEsuperfamily=as.character(args[10])

outprefix=args[11]



print("loading centroids BED...")
centroids.gr=import(centroidsBed)
names(centroids.gr)=centroids.gr %>% mcols %>% as.data.frame %>% .$name %>% as.character

centroids.seqnames.v=centroids.gr %>% as.data.frame %>% as.tibble(rownames="centroidID") %>% pull(seqnames,centroidID)
centroids.start.v=centroids.gr %>% as.data.frame %>% as.tibble(rownames="centroidID") %>% pull(start,centroidID)
centroids.end.v=centroids.gr %>% as.data.frame %>% as.tibble(rownames="centroidID") %>% pull(end,centroidID)
centroids.strand.v=centroids.gr %>% as.data.frame %>% as.tibble(rownames="centroidID") %>% pull(strand,centroidID)


# fixing ORF coords
print("loading ORF BED...")
orfsBed=orfsBed %>% mutate(seqnames=centroids.seqnames.v[centroidID], start.centroid=centroids.start.v[centroidID],end.centroid=centroids.end.v[centroidID],strand.centroid=centroids.strand.v[centroidID]) %>% mutate(start.asm=start+start.centroid, end.asm=end+start.centroid)
orfs.start.v=orfsBed %>% pull(start.asm,ORF_ID)
orfs.end.v=orfsBed %>% pull(end.asm,ORF_ID)
orfs.strand.v=orfsBed %>% pull(strand,ORF_ID)



# fixing CDD domains coords
print("loading CDD results...")
rps=rps %>% mutate(start.orf=orfs.start.v[ORF_ID],end.orf=orfs.end.v[ORF_ID],strand.orf=orfs.strand.v[ORF_ID]) %>% mutate(domain.len=abs(to-from)) %>%
	mutate(start.domain.tmp=ifelse(strand.orf=="+",start.orf+from*3, end.orf-from*3)) %>% mutate(end.domain.tmp=ifelse(strand.orf=="+",
															   start.orf+to*3,
															   end.orf-to*3)) %>%
															   #start.domain.tmp+(domain.len*3),
															   #start.domain.tmp-(domain.len*3))) %>%
	mutate(strand.domain=ifelse(start.domain.tmp>end.domain.tmp,"-","+")) %>%
	mutate(start.domain=ifelse(start.domain.tmp>end.domain.tmp,end.domain.tmp,start.domain.tmp)) %>%
	mutate(end.domain=ifelse(start.domain.tmp>end.domain.tmp,start.domain.tmp,end.domain.tmp))
	#mutate(end.domain=end.domain*3)

rps$start.domain=as.integer(rps$start.domain)
rps$end.domain=as.integer(rps$end.domain)


# fixing BLAST hits coords
print("loading BLAST results...")
orfTocentroidID=orfsBed %>% pull(centroidID,ORF_ID)
bls$centroidID=orfTocentroidID[bls$ORF_ID]
bls$orfStrand=orfs.strand.v[bls$ORF_ID]

bls=bls %>% mutate(start.orf=orfs.start.v[ORF_ID],end.orf=orfs.end.v[ORF_ID]) %>%
	mutate(start.hit=start.orf+qstart*3,end.hit=start.orf+qend*3)

# fixing RT hits coords
print("loading HMM results...")
hmm$centroidID=orfTocentroidID[hmm$ORF_ID]
hmm$orfStrand=orfs.strand.v[hmm$ORF_ID]
hmm=hmm %>% mutate(start.orf=orfs.start.v[ORF_ID],end.orf=orfs.end.v[ORF_ID]) %>% mutate(start.hit=start.orf+aliFrom*3,end.hit=start.orf+aliTo*3)
# visualization window size
w.size=max(centroids.gr %>% width)+max(centroids.gr %>% width)/2



# Generamos un tibble para exportarlo a nivel de centroide
print("Generating centroid summary table...")
centroids.tibble=centroids.gr %>% as.data.frame %>% as.tibble(rownames="centroidID")

nORFs.v=orfsBed %>% count(centroidID) %>% pull(n,centroidID)
print("ORFs ok")

#nFrags.v=left_join(famCentroidEq,rmgr.tibble,by="TEinsertion") %>% count(centroidID) %>% pull(n,centroidID)
#print("Frags ok")

#nBPs.v=left_join(famCentroidEq,rmgr.tibble,by="TEinsertion") %>% group_by(centroidID) %>% summarise(nBPs = sum(width)) %>% pull(nBPs,centroidID)
#print("BPs ok")

fam.v=left_join(famCentroidEq,rmgr.tibble,by="TEinsertion") %>% pull(TEfam,centroidID)
print("fams ok")

print("loading seqnames sizes...")
genomeSize=seqnamesSizes %>% pull(width) %>% sum
superfamilyTotalContrib=rmgr.tibble %>% filter(TEsuperfam==TEsuperfamily) %>% pull(width) %>% sum


# txs
print("Searching for overlaps with assembled transcripts...")
overlappingTx.v=tibble(centroidID=centroids.tibble$centroidID,value="no-evaluated") %>% pull(value,centroidID)
if(transcriptsAnnotationPath!="none"){

	txs.overlaps=findOverlapPairs(TEinsertionsCentroidEq.gr,transcriptsAnnotation)
	l=txs.overlaps %>% length
	if(l > 0){
		print("There are overlaps :)")
		txs.gr=txs.overlaps@second
		txs.gr$TEinsertion=txs.overlaps@first %>% names
	}else{
		print("No overlaps found between assembled transcripts and centroid insertions.")
	}
	overlappingTx.v=txs.gr %>% as.data.frame %>% as_tibble() %>% pull(name,TEinsertion)
		

}

print("Starting to fill consolidated Table...")
centroids.tibble=centroids.tibble %>% mutate(	TEfam=fam.v[centroidID],
					     	n_ORFs=nORFs.v[centroidID],
						overlappingTx=ifelse(is.na(overlappingTx.v[centroidID])==T,F,T))


# diccionario para saber si centroide sobrelapa transcrito
centroidsOverlappingTx.v=centroids.tibble %>% pull(overlappingTx,centroidID)





print("working on Gviz...")



# palette
trackColours=c("#636363", pal_futurama("planetexpress",alpha = 0.6)(12))
trackColoursDark=c("#636363", pal_futurama("planetexpress")(12))
featureBorderColor="black"


#dicts for completing centroids.tibble
allCentroids=centroids.gr %>% mcols %>% as.data.frame %>% .$name %>% as.character


n_fragments_overlapping_txs.v=vector(mode="list",length(allCentroids))
names(n_fragments_overlapping_txs.v)=allCentroids

bps_overlapping_txs.v=vector(mode="list",length(allCentroids))
names(bps_overlapping_txs.v)=allCentroids

total_bps.v=vector(mode="list",length(allCentroids))
names(total_bps.v)=allCentroids

total_frags.v=vector(mode="list",length(allCentroids))
names(total_frags.v)=allCentroids

nRTs.v=rep(0,length(allCentroids))
names(nRTs.v)=allCentroids

if(transcriptsAnnotationPath!="none"){
	print("removing transcript redundancy on overlaps...")
	transcriptsAnnotation=GenomicRanges::reduce(transcriptsAnnotation)
}
# iterate through centroids for plotting
pdf(paste0(outprefix,".selected.centroids.pdf"),height=28,width=40)  #14x20 # antes era 7x14
for(c in allCentroids){

	#Centroids
	print(paste0("working on ",c))
	axisTrack=GenomeAxisTrack(col="black")
	c.gr=centroids.gr[c]
	print("centroid GRanges done...")



	# obtain fragments of this centroid (clustering based) and calculate totals
	allTEfrags.gr=famCentroidEq %>% filter(centroidID==c) %>% left_join(.,rmgr.tibble,by="TEinsertion") %>% makeGRangesFromDataFrame(keep.extra.columns=T)
	total_frags=allTEfrags.gr %>% length
	total_frags.v[c]=total_frags
	total_bps=GenomicRanges::reduce(allTEfrags.gr) %>% width %>% sum
	total_bps.v[c]=total_bps
	
	# calculating overlaps of all fragments of this centroid against transcripts
	if(transcriptsAnnotationPath!="none"){

		oPairs=findOverlapPairs(allTEfrags.gr, transcriptsAnnotation)
		n_frags_ovelapping_txs=oPairs@first$TEinsertion %>% unique %>% length
		n_fragments_overlapping_txs.v[c]=n_frags_ovelapping_txs

		bps_ovelapping_txs=GenomicRanges::intersect(GenomicRanges::reduce(allTEfrags.gr), transcriptsAnnotation) %>% width %>% sum
		if(length(bps_ovelapping_txs)==0){
			bps_ovelapping_txs=0
		}
		bps_overlapping_txs.v[c]=bps_ovelapping_txs

	}

	w.gr=resize(c.gr, w.size, fix="center", use.names=TRUE, ignore.strand=FALSE)
	print("window GRanges done...")

	centroidTrack=AnnotationTrack(start = start(c.gr),
				      end = end(c.gr),
				      strand = strand(c.gr),
				      id = as.character(c), 
				      fill= trackColours[1],
				      col=trackColoursDark[1],
				      name = as.character(c))

	print("Centroid track done...")

	#creamos una lista con los tracks para GViz
	trackListDefault=list()



	#Transcripts
	if(transcriptsAnnotationPath!="none"){
		if(any(centroidsOverlappingTx.v[c])){
			print(paste0("the centroid ",c," is overlapping assembled transcripts"))
			tx.gr=txs.gr[txs.gr$TEinsertion==c]
			txName=tx.gr$name
			txsTrack=AnnotationTrack(     start = start(tx.gr),
						 end = end(tx.gr),
						 strand = strand(tx.gr),
						 #chromosome = seqnames(tx.gr),
						 id = as.character(txName), 
						 fill= trackColours[6],
						 col=trackColoursDark[6],
						 name = as.character(paste0(c," Txs")))

			trackListDefault=append(trackListDefault,txsTrack)
			print("Transcripts track done...")
		}
	}else{
		print("No Transcripts track.")
	}

	#ORFs
	orfsInCentroid=orfsBed %>% filter(centroidID==as.character(c))
	if(NROW(orfsInCentroid)!=0){									# <- IF DE ORF
		orf.gr=GRanges(seqnames=orfsInCentroid$seqnames,
			       ranges=IRanges(start=orfsInCentroid$start.asm,end=orfsInCentroid$end.asm, name=orfsInCentroid$name),
			       strand=orfsInCentroid$strand,
			       group=orfsInCentroid$name)
		print("ORFs GRanges done...")


		# definimos los ORFs que corresponden al centroide actual
		orfsInCentroid.names=orfsInCentroid %>% pull(ORF_ID) %>% unique %>% as.character


		orfsTrack=AnnotationTrack( start=start(orf.gr),
					  end = end(orf.gr),
					  strand = strand(orf.gr),
					  id = as.character(names(orf.gr)), 
					  fill= trackColours[2],
					  col=trackColoursDark[2],
					  #group=strand(orf.gr),
					  name = as.character("ORFs"))

		trackListDefault=append(trackListDefault,list(centroidTrack,orfsTrack))
		print("ORFs track done...")



		# busqueda de dominios/hits/rt

		#Domains

		# filtramos dominios en funcion de los ORF en el centroide actual
		rpsInCentroid=rps %>% filter(ORF_ID %in% orfsInCentroid.names)
		if(NROW(rpsInCentroid)!=0){
			rps.gr=GRanges(seqnames=orfsInCentroid$seqnames %>% unique %>% as.character,
				       ranges=IRanges(start=rpsInCentroid$start.domain,end=rpsInCentroid$end.domain,name=rpsInCentroid$accession),strand=rpsInCentroid$strand.domain,
				       shortName=rpsInCentroid$short_name)
			print("Domains GRanges done...")


			cddTrack=AnnotationTrack( start=start(rps.gr),
						 end = end(rps.gr),
						 strand = strand(rps.gr),
						 fill= trackColours[3],
						 col= trackColoursDark[3],
						 id = as.character(rps.gr$shortName), 
						 name = as.character("Domains"))

			trackListDefault=append(trackListDefault,cddTrack)
			print("CDD track done...")
		}else{
			print("No CDD track.")
		}


		# también con los hits BLASTP
		blsInCentroid=bls %>% filter(ORF_ID %in% orfsInCentroid.names)

		if(NROW(blsInCentroid)!=0){
			bls.gr=GRanges(seqnames=orfsInCentroid$seqnames %>% unique %>% as.character,
				       ranges=IRanges(start=blsInCentroid$start.hit,end=blsInCentroid$end.hit,name=blsInCentroid$hitAccession),
				       strand=blsInCentroid$orfStrand,
				       shortName=blsInCentroid$hitAccession)
			print("BLAST hits GRanges done...")
			blsTrack=AnnotationTrack( start=start(bls.gr),
						 end = end(bls.gr),
						 strand = strand(bls.gr),
						 id = as.character(bls.gr$shortName), 
						 fill= trackColours[4],
						 col= trackColoursDark[4],
						 #				  stackHeight = 0.3,
						 name = as.character("BLAST hits"))
			trackListDefault=append(trackListDefault,blsTrack)
			print("BLAST track done...")
		}else{
			print("No BLAST track.")
		}
		# y con los hits HMMER contra RT	

		hmmInCentroid=hmm %>% filter(ORF_ID %in% orfsInCentroid.names)
		nRTs.v[c]=hmmInCentroid %>% pull(ORF_ID) %>% length

		if(NROW(hmmInCentroid)!=0){
			hmm.gr=GRanges(seqnames=orfsInCentroid$seqnames %>% unique %>% as.character,
				       ranges=IRanges(start=hmmInCentroid$start.hit,end=hmmInCentroid$end.hit,name=hmmInCentroid$RT_ID),
				       strand=hmmInCentroid$orfStrand,
				       shortName=hmmInCentroid$RT_ID)
			print("HMMER RT hits GRanges done...")
			hmmTrack=AnnotationTrack( start=start(hmm.gr),
						 end = end(hmm.gr),
						 strand = strand(hmm.gr),
						 id = as.character(hmmInCentroid$RT_ID), 
						 fill= trackColours[5],
						 col= trackColoursDark[5],
						 name = as.character("RT domain"))
			trackListDefault=append(hmmTrack,trackListDefault)
			print("RT track done...")
		}else{
			print("F**k!, no HMMER RT hits!")
		}



	# ELSE ORF
	}else{
		trackListDefault=append(trackListDefault,centroidTrack)
		print("No ORFs track.")
	}


	print("I'll start plotting...")

	# The real plotting
	plotTracks(append(trackListDefault,axisTrack),to=start(w.gr),from=end(w.gr),scale=0.2, main=as.character(c), cex.main=1.5) #,featureAnnotation = "id")

}


print("Last details while consolidating Table...")
total_frags.v=total_frags.v %>% unlist
total_bps.v=total_bps.v %>% unlist

if(transcriptsAnnotationPath!="none"){
	n_fragments_overlapping_txs.v=n_fragments_overlapping_txs.v %>% unlist
	bps_overlapping_txs.v=bps_overlapping_txs.v %>% unlist

	centroids.tibble=centroids.tibble %>% mutate(n_RT_domains=nRTs.v[centroidID],
						     n_fragments=total_frags.v[centroidID],
						     n_fragments_overlapping_txs=n_fragments_overlapping_txs.v[centroidID],
						     perc_fragments_overlapping_txs=n_fragments_overlapping_txs.v[centroidID]*100/total_frags.v[centroidID],	
						     n_bps=total_bps.v[centroidID],
						     n_bps_overlapping_txs=bps_overlapping_txs.v[centroidID],
						     perc_bps_overlapping_txs=bps_overlapping_txs.v[centroidID]*100/total_bps.v[centroidID])
}else{
	centroids.tibble=centroids.tibble %>% mutate(n_RT_domains=nRTs.v[centroidID],
						     n_fragments=total_frags.v[centroidID],
						     n_fragments_overlapping_txs="no-data",
						     perc_fragments_overlapping_txs="no-data",	
						     n_bps=total_bps.v[centroidID],
						     n_bps_overlapping_txs="no-data",
						     perc_bps_overlapping_txs="no-data")
}

dev.off()
save(centroids.tibble,file=paste0(outprefix,".selected.centroids.table.RData"))
save.image(paste0(outprefix,".selected.centroids.RData"))
