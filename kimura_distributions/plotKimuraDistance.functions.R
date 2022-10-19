buildTEpaletteKimura=function(nonRedundantTElist){
	library(RColorBrewer)
	nTEs=length(nonRedundantTElist)
	pal=colorRampPalette(brewer.pal(9, "Set1"))(nTEs)
	names(pal)=nonRedundantTElist
	return(pal)
}


plotKimuraCompareBetweenSpecies=function(d.bed.lst, TEclassification=NULL,classToShow="LINE",pal=NULL){
	for(sp in names(d.bed.lst)){
		d.bed.lst[[sp]]$Specie=as.character(sp)
	}
	d=do.call("rbind",d.bed.lst)
	if(!is.null(TEclassification)){
		TEorder.v=TEclassification %>% pull(TEorder,id)
		TEsuperfam.v=TEclassification %>% pull(TEsuperfam,id)
		d$TEorder=TEorder.v[d$name]
		d$TEsuperfam=TEsuperfam.v[d$name]
	}
	return(d)
}





plotKimura=function(d.bed,TEfam.lst=NULL,title="Kimura distance distribution",xmax=50,ymax=NULL,TEorders2plot=unique(d.bed$TEorder),palette=NULL){
	`%notin%` <<- Negate(`%in%`)
	TEorders2plot=TEorders2plot[TEorders2plot %notin% c("Unknown","")]



	if(!is.null(TEfam.lst)){
		if(is.null(palette)){
			palette=colorRampPalette(brewer.pal(9, "Set1"))(length(TEfam.lst))

			# for automatic palette
			if(is.null(names(palette))){
				names(palette)=TEfam.lst
			}
		}

		exclude.lst=TEorders2plot
		TEfam.lst=TEfam.lst[TEfam.lst %notin% exclude.lst]


		# dejamos fuera las TEsuperfams que se llaman igual que el orden, xq genera conflicto
		palette=palette[names(palette) %notin% TEorders2plot]

		orderColor="#666666"
		orderPalette=rep(orderColor,length(TEorders2plot))
		names(orderPalette)=TEorders2plot
		palette=c(palette,orderPalette)
	}



	dict=vector(mode="list",length=length(TEorders2plot))
	names(dict)=TEorders2plot

	if(!is.null(TEfam.lst)){
		TEfams=d.bed[d.bed$TEsuperfam %in% TEfam.lst,]
	}
	else{
		TEfams=d.bed
	}
	for(Order in TEorders2plot){
		presence.boolean.order = as.character(Order) == d.bed$TEorder
		if(any(presence.boolean.order)){
			print(as.character(Order))
			TEfams.order=TEfams[TEfams$TEorder == as.character(Order),]

			# global
			kimura.freq.order=table(as.integer(d.bed[d.bed$TEorder == as.character(Order),]$kimura))

			order.df=data.frame(kimura.freq.order,as.character(Order))
			colnames(order.df)=c("Kimura","Freq","Label")
			dict[[Order]]=rbind(dict[[Order]],order.df)

			if(!is.null(TEfam.lst)){
				for(TEfam in TEfam.lst){
					presence.boolean = as.character(TEfam) == TEfams.order$TEsuperfam

					if(any(presence.boolean)){
						TEfam.bed=TEfams.order[TEfams.order$TEsuperfam == as.character(TEfam),]

						kimura.freq.tefam=table(as.integer(TEfam.bed$kimura))
						TEfam.df=data.frame(kimura.freq.tefam,as.character(TEfam))
						colnames(TEfam.df)=c("Kimura","Freq","Label")


						dict[[Order]]=rbind(dict[[Order]],TEfam.df)

					}
					else{
						warning(paste(as.character(TEfam)," was excluded from the plot since it wasn't found in the dataset (d.bed)"))
					}
				}# cierre del for de TEfams

				dict[[Order]]$Kimura=as.numeric(dict[[Order]]$Kimura)
				dict[[Order]]$Freq=as.numeric(dict[[Order]]$Freq)

			}else{
				warning(paste("Order (",as.character(Order),") was excluded from the plot since it wasn't found in the dataset (d.bed)",sep=""))
				dict[Order]=NULL
			}

		} # cierre if TEfam.lst no es nulo
	}# cierre del for de ordenes



	buildGridPlot(title=title,dict,xmax,ymax,pal=palette)
}# cierre de la funcion






buildGridPlot=function(title,dict,xmax=50,ymax=NULL,pal){
        library(ggplot2)
        library(grid)
        library(gridExtra)


	if(is.null(pal)){
		pal=colorRampPalette(brewer.pal(9, "Set1"))(length(TEfam.lst))
	}

	p=vector(mode="list",length=length(names(dict)))
	names(p)=names(dict)

	for(Order in names(p)){
		if(is.numeric(ymax)){
        		p[[Order]]=ggplot(dict[[Order]],aes(x=Kimura,y=Freq,group=Label))+geom_line(data=dict[[Order]],aes(colour=Label),size=0.8)+ylim(0,ymax)+xlim(0,xmax) + ggtitle(as.character(Order)) + scale_color_manual(values=pal) + theme(legend.position = c(0.8, 0.8))
		}
		else{
		
        		p[[Order]]=ggplot(dict[[Order]],aes(x=Kimura,y=Freq,group=Label))+geom_line(data=dict[[Order]],aes(colour=Label),size=0.8)+xlim(0,xmax) + ggtitle(as.character(Order))+ scale_color_manual(values=pal) + theme(legend.position = c(0.8, 0.8))
		}
	}

	title.grob=textGrob(as.character(title),gp=gpar(fontsize=20,font=3))
	grid=do.call("arrangeGrob",c(p, ncol=2))
	grid.arrange(grid,top=title.grob)
	
	#grid   top=title.grob))
	#do.call("grid.arrange",c(p, ncol=2, top=textGrob(as.character(title),gp=gpar(fontsize=20,font=3))))

#        grid.arrange(line.p,sine.p,ltr.p,dna.p, nrow = 2, ncol=2,  top = textGrob(as.character(title),gp=gpar(fontsize=20,font=3)))
#        grid.arrange(p,  top = textGrob(as.character(title),gp=gpar(fontsize=20,font=3)))

}






# funcion modificada para plotear familias especificas



plotKimuraTEfam=function(d.bed,TEfam.lst=NULL,title="Kimura distance distribution",xmax=50,ymax=NULL,TEorders2plot=unique(d.bed$TEorder),palette=NULL){
	`%notin%` <<- Negate(`%in%`)
	TEorders2plot=TEorders2plot[TEorders2plot %notin% c("Unknown","")]



	if(!is.null(TEfam.lst)){
		if(is.null(palette)){
			palette=colorRampPalette(brewer.pal(9, "Set1"))(length(TEfam.lst))

			# for automatic palette
			if(is.null(names(palette))){
				names(palette)=TEfam.lst
			}
		}

		# dejamos fuera las TEsuperfams que se llaman igual que el orden, xq genera conflicto
		exclude.lst=TEorders2plot
		TEfam.lst=TEfam.lst[TEfam.lst %notin% exclude.lst]

		palette=palette[names(palette) %notin% TEorders2plot]

		orderColor="#666666"
		orderPalette=rep(orderColor,length(TEorders2plot))
		names(orderPalette)=TEorders2plot
		palette=c(palette,orderPalette)
	}



	dict=vector(mode="list",length=length(TEorders2plot))
	names(dict)=TEorders2plot

	# filtra el data.frame para quedarse con las TEfams de interes
	if(!is.null(TEfam.lst)){
		TEfams.df=d.bed[d.bed$TEsuperfam %in% TEfam.lst,]
	}
	else{
		TEfams.df=d.bed
	}



	for(Order in TEorders2plot){
		presence.boolean.order = as.character(Order) == d.bed$TEorder
		if(any(presence.boolean.order)){
			print(as.character(Order))
			TEfams.order=TEfams[TEfams$TEorder == as.character(Order),]

			# global
			kimura.freq.order=table(as.integer(d.bed[d.bed$TEorder == as.character(Order),]$kimura))

			order.df=data.frame(kimura.freq.order,as.character(Order))
			colnames(order.df)=c("Kimura","Freq","Label")
			dict[[Order]]=rbind(dict[[Order]],order.df)

			if(!is.null(TEfam.lst)){
				for(TEfam in TEfam.lst){
					presence.boolean = as.character(TEfam) == TEfams.order$TEsuperfam

					if(any(presence.boolean)){
						TEfam.bed=TEfams.order[TEfams.order$TEsuperfam == as.character(TEfam),]

						kimura.freq.tefam=table(as.integer(TEfam.bed$kimura))
						TEfam.df=data.frame(kimura.freq.tefam,as.character(TEfam))
						colnames(TEfam.df)=c("Kimura","Freq","Label")


						dict[[Order]]=rbind(dict[[Order]],TEfam.df)

					}
					else{
						warning(paste(as.character(TEfam)," was excluded from the plot since it wasn't found in the dataset (d.bed)"))
					}
				}# cierre del for de TEfams

				dict[[Order]]$Kimura=as.numeric(dict[[Order]]$Kimura)
				dict[[Order]]$Freq=as.numeric(dict[[Order]]$Freq)

			}else{
				warning(paste("Order (",as.character(Order),") was excluded from the plot since it wasn't found in the dataset (d.bed)",sep=""))
				dict[Order]=NULL
			}

		} # cierre if TEfam.lst no es nulo
	}# cierre del for de ordenes



	buildGridPlot(title=title,dict,xmax,ymax,pal=palette)
}# cierre de la funcion









plotTEsuperfamilyKimura=function(TEsuperfamily,rm.gr,size=3,shape=23,palette="Paired"){
	# dataframe filtrado por TEsuperfam
	rm.gr.filtered=rm.gr %>% filter(TEsuperfam == TEsuperfamily)

	# conteos
	TEfam.count=rm.gr.filtered %>% count(TEfam,Kimura)
#	TEsuperfam.count=rm.gr.filtered %>% count(TEsuperfam,Kimura)

	colnames(TEfam.count)=c("label","Kimura","Freq")
#	colnames(TEsuperfam.count)=c("label","Kimura","Freq")

#	d=rbind(TEfam.count,TEsuperfam.count)
	d=TEfam.count
	d$Kimura=as.numeric(d$Kimura)
	library(RColorBrewer)


	colourCount = length(unique(d$label))
	getPalette = colorRampPalette(brewer.pal(8, palette))


	p=ggplot(d,aes(x=Kimura,y=Freq,colour=label))+geom_point(size=size,shape=shape)+scale_color_manual(values = getPalette(colourCount))+theme_minimal()

		#scale_color_brewer(palette=palette)
	return(p)

}







