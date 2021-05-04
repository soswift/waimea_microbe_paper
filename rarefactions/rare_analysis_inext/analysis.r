#!/usr/bin/env Rscript

# nohup ./rare_analysis_mm_fulldata.r >status.txt &
# cd ~/Documents/coauthor/anthony_waimea/rare_fig/rare_analysis_SAR/

# greet
	message("it's started!")

# load stuff
	library(iNEXT)
	library(data.table)
	library(parallel)

# read in data
	otutable_fp <- "../../data/raw/sequence/waimea_16S_v1/abundance_table_100.shared"
	metadata_fp <- "../../data/raw/sequence/waimea_16S_v1/combined_prep_sample_metadata_16s.csv"
	taxonomy_fp <- "../../data/raw/sequence/waimea_16S_v1/annotations_100.taxonomy"
	otutable <- fread(otutable_fp, sep="\t")
	metadata <- read.delim(metadata_fp, sep=",", header=T, stringsAsFactors=FALSE)
	taxonomy <- read.delim(taxonomy_fp, sep="\t", header=T, stringsAsFactors=FALSE)

# convert otutable to matrix
	otutable <- as.matrix(otutable)
	rownames(otutable) <- otutable[,colnames(otutable)=="Group"]
	# remove first 3 columns
	otutable <- otutable[,4:ncol(otutable)]

# get rid of stuff
	metadata <- metadata[metadata$habitat %in% c("Marine", "Riverine", "Terrestrial"), ]
	# remove samples from otutable not in metadata
	otutable <- otutable[rownames(otutable) %in% metadata$sequencing_id,]
	# remove otus that are not present in any samples
	otutable <- otutable[,colSums(otutable) > 0]
	# get rid of useless OTUs
	chlorOTUs <- taxonomy$OTU[grepl(x=taxonomy$Taxonomy, pattern="chloroplast", ignore.case=T)]
	mitoOTUs <- taxonomy$OTU[grepl(x=taxonomy$Taxonomy, pattern="mitochond", ignore.case=T)]
	archaOTUs <- taxonomy$OTU[grepl(x=taxonomy$Taxonomy, pattern="archaea", ignore.case=T)]
	nonbacOTUs <- taxonomy$OTU[!grepl(taxonomy$Taxonomy, pattern="Bacteria;")]
	otutable <- otutable[, ! colnames(otutable) %in% c(chlorOTUs, mitoOTUs, archaOTUs, nonbacOTUs)]

# transpose otutable to "samples on top" format
	otutable <- t(otutable)

# align otutable and metadata
	otutable <- otutable[,colnames(otutable) %in% metadata$sequencing_id]
	metadata <- metadata[metadata$sequencing_id %in% colnames(otutable), ]
	otutable <- otutable[,order(colnames(otutable))]
	metadata <- metadata[order(metadata$sequencing_id),]
	message("Next line better be TRUE (checking for metadata/otutable sort):")
	all(colnames(otutable) == metadata$sequencing_id)

# collapse "Omnivore", "Herbivore", "Detritivore", "Carnivore" into "Consumer"
	consumers <- c("Omnivore", "Herbivore", "Detritivore", "Carnivore")
	metadata$trophic[metadata$trophic %in% consumers] <- "Consumer"

# status update
	message("Data load and organization complete")




# function to get vector from otutable given hab+tro
	otb_habtro <- function(hab, tro, otb=otutable, md=metadata){
		if(hab %in% c("All", "all")){ hab <- c("Marine", "Riverine", "Terrestrial") }
		if(tro %in% c("All", "all")){ tro <- c("Environmental", "PrimaryProducer", "Consumer") }
		samps2get <- md$sequencing_id[ md$habitat%in%hab & md$trophic%in%tro ]
		otb[,colnames(otb) %in% samps2get]
	}

# function to rowsum otutable then prune rel abund vector by 
	# removing zeroes
	rs_otb <- function(otb){
		out <- rowSums(otb)
		out <- out[out>0]
		return(out)
	}

# iNEXT helper functions
	xyrange_inext <- function(x){
		xrange <- range(x$iNextEst[[1]]$m)
		yrange <- range(x$iNextEst[[1]]$qD)
		return(cbind(xrange, yrange))
	}
	addcurve_inext <- function(x, lcol="blue", lpch=15){
		xdf <- x$iNextEst[[1]]
		xvals_intrp <- xdf$m[xdf$method=="interpolated"]
		yvals_intrp <- xdf$qD[xdf$method=="interpolated"]
		xvals_extrp <- xdf$m[xdf$method=="extrapolated"]
		yvals_extrp <- xdf$qD[xdf$method=="extrapolated"]
		xobs <- xdf$m[xdf$method=="observed"]
		yobs <- xdf$qD[xdf$method=="observed"]
		points(xvals_intrp, yvals_intrp, type="l", col=lcol, lty=1)
		points(xvals_extrp, yvals_extrp, type="l", col=lcol, lty=2)
		points(xobs, yobs, pch=lpch, col=lcol, cex=2)
	}

	inextwrap <- function(hab, tro ){
		iNEXT(list(rs_otb(otb_habtro(hab, tro))), q=0, datatype="abundance")
	}

	plot_troslist <- function(tros, cols, pchs){
		ranges <- do.call("rbind", lapply(X=tros, FUN=xyrange_inext))
		plot(0, type="n", xlim=range(ranges[,1]), ylim=range(ranges[,2]),
			xlab="Sampling effort (# of sequences)", ylab="Richness")
		for(i in 1:length(tros)){
			addcurve_inext(tros[[i]], cols[i], pchs[i])
		}
	}
	# status update
	message("Functions defined")


	# tros            trocols 
	# Environmental   #fc8d62 
	# PrimaryProducer #66c2a5 
	# Consumer        #756bb1 

	# habs            habcols
	# Marine          #1f78b4
	# Riverine        #a6cee3
	# Terrestrial     #a6611a



	habs <- c("Marine", "Riverine", "Terrestrial", "All")
	habcols <- c("#1f78b4", "#a6cee3", "#a6611a", "black")

	tros <- c("Environmental", "PrimaryProducer", "Consumer", "All")
	trocols <- c("#fc8d62", "#66c2a5", "#756bb1", "black")
	

	plotdata_byTro <- mclapply(
		X=tros,
		FUN=function(tro){inextwrap("All", tro)},
		mc.cores=4
	)
	saveRDS(plotdata_byTro, file="plotdata_byTro_All.RDS")


	plotdata_byHab <- mclapply(
		X=habs,
		FUN=function(hab){inextwrap(hab, "All")},
		mc.cores=4
	)
	saveRDS(plotdata_byHab, file="plotdata_byHab_All.RDS")

	


# plot
	# By Trophic
	pdf("rare_curves_full_trophic.pdf")
	plot_troslist(plotdata_byTro, cols=trocols, pchs=15:18)
	dev.off()


	# Percentage of estimated diversity
	tros_percentages <- data.frame(
		Trophic=tros,
		ObservedRichness=sapply(X=plotdata_byTro, FUN=function(x){x$AsyEst$Observed[1]}),
		ChaoEstRichness=sapply(X=plotdata_byTro, FUN=function(x){x$AsyEst$Estimator[1]}),
		PercentEst=sapply( X=plotdata_byTro, FUN=function(x){ (x$AsyEst$Observed[1] / x$AsyEst$Estimator[1]) * 100})
	)
	write.table(tros_percentages, file="tros_percentages.csv", sep='\t', row.names=FALSE, quote=FALSE)
	
	# By Habitat
	pdf("rare_curves_full_habitat.pdf")
	plot_troslist(plotdata_byHab, cols=habcols, pchs=15:18)
	dev.off()


	# Percentage of estimated diversity
	habs_percentages <- data.frame(
		Habitat=habs,
		ObservedRichness=sapply(X=plotdata_byHab, FUN=function(x){x$AsyEst$Observed[1]}),
		ChaoEstRichness=sapply(X=plotdata_byHab, FUN=function(x){x$AsyEst$Estimator[1]}),
		PercentEst=sapply( X=plotdata_byHab, FUN=function(x){ (x$AsyEst$Observed[1] / x$AsyEst$Estimator[1]) * 100})
	)
	write.table(habs_percentages, file="habs_percentages.csv", sep='\t', row.names=FALSE, quote=FALSE)
	
