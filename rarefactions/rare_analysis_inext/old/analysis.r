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
	otutable_fp <- "../../rawdata/waimea_full/abundance_table_100.shared"
	metadata_fp <- "../../rawdata/waimea_full/combined_prep_sample_metadata_16s.csv"
	taxonomy_fp <- "../../rawdata/waimea_full/annotations_100.taxonomy"
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

	plot_troslist <- function(tros, cols, pchs, ...){
		ranges <- do.call("rbind", lapply(X=tros, FUN=xyrange_inext))
		plot(0, type="n", xlim=range(ranges[,1]), ylim=range(ranges[,2]),
			xlab="Sampling effort (# of sequences)", ylab="Richness", ...)
		for(i in 1:length(tros)){
			addcurve_inext(tros[[i]], cols[i], pchs[i])
		}
	}
	# status update
	message("Functions defined")



	habs <- c("Marine", "Riverine", "Terrestrial")
	tros <- c("Environmental", "PrimaryProducer", "Consumer")
	trocols <- c("#fc8d62", "#66c2a5", "#756bb1")

	plotdata <- list()

	for(i in 1:length(habs)){
		plotdata[[i]] <- mclapply(
			X=tros,
			FUN=function(x){
				out <- inextwrap(habs[i], x)
				message(paste0("...", habs[i], "/", x, " done"))
				return(out)
			},
			mc.cores=3
		)
		names(plotdata[[i]]) <- habs[i]
	}

	save(list="plotdata", file="plotdata.rdata")

	# for each hab, use plot_troslist() to make a plot. do them in a 2x2
	pdf("inext_curves.pdf")
	par(mfrow=c(2,2))
	for(i in 1:length(habs)){
		plot_troslist(plotdata[[i]], cols=trocols, pchs=c(15:17), main=habs[i])
	}
	plot(x=c(1,1,1), y=c(1,2,3), pch=15:17, col=trocols, ylim=c(0,4), xlim=c(0,4), cex=3,
		xaxt="n", yaxt="n", bty="n", xlab="", ylab="")
	text(x=c(2,2,2), y=c(1,2,3), labels=tros, adj=c(0, 0.5))
	dev.off()



