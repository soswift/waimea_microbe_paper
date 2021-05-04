##!/usr/bin/env Rscript

# greet
	message("it's started!")

# load stuff
	library("data.table")
	library("parallel")
	library(Rcpp)
	source("rare_funs_endos.r")

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

# how many samples of each category are there?
	samp_tab <- list()
	i <- 0
	for(hab in unique(metadata$habitat)){
		for(tro in unique(metadata$trophic)){
			n <- sum((metadata$habitat == hab) & (metadata$trophic==tro))
			i <- i + 1
			samp_tab[[i]] <- data.frame(
				habitat=hab,
				trophic=tro,
				nsamp=n
			)
		}
	}
	samp_tab <- do.call("rbind", samp_tab)
	write.table(samp_tab, file="samp_tab.csv", sep=",", row.names=F, quote=F)

# function to get matching samples from otutable given hab+tro
	otb_habtro <- function(hab, tro, otb=otutable, md=metadata){
		if(hab %in% c("All", "all")){ hab <- c("Marine", "Riverine", "Terrestrial") }
		if(tro %in% c("All", "all")){ tro <- c("Environmental", "PrimaryProducer", "Consumer") }
		samps2get <- md$sequencing_id[ md$habitat%in%hab & md$trophic%in%tro ]
		otb[,colnames(otb) %in% samps2get]
	}

# vars
	maxdepth <- 70000
	bindepth <- 60000
	stepsize <- 300


# split by habitat
	remake <- TRUE
	colcurves_fp <- "colcurves_tros.rdata"
	tros <- c("Environmental", "PrimaryProducer", "Consumer")

	if(remake){
		colcurves_list <- list()
		for(i in 1:length(tros)){
			print(paste("doing", tros[i]))
			otb <- otb_habtro(hab="All", tro=tros[i])
			otb <- otb[rowSums(otb) > 0, ]
			otb <- otb[, colSums(otb) > bindepth ]
			colcurves_list[[i]] <- multiple_collectors_curves_parallel(
				otb, maxd=maxdepth, stepsize=stepsize, ncores=10)
		}
		save(list="colcurves_list", file=colcurves_fp)
	}else{
		load(colcurves_fp)
	}

	pdf("tros_plot.pdf")
	par(mfrow=c(2,2))
	for(i in 1:length(tros)){
		plot_multiple_binned_collectors_curves(colcurves_list[[i]], n_bins=5,
			bin_d=20000, main=tros[i], ylab="ASV Richness")
	}
	dev.off()



# split by habitat
	remake <- TRUE
	colcurves_fp <- "colcurves_habs.rdata"
	habs <- c("Marine", "Riverine", "Terrestrial")

	if(remake){
		colcurves_list <- list()
		for(i in 1:length(habs)){
			print(paste("doing", habs[i]))
			otb <- otb_habtro(hab=habs[i], tro="All")
			otb <- otb[rowSums(otb) > 0, ]
			otb <- otb[, colSums(otb) > bindepth ]
			colcurves_list[[i]] <- multiple_collectors_curves_parallel(
				otb, maxd=maxdepth, stepsize=stepsize, ncores=10)
		}
		save(list="colcurves_list", file=colcurves_fp)
	}else{
		load(colcurves_fp)
	}

	pdf("habs_plot.pdf")
	par(mfrow=c(2,2))
	for(i in 1:length(habs)){
		plot_multiple_binned_collectors_curves(colcurves_list[[i]], n_bins=5,
			bin_d=40000, main=habs[i], ylab="ASV Richness")
	}
	dev.off()



