##!/usr/bin/env Rscript

# nohup ./rare_analysis_mm_fulldata.r >status.txt &
# cd ~/Documents/coauthor/anthony_waimea/rare_fig/rare_analysis_SAR/

# greet
	message("it's started!")

# load stuff
	library("data.table")
	library("parallel")
	library("drc")
	library("shape")
	source("rare_analysis_sar_funs.r")

# read in data
	otutable_fp <- "../../../data/raw/sequence/waimea_16S_v1/abundance_table_100.shared"
	metadata_fp <- "../../../data/raw/sequence/waimea_16S_v1/combined_prep_sample_metadata_16s.csv"
	taxonomy_fp <- "../../../data/raw/sequence/waimea_16S_v1/annotations_100.taxonomy"
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

# function to make color semitransparent
	col2alpha <- function(col, al01=0.5){
		c2rgb <- col2rgb(col)
		al <- al01 * 255
		rgb(red=c2rgb[1,1], green=c2rgb[2,1], blue=c2rgb[3,1], alpha=al, maxColorValue=255)
	}


# curves by tro, all habs
	remake <- FALSE
	tros <- c("Environmental", "PrimaryProducer", "Consumer", "All")
	cols <-  c("#fc8d62", "#66c2a5", "#756bb1", "black")

	pchs <- 15:18

	tablelist <- list()

	if(remake){
		tros_rbrtabs <- list()
	}else{
		load("tros_rbrtabs.rdata")
	}
	pdf("tros_plots.pdf")
	par(mfrow=c(2,2))
	for(i in 1:length(tros)){
		message(paste0("Starting analysis of ", tros[i]))

		if(remake){
			otb <- otb_habtro("All", tros[i])
			otb <- otb[rowSums(otb) > 0, ]
			set.seed(12345)
			colectcurves2 <- mclapply(X=1:50, FUN=function(x){
				data.frame(
					rich=collect_samps_otb(otb, 0),
					nsam=1:ncol(otb)
				)
			}, mc.cores=20)
			tros_rbrtabs[[i]] <- do.call("rbind", colectcurves2)
		}

		mm_mod <- drm(rich ~ nsam, data=tros_rbrtabs[[i]], fct=MM.3())

		ee <- coef(mm_mod)[3]
		dd <- coef(mm_mod)[2] # expected total richness
		cc <- coef(mm_mod)[1]
		dd_lci <- confint(mm_mod)[2,1]
		dd_hci <- confint(mm_mod)[2,2]


		exp_nsamp_90percent <- (ee*(0.90 * dd - cc)) / (dd - 0.90 * dd)

		interpolated <- data.frame( nsam = unique(tros_rbrtabs[[i]]$nsam) )
		interpolated$rich <- predict(mm_mod, interpolated)

		extrapolated <- data.frame( nsam = seq(from=max(tros_rbrtabs[[i]]$nsam), to=exp_nsamp_90percent, length.out=100) )
		extrapolated$rich <- predict(mm_mod, extrapolated)

		obs_total_rich <- max(tros_rbrtabs[[i]]$rich)
		percent_est <- round((obs_total_rich / dd) * 100, 2)
		percent_lci <- round((obs_total_rich / dd_hci) * 100, 2)
		percent_hci <- round((obs_total_rich / dd_lci) * 100, 2)

		mm3 <- function(x, cc, dd, ee){ cc + (dd - cc) / (1 + ee/x)	}
		extrap_d_95ci <- data.frame(
			nsam=extrapolated$nsam,
			lci=mm3(extrapolated$nsam, cc, dd_lci, ee),
			hci=mm3(extrapolated$nsam, cc, dd_hci, ee)
		)

		plot(0, type="n", 
			xlim=c(1, max(extrapolated$nsam)), 
			ylim=c(min(tros_rbrtabs[[i]]$rich), max(extrapolated$rich)),
			xlab="Number of Samples", ylab="Richness", main=tros[i]
		)

		poly_y <- c(extrap_d_95ci$hci, rev(extrap_d_95ci$lci))
		poly_x <- c(extrap_d_95ci$nsam, rev(extrap_d_95ci$nsam))
		polygon(x=poly_x, y=poly_y, border=NA, col=col2alpha(cols[i], 0.25))

		points(x=tros_rbrtabs[[i]]$nsam, y=tros_rbrtabs[[i]]$rich, pch=20, cex=0.1, col=cols[i])

		# points(x=interpolated$nsam, y=interpolated$rich, col=cols[i], lwd=3, type="l")
		points(x=extrapolated$nsam, y=extrapolated$rich, col=cols[i], lwd=3, type="l", lty=2)
		points(x=max(interpolated$nsam), y=max(interpolated$rich), col=cols[i], pch=pchs[i], cex=3)

		txt_x <- max(extrapolated$nsam) / 2
		ymin <- min(interpolated$rich)
		yrange <- max(extrapolated$rich) - min(interpolated$rich)
		txt_y1 <- ymin + yrange * (1/8)
		txt_y2 <- ymin + yrange * (2/8)
		txt_y3 <- ymin + yrange * (3/8)
		fcex <- 1.3
		text(txt_x, txt_y3, paste0("HCI: ", percent_hci, "%"), 
			adj=c(0, 0.5), col=cols[i], family="mono", cex=fcex)
		text(txt_x, txt_y2, paste0("EST: ", percent_est, "%"),
			adj=c(0, 0.5), col=cols[i], family="mono", cex=fcex, font=2)
		text(txt_x, txt_y1, paste0("LCI: ", percent_lci, "%"),
			adj=c(0, 0.5), col=cols[i], family="mono", cex=fcex)

		tablelist[[i]] <- data.frame(x=tros[i], total_observed_rich=max(tros_rbrtabs[[i]]$rich), mm_rich=dd, nsamp=max(tros_rbrtabs[[i]]$nsam))
	}
	dev.off()

	if(remake){
		save(list="tros_rbrtabs", file="tros_rbrtabs.rdata")
	}
	trostab <- do.call("rbind", tablelist)
	write.table(file="tros_table.txt", trostab, sep="\t", row.names=FALSE, quote=FALSE)





# curves by hab, all tros
	remake <- FALSE
	tablelist <- list()
	tros <- c("Marine", "Riverine", "Terrestrial", "All")
	cols <- c("#1f78b4", "#a6cee3", "#fed9a6", "black")


	pchs <- 11:14 #15:18

	if(remake){
		tros_rbrtabs <- list()
	}else{
		load("habs_rbrtabs.rdata")
	}
	pdf("habs_plots.pdf")
	par(mfrow=c(2,2))
	for(i in 1:length(tros)){
		message(paste0("Starting analysis of ", tros[i]))

		if(remake){
			otb <- otb_habtro(tros[i], "All")
			otb <- otb[rowSums(otb) > 0, ]
			set.seed(12345)
			colectcurves2 <- mclapply(X=1:50, FUN=function(x){
				data.frame(
					rich=collect_samps_otb(otb, 0),
					nsam=1:ncol(otb)
				)
			}, mc.cores=20)
			tros_rbrtabs[[i]] <- do.call("rbind", colectcurves2)
		}

		mm_mod <- drm(rich ~ nsam, data=tros_rbrtabs[[i]], fct=MM.3())

		ee <- coef(mm_mod)[3]
		dd <- coef(mm_mod)[2] # expected total richness
		cc <- coef(mm_mod)[1]
		dd_lci <- confint(mm_mod)[2,1]
		dd_hci <- confint(mm_mod)[2,2]


		exp_nsamp_90percent <- (ee*(0.90 * dd - cc)) / (dd - 0.90 * dd)

		interpolated <- data.frame( nsam = unique(tros_rbrtabs[[i]]$nsam) )
		interpolated$rich <- predict(mm_mod, interpolated)

		extrapolated <- data.frame( nsam = seq(from=max(tros_rbrtabs[[i]]$nsam), to=exp_nsamp_90percent, length.out=100) )
		extrapolated$rich <- predict(mm_mod, extrapolated)

		obs_total_rich <- max(tros_rbrtabs[[i]]$rich)
		percent_est <- round((obs_total_rich / dd) * 100, 2)
		percent_lci <- round((obs_total_rich / dd_hci) * 100, 2)
		percent_hci <- round((obs_total_rich / dd_lci) * 100, 2)

		mm3 <- function(x, cc, dd, ee){ cc + (dd - cc) / (1 + ee/x)	}
		extrap_d_95ci <- data.frame(
			nsam=extrapolated$nsam,
			lci=mm3(extrapolated$nsam, cc, dd_lci, ee),
			hci=mm3(extrapolated$nsam, cc, dd_hci, ee)
		)

		plot(0, type="n", 
			xlim=c(1, max(extrapolated$nsam)), 
			ylim=c(min(tros_rbrtabs[[i]]$rich), max(extrapolated$rich)),
			xlab="Number of Samples", ylab="Richness", main=tros[i]
		)

		poly_y <- c(extrap_d_95ci$hci, rev(extrap_d_95ci$lci))
		poly_x <- c(extrap_d_95ci$nsam, rev(extrap_d_95ci$nsam))
		polygon(x=poly_x, y=poly_y, border=NA, col=col2alpha(cols[i], 0.25))

		points(x=tros_rbrtabs[[i]]$nsam, y=tros_rbrtabs[[i]]$rich, pch=20, cex=0.1, col=cols[i])

		# points(x=interpolated$nsam, y=interpolated$rich, col=cols[i], lwd=3, type="l")
		points(x=extrapolated$nsam, y=extrapolated$rich, col=cols[i], lwd=3, type="l", lty=2)
		points(x=max(interpolated$nsam), y=max(interpolated$rich), col=cols[i], pch=pchs[i], cex=3)

		txt_x <- max(extrapolated$nsam) / 2
		ymin <- min(interpolated$rich)
		yrange <- max(extrapolated$rich) - min(interpolated$rich)
		txt_y1 <- ymin + yrange * (1/8)
		txt_y2 <- ymin + yrange * (2/8)
		txt_y3 <- ymin + yrange * (3/8)
		fcex <- 1.3
		text(txt_x, txt_y3, paste0("HCI: ", percent_hci, "%"), 
			adj=c(0, 0.5), col=cols[i], family="mono", cex=fcex)
		text(txt_x, txt_y2, paste0("EST: ", percent_est, "%"),
			adj=c(0, 0.5), col=cols[i], family="mono", cex=fcex, font=2)
		text(txt_x, txt_y1, paste0("LCI: ", percent_lci, "%"),
			adj=c(0, 0.5), col=cols[i], family="mono", cex=fcex)

		tablelist[[i]] <- data.frame(x=tros[i], total_observed_rich=max(tros_rbrtabs[[i]]$rich), mm_rich=dd, nsamp=max(tros_rbrtabs[[i]]$nsam))

	}
	dev.off()

	if(remake){
		save(list="tros_rbrtabs", file="habs_rbrtabs.rdata")
	}
	trostab <- do.call("rbind", tablelist)
	write.table(file="habs_table.txt", trostab, sep="\t", row.names=FALSE, quote=FALSE)








