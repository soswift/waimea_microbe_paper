#!/usr/bin/env Rscript

# goal is to compare paired beta-div between TERRESTRIAL and RIVERINE
	# samples at each distance from shore


# greet
	message("it's started!")

# load stuff
	library("data.table")
	library("vegan")
	library("parallel")
	library("Rcpp")
	library("fields")
	library("ggplot2")
	library("compositions")

	#source("../../rare_analysis_SAR/rare_analysis_sar_funs.r")

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


# (end standard-ish copypasted data organization)


# subset for only riverine and terrestrial samples
	message("subset for only riverine and terrestrial samples")
	otutable <- otutable[, metadata$habitat %in% c("Riverine", "Terrestrial")]
	metadata <- metadata[metadata$habitat %in% c("Riverine", "Terrestrial"), ]
	# double check that worked
	message("Next line better be TRUE (checking for metadata/otutable sort):")
	all(colnames(otutable) == metadata$sequencing_id)
	# get rid of empty otus
	message("Getting rid of empty OTUs: dim(otutable) before -> after")
	dim(otutable)
	otutable <- otutable[rowSums(otutable) > 0, ]
	dim(otutable)
	message("transform otutable to prop abund")
	otutable <- apply(X=otutable, MAR=2, FUN=function(x){x/sum(x)})

# functions for calculating bray
	message("calculating beta div matrix")
	cppFunction("
		float bray2samp(const NumericVector& i, const NumericVector& j){
			float out = 0;
			float c = 0;
			int n = i.size();
			float sumi = 0;
			float sumj = 0;
			// calculate c
			for(int k = 0; k<n; k++){
				sumi += i[k];
				sumj += j[k];
				if(i[k] < j[k]){
					c += i[k];
				}else{
					c += j[k];
				}
			}
			out = 1 - (2*c) / (sumi + sumj);
			return(out);
		}
	")

	parallel_bray_between <- function(rowsamps, colsamps, data, cols="species", ncores=10){
		require(parallel)
		require(vegan)
		if(cols=="species"){
			# orientation OK
			if(! all(rowsamps %in% rownames(data))){stop("Not all rowsamps in data")}
			if(! all(colsamps %in% rownames(data))){stop("Not all rowsamps in data")}
		}else if(cols=="samples"){
			if(! all(rowsamps %in% colnames(data))){stop("Not all rowsamps in data")}
			if(! all(colsamps %in% colnames(data))){stop("Not all rowsamps in data")}
			# transpose:
			data <- t(data)
		}else{
			stop("cols must be either 'species' or 'samples'")
		}
		# create list of sample combinations for which to compute beta div
		# this list is in COLUMN order, the same way R fills matrices by default
		# each item in list is a tuple of c(row_sample, col_sample).
		combo_list <- lapply(X=colsamps, FUN=function(cs){ lapply(X=rowsamps, FUN=function(rs){ c(rs, cs) }) })
		combo_list <- do.call("c", combo_list)
		# transform combo list into INDICES instead of names
		combo_list <- lapply(X=combo_list, FUN=function(samps){c(
			which(rownames(data)==samps[1]), which(rownames(data)==samps[2]))})

		if(ncores > 1){
			betas <- mclapply(
				X=combo_list, 
				FUN=function(samps){bray2samp(data[samps[1],], data[samps[2],])},
				mc.cores=ncores
			)
		}else{
			betas <- lapply(
				X=combo_list, 
				FUN=function(samps){bray2samp(data[samps[1],], data[samps[2],])}
			)
		}

		output <- matrix(unlist(betas), ncol=length(colsamps))
		rownames(output) <- rowsamps
		colnames(output) <- colsamps

		return(output)
	}

# get beta div mats by site, for each of the three trophic levels

	troph <- "Consumer"                          # edit this line
		betas <- list()
		sites <- unique(metadata$site_name)
		for(i in 1:length(sites)){
			site <- sites[i]
			rs <- metadata$sequencing_id[
					metadata$site_name==site & 
					metadata$habitat=="Riverine" &
					metadata$trophic==troph]
			cs <- metadata$sequencing_id[
					metadata$site_name==site & 
					metadata$habitat=="Terrestrial" &
					metadata$trophic==troph]
			if(length(rs) > 2 && length(cs) > 2){
				betas[[i]] <- parallel_bray_between(
					rowsamps = rs,
					colsamps = cs,
					data = otutable,
					cols = "samples",
					ncores = 10
				)
			}else{
				betas[[i]] <- NA
			}
			message(paste0(i, "/", length(sites)))
		}
		names(betas) <- sites
		consumer_betas <- betas                       # edit this line

	troph <- "Environmental"                          # edit this line
		betas <- list()
		sites <- unique(metadata$site_name)
		for(i in 1:length(sites)){
			site <- sites[i]
			rs <- metadata$sequencing_id[
					metadata$site_name==site & 
					metadata$habitat=="Riverine" &
					metadata$trophic==troph]
			cs <- metadata$sequencing_id[
					metadata$site_name==site & 
					metadata$habitat=="Terrestrial" &
					metadata$trophic==troph]
			if(length(rs) > 2 && length(cs) > 2){
				betas[[i]] <- parallel_bray_between(
					rowsamps = rs,
					colsamps = cs,
					data = otutable,
					cols = "samples",
					ncores = 10
				)
			}else{
				betas[[i]] <- NA
			}
			message(paste0(i, "/", length(sites)))
		}
		names(betas) <- sites
		environmental_betas <- betas                       # edit this line

	troph <- "PrimaryProducer"                         # edit this line
		betas <- list()
		sites <- unique(metadata$site_name)
		for(i in 1:length(sites)){
			site <- sites[i]
			rs <- metadata$sequencing_id[
					metadata$site_name==site & 
					metadata$habitat=="Riverine" &
					metadata$trophic==troph]
			cs <- metadata$sequencing_id[
					metadata$site_name==site & 
					metadata$habitat=="Terrestrial" &
					metadata$trophic==troph]
			if(length(rs) > 2 && length(cs) > 2){
				betas[[i]] <- parallel_bray_between(
					rowsamps = rs,
					colsamps = cs,
					data = otutable,
					cols = "samples",
					ncores = 10
				)
			}else{
				betas[[i]] <- NA
			}
			message(paste0(i, "/", length(sites)))
		}
		names(betas) <- sites
		pripro_betas <- betas                       # edit this line

	# save to save time later
	save(list=c("consumer_betas", "environmental_betas", "pripro_betas"), 
		file="betasbytroph.rdata")
	#load("betasbytroph.rdata")

# function to get distances between each site and the beach
	# d2b = dist 2 beach
	d2b <- function(sitename){
		sitename_lat <- metadata$latitude_deg[metadata$site_name == sitename][1]
		sitename_lng <- metadata$longitude_deg[metadata$site_name == sitename][1]
		beach_lat <- metadata$latitude_deg[metadata$site_name == "Beach"][1]
		beach_lng <- metadata$longitude_deg[metadata$site_name == "Beach"][1]
		# dist( matrix(c(sitename_lat, beach_lat, sitename_lng, beach_lng), ncol=2) )
		out <- fields::rdist.earth(
			x1 <- matrix(c(sitename_lng, sitename_lat), ncol=2),
			x1 <- matrix(c(beach_lng, beach_lat), ncol=2),
			miles=FALSE
		)
		return(as.numeric(out))
	}

# arrange all the data into one COLOSSAL data.frame
	betas2df <- function(betas, tr){
		lens <- sapply(X=betas, FUN=length)
		data.frame(
			beta=unlist(betas),
			site=unlist(mapply(FUN=rep, sites, lens, USE.NAMES=F)),
			trophic=tr,
			stringsAsFactors=FALSE
		)
	}
	plotdf <- rbind(
		betas2df(consumer_betas, "Consumer"),
		betas2df(environmental_betas, "Environmental"),
		betas2df(pripro_betas, "PrimaryProducer")
	)
	plotdf$dist2beach <- sapply(X=plotdf$site, FUN=d2b)
	plotdf$dist2beachlog10 <- log10( plotdf$dist2beach + 0.1 )

# aggregate means for simple plot
	plotdf_mean <- aggregate(plotdf$beta, by=list(plotdf$site, plotdf$trophic), FUN=mean)
	colnames(plotdf_mean) <- c("site", "trophic", "betamean")
	plotdf_mean$dist2beach <- sapply(X=plotdf_mean$site, FUN=d2b)
	plotdf_mean$dist2beachlog10 <- log10(plotdf_mean$dist2beach+0.1)
	lcifun <- function(x){if(length(x) > 5){ sort(x)[round(length(x) * 0.25)]}else{NA}}
	hcifun <- function(x){if(length(x) > 5){ sort(x)[round(length(x) * 0.75)]}else{NA}}
	plotdf_mean$lci <- aggregate(plotdf$beta, by=list(plotdf$site, plotdf$trophic), FUN=lcifun)$x
	plotdf_mean$hci <- aggregate(plotdf$beta, by=list(plotdf$site, plotdf$trophic), FUN=hcifun)$x


	# save so gamlss can be started from here if needed
	save(list=c("plotdf"), file="stuff_for_gamlss.rdata")


	# beta regression using gamlss
		# gamlss manual:
		# http://www.gamlss.com/wp-content/uploads/2013/01/gamlss-manual.pdf
		# one-inflated beta distribution reference:
		# http://www.gamlss.com/wp-content/uploads/2018/01/InflatedDistributioninR.pdf
	# control settings for gamlss
	library(gamlss)
	i_control_sett <- glim.control()
	i_control_sett$"c.cyc" <- 300
	i_control_sett$"bf.cyc" <- 300
	g_control_sett <- gamlss.control()
	g_control_sett$"n.cyc" <- 300

	if(FALSE){ # old parameterization
		mgam <- gamlss(
			formula=beta~dist2beachlog10*factor(trophic),
			sigma.formula=~dist2beachlog10*factor(trophic),
			nu.formula=~trophic,
			family="BEINF1", data=na.omit(plotdf),
			i.control=i_control_sett, control=g_control_sett
		)
	}

	plotdf$isConsumer <- plotdf$trophic == "Consumer"
	plotdf$isPrimaryProducer <- plotdf$trophic == "PrimaryProducer"

	mgam <- gamlss(
		formula=beta~dist2beachlog10*isConsumer+dist2beachlog10*isPrimaryProducer,
		sigma.formula=~dist2beachlog10*isConsumer+dist2beachlog10*isPrimaryProducer,
		nu.formula=~trophic,
		family="BEINF1", data=na.omit(plotdf),
		i.control=i_control_sett, control=g_control_sett
	)




	# write out gamlss results
	capture.output(summary(mgam), file="gamlss_results.txt")
	capture.output(Rsq(mgam), file="gamlss_results.txt", append=TRUE)

	

	# draw a simple plot 
	# ylims <- range(c(plotdf_mean$lci, plotdf_mean$hci), na.rm=TRUE)
	pdf("gam.pdf", useDingbats=FALSE)
	ylims <- range(c(plotdf_mean$lci, plotdf_mean$hci), na.rm=T)
	xlims <- range(plotdf_mean$dist2beachlog10, na.rm=TRUE)
	xfudge <- dist(xlims) * 0.01
	plot(0, type="n", xlim=xlims, ylim=ylims, xlab="Distance from beach (km)", 
		ylab="Marine vs terrestrial beta diversity", xaxt="n")
	axistix <- c(-1, -0.5, 0, 0.5, 1)
	axis(side=1, at=axistix, labels=round(10^axistix, 1))
	tros <- c("Environmental", "PrimaryProducer", "Consumer")
	trocols <- c("#fc8d62", "#66c2a5", "#756bb1")
	tropchs <- 15:17
	for(i in 1:3){
		xf <- (i - 1) * xfudge
		points(
			x=plotdf_mean$dist2beachlog10[plotdf_mean$trophic==tros[i]] + xf,
			y=plotdf_mean$betamean[plotdf_mean$trophic==tros[i]],
			col=trocols[i], pch=tropchs[i], cex=2
		)
		# segments(
		# 	x0=plotdf_mean$dist2beachlog10[plotdf_mean$trophic==tros[i]] + xf,
		# 	y0=plotdf_mean$lci[plotdf_mean$trophic==tros[i]],
		# 	y1=plotdf_mean$hci[plotdf_mean$trophic==tros[i]],
		# 	col=trocols[i], lwd=1
		# )


	}

	# make predictions for each trophic lvl and plot them
	predd2b <- seq(from=min(plotdf$dist2beachlog10, na.rm=T), to=max(plotdf$dist2beachlog10, na.rm=T), length.out=100)
	invlogit <- function(x){1/(1+exp(-1*x))}
	for(i in 1:3){
		# make pred df
		preddf <- data.frame(dist2beachlog10=predd2b, trophic=rep(tros[i], length(predd2b)))
		preddf$isConsumer <- ( preddf$trophic == "Consumer" )
		preddf$isPrimaryProducer <- ( preddf$trophic == "PrimaryProducer" )

		preds <- invlogit( predict(mgam, what="mu", newdata=preddf) )
		preds_sig <- invlogit( predict(mgam, what="sigma", newdata=preddf) )
		points(preds~predd2b, type="l", lwd=3, col=trocols[i])
	}
	dev.off()




	# redo the whole plot on logit scale
	# draw a simple plot 
	# ylims <- range(c(plotdf_mean$lci, plotdf_mean$hci), na.rm=TRUE)
	pdf("gam_logit.pdf", useDingbats=F)
	logit <- function(p){log(p/(1-p))}
	plotdf_mean$logitbeta <- logit(plotdf_mean$beta)
	ylims <- range(plotdf_mean$logitbeta, na.rm=T)
	xlims <- range(plotdf_mean$dist2beachlog10, na.rm=TRUE)
	xfudge <- dist(xlims) * 0.01
	plot(0, type="n", xlim=xlims, ylim=ylims, xlab="Distance from beach (km)", 
		ylab="logit(marine vs terrestrial beta diversity)", xaxt="n")
	axistix <- c(-1, -0.5, 0, 0.5, 1)
	axis(side=1, at=axistix, labels=round(10^axistix, 1))
	tros <- c("Environmental", "PrimaryProducer", "Consumer")
	trocols <- c("#fc8d62", "#66c2a5", "#756bb1")
	tropchs <- 15:17
	for(i in 1:3){
		xf <- (i - 1) * xfudge
		points(
			x=plotdf_mean$dist2beachlog10[plotdf_mean$trophic==tros[i]] + xf,
			y=plotdf_mean$logitbeta[plotdf_mean$trophic==tros[i]],
			col=trocols[i], pch=tropchs[i], cex=2
		)
		if(FALSE){
			segments(
				x0=plotdf_mean$dist2beachlog10[plotdf_mean$trophic==tros[i]] + xf,
				y0=logit(plotdf_mean$lci[plotdf_mean$trophic==tros[i]]),
				y1=logit(plotdf_mean$hci[plotdf_mean$trophic==tros[i]]),
				col=trocols[i], lwd=1
			)
		}


	}

	# plot logit models just for fun
	for(i in 1:3){ 
		preddf <- data.frame(dist2beachlog10=predd2b, trophic=rep(tros[i], length(predd2b)))
		preddf$isConsumer <- ( preddf$trophic == "Consumer" )
		preddf$isPrimaryProducer <- ( preddf$trophic == "PrimaryProducer" )

		points(
			x=predd2b,
			y=predict(mgam, what="mu", newdata=preddf),
			col=trocols[i], lwd=3, type="l")
	}
	dev.off()












# maybe useful for something else
	drawviolin <- function(data, xpos, scale=1, ...){
		if(length(data) > 5){
			dens <- density(data)
			# scale density to xwid/2
			dens$y <- dens$y / (2*sum(dens$y))
			dens$y <- dens$y * scale
			polyx <- c((xpos-dens$y), rev(xpos+dens$y))
			polyy <- c(dens$x, rev(dens$x))
			polygon(x=polyx, y=polyy, ...)
		}
	}
	drawviolins <- function(x, y, scale=1, ...){
		xs <- unique(x)
		ys <- lapply(xs, function(g){y[x==g]})
		trash <- mapply(drawviolin, ys, xs, scale, ...)
	}