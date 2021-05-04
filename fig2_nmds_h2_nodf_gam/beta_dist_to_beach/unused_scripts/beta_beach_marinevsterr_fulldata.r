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

	source("../../rare_analysis_SAR/rare_analysis_sar_funs.r")

# read in data
	otutable_fp <- "../rawdata/waimea_full/abundance_table_100.shared"
	metadata_fp <- "../rawdata/waimea_full/combined_prep_sample_metadata_16s.csv"
	taxonomy_fp <- "../rawdata/waimea_full/annotations_100.taxonomy"
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

# calculate beta div matrix

	# just for example

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
	# needs to be done item-wise, because otherwise you get "can't allocate vector" error.
	if(FALSE){
		site <- "Confluence"
		rowsamps <- metadata$sequencing_id[metadata$site_name==site & metadata$habitat=="Riverine"]
		colsamps <- metadata$sequencing_id[metadata$site_name==site & metadata$habitat=="Terrestrial"]
		data <- otutable
		cols <- "samples"
		ncores <- 10
	}


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

	if(FALSE){
		# test
		set.seed(12345)
		testmat <- otutable[,sample(1:ncol(otutable), 12)]

		vegdist_mat <- as.matrix(vegdist(t(testmat), method="bray"))
		jldcpp_mat <- parallel_bray_between(
			rowsamps=colnames(testmat),
			colsamps=colnames(testmat),
			data=t(testmat), ncores=1
		)
		plot(vegdist_mat ~ jldcpp_mat)
		abline(a=0, b=1, col="red")

	}

# get beta div mats by site
	betas <- list()
	sites <- unique(metadata$site_name)
	for(i in 1:length(sites)){
		site <- sites[i]
		betas[[i]] <- parallel_bray_between(
			rowsamps = metadata$sequencing_id[metadata$site_name==site & metadata$habitat=="Riverine"],
			colsamps = metadata$sequencing_id[metadata$site_name==site & metadata$habitat=="Terrestrial"],
			data = otutable,
			cols = "samples",
			ncores = 10
		)
		message(paste0(i, "/", length(sites)))
	}
	names(betas) <- sites
	save(list="betas", file="fulldata_betas.rdata")
	load("fulldata_betas.rdata")

# get distances between each site and the beach
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
	dists2beach <- sapply(X=sites, FUN=d2b)

# aggregate all the data into one COLOSSAL data.frame
	lens <- sapply(X=betas, FUN=length)
	plotdf <- data.frame(
		beta=unlist(betas),
		site=unlist(mapply(FUN=rep, sites, lens, USE.NAMES=F)),
		dist2beach=unlist(mapply(FUN=rep, dists2beach, lens, USE.NAMES=F)),
		stringsAsFactors=FALSE
	)






# plot mean and 95% CIs
	lcifun <- function(x){ sort(x)[round(length(x) * 0.025)] }
	hcifun <- function(x){ sort(x)[round(length(x) * 0.975)] }

	plotdf_agg <- data.frame(
		betamean = sapply(X=betas, FUN=mean),
		betalci = sapply(X=betas, FUN=lcifun),
		betahci = sapply(X=betas, FUN=hcifun),
		site = names(betas),
		dist2beach = sapply(X=names(betas), FUN=d2b)
	)

	m_all <- lm(beta~dist2beach, data=plotdf)
	m_all_summ <- summary(m_all)

	pdf("fulldata_beta_beach_allpairs_scatter.pdf")
	pcol <- "darkblue"
	plot(x=plotdf_agg$dist2beach, y=plotdf_agg$betamean, pch=20, 
		xlab="Distance to beach (km)", cex=3, col=pcol,
		ylab="Terrestrial vs. riverine Bray-Curtis dissimilarity",
		ylim=range(c(plotdf_agg$betahci, plotdf_agg$betalci)),
		main=paste0("Rsq=", round(m_all_summ$r.squared, 3)) )
	segments(x0=plotdf_agg$dist2beach, x1=plotdf_agg$dist2beach, 
		y0=plotdf_agg$betamean, y1=plotdf_agg$betahci, col=pcol)
	segments(x0=plotdf_agg$dist2beach, x1=plotdf_agg$dist2beach, 
		y0=plotdf_agg$betamean, y1=plotdf_agg$betalci, col=pcol)
	points(plotdf_agg$betahci ~ plotdf_agg$dist2beach, col=pcol, pch=20, cex=0.7)
	points(plotdf_agg$betalci ~ plotdf_agg$dist2beach, col=pcol, pch=20, cex=0.7)
	abline(m_all, lty=2)
	dev.off()


# ggplot style with violins

# make a plot
	p <- ggplot(plotdf, aes(x=dist2beach, y=beta, fill=factor(dist2beach))) +
		geom_violin() +
		scale_x_continuous(breaks=dists2beach, labels=round(dists2beach,1) )
	p <- p + scale_fill_manual(values=rep("gray", length(sites)))
	p <- p + ylim(range(c(plotdf_agg$betahci, plotdf_agg$betalci)))
	p <- p +  xlab("Distance to Beach (km)") + ylab("Terrestrial vs. riverine Bray-Curtis dissimilarity")
	p <- p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
	p <- p + theme(legend.position = "none")
	p <- p + geom_point(data = plotdf_agg, aes(x=dist2beach, y=betamean), size=3, color="darkblue")
	p <- p + geom_abline(intercept=coef(m_all)[1], slope=coef(m_all)[2], linetype="dashed")
	pdf("fulldata_beta_beach_allpairs_violin.pdf")
	p
	dev.off()


