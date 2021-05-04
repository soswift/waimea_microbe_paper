#!/usr/bin/env Rscript

# nohup ./rich.r >status.txt &

# greet
	message("it's started!")

# load stuff
	library("data.table")
	library("parallel")
	library("ggplot2")
	library("gridExtra")
	library("car")

	source("../../rare_analysis_SAR/rare_analysis_sar_funs.r")

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

# function to get richness
	library(iNEXT)
	library(parallel)
	get_rich <- function(x){
		x <- x[x > 0]
		a <- ChaoRichness(x, datatype="abundance")
		return(c(observed=a$Observed, estimator=a$Estimator))
	}
	otutable_list <- as.data.frame((otutable))

	richs <- mclapply(X=otutable_list, FUN=get_rich, mc.cores=20)
	richs <- do.call("rbind", richs)

# function to get median and lower+upper quartiles for data
	medquart <- function(x){
		qs <- quantile(x)
		return(data.frame(
			y=median(x), 
			ymin=qs[2],
			ymax=qs[4]
		))
	}

# make plot - observed richness, aggregated by Habitat
	ggdata <- data.frame(Habitat=metadata$habitat, Trophic=metadata$trophic, richs, stringsAsFactors=F)
	plotslist <- list()
	habs <- c(unique(ggdata$Habitat), "All")
	for(i in 1:length(habs)){
		hab <- habs[i]
		if(hab == "All"){
			habdata <- ggdata
		}else{
			habdata <- ggdata[ggdata$Habitat == hab,]
		}

		p <- ggplot(habdata, aes(x=Trophic, y=observed, fill=Trophic)) + geom_violin()
		p <- p + geom_jitter(shape=16, position=position_jitter(0.2), size=0.30)
		p <- p + stat_summary(fun.data=medquart, geom="errorbar", width=0.15, color="white")
		p <- p + stat_summary(fun.data=medquart, geom="point", color="white")
		# p <- p + geom_dotplot(binaxis="y", stackdir="center", dotsize=0.30)
		# order:                            consumer    env        pp
		p <- p + scale_fill_manual(values=c("#756bb1", "#fc8d62", "#66c2a5"))
		p <- p + theme(legend.position="none", axis.title.x=element_blank(), axis.text=element_text(size=8))
		p <- p + labs(y="Observed Richness")
		p <- p + ggtitle(hab)
		plotslist[[i]] <- p
		names(plotslist)[i] <- hab
	}
	pdf("richness_by_habitat.pdf")
	grid.arrange(plotslist[[1]], plotslist[[2]], plotslist[[3]], plotslist[[4]], nrow=2)
	dev.off()


	# plot grouped
	ggdata2 <- ggdata
	ggdata2$Habitat <- "All"
	ggdata2 <- rbind(ggdata, ggdata2)
	ggdata2$Habitat <- factor(ggdata2$Habitat, levels=c("Marine", "Riverine", "Terrestrial", "All"))
	ddg <- 0.6
	p <- ggplot(ggdata2, aes(x=Habitat, y=observed, fill=Trophic))
	p <- p +  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
	p <- p + geom_violin(trim = TRUE, position = position_dodge(ddg))
	p <- p + geom_point(pch = ".", position = position_jitterdodge(dodge.width=ddg, jitter.width=0.20), size=0.7)
	p <- p + geom_boxplot(width = 0.15, position = position_dodge(ddg), outlier.shape=NA, coef = 0, alpha=0, color="white")
	p <- p + theme(legend.position = c(0.65, 0.8))
	p <- p + theme(legend.background=element_rect(fill=NA, color=NA))
	p <- p + theme(legend.title = element_blank())
	p <- p + theme(aspect.ratio=0.6)
	# order:                            consumer    env        pp
	p <- p + scale_fill_manual(values=c("#756bb1", "#fc8d62", "#66c2a5"))
	p <- p + theme(axis.title.x=element_blank(), axis.text=element_text(size=8))
	p <- p + labs(y="Observed Richness")
	pdf("richness_by_habitat_grouped.pdf")
	p
	dev.off()

	# do an anova
	get_eta2 <- function(res){
		res <- as.data.frame(res)
		res$eta2 <- res$"Sum Sq" / sum(res$"Sum Sq")
		return(res)
	}
	m1 <- aov(log10(ggdata$observed) ~ ggdata$Habitat + ggdata$Trophic)
	res <- car::Anova(m1)
	res <- get_eta2(res)
	capture.output(res, file="anova_results.txt")
	pdf("anova_varexp_piechart.pdf")
	pie(x=res$eta2, labels=c("Habitat***", "Trophic***", "Residuals"))
	dev.off()


# are observed and chao1 just the same thing because sample depth is so good?
	pdf("richness_deficit.pdf")
	hist(ggdata$estimator - ggdata$observed, col="lightblue", xlab="Richness deficit (Chao1)", main="")
	dev.off()


# make plot - chao1 richness, aggregated by Trophic Level
	ggdata <- data.frame(Habitat=metadata$habitat, Trophic=metadata$trophic, richs, stringsAsFactors=F)
	plotslist <- list()
	tros <- c(unique(ggdata$Trophic), "All")
	for(i in 1:length(tros)){
		tro <- tros[i]
		if(tro == "All"){
			trodata <- ggdata
		}else{
			trodata <- ggdata[ggdata$Trophic == tro,]
		}

		p <- ggplot(trodata, aes(x=Habitat, y=observed, fill=Habitat)) + geom_violin()
		p <- p + geom_jitter(shape=16, position=position_jitter(0.2), size=0.30)
		p <- p + stat_summary(fun.data=medquart, geom="errorbar", width=0.15, color="white")
		p <- p + stat_summary(fun.data=medquart, geom="point", color="white")
		# p <- p + geom_dotplot(binaxis="y", stackdir="center", dotsize=0.30)
		# order:                             marine     river      terr
		p <- p + scale_fill_manual(values=c("#1f78b4", "#a6cee3", "#fed9a6"))
		p <- p + theme(legend.position="none", axis.title.x=element_blank(), axis.text=element_text(size=8))
		p <- p + labs(y="Observed Richness")
		p <- p + ggtitle(tro)
		plotslist[[i]] <- p
		names(plotslist)[i] <- tro
	}
	pdf("richness_by_trophic.pdf")
	grid.arrange(plotslist[[1]], plotslist[[2]], plotslist[[3]], plotslist[[4]], nrow=2)
	dev.off()

	# plot grouped
	ggdata2 <- ggdata
	ggdata2$Trophic <- "All"
	ggdata2 <- rbind(ggdata, ggdata2)
	ggdata2$Trophic <- factor(ggdata2$Trophic, levels=c("Environmental", "PrimaryProducer", "Consumer", "All"))
	ddg <- 0.6
	p <- ggplot(ggdata2, aes(x=Trophic, y=observed, fill=Habitat))
	p <- p +  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
	p <- p + geom_violin(trim = TRUE, position = position_dodge(ddg))
	p <- p + geom_point(pch = ".", position = position_jitterdodge(dodge.width=ddg, jitter.width=0.20), size=0.7)
	p <- p + geom_boxplot(width = 0.15, position = position_dodge(ddg), outlier.shape=NA, coef = 0, alpha=0, color="white")
	p <- p + theme(legend.position = c(0.65, 0.8))
	p <- p + theme(legend.background=element_rect(fill=NA, color=NA))
	p <- p + theme(legend.title = element_blank())
	p <- p + theme(aspect.ratio=0.6)
	# order:                             marine     river      terr
	p <- p + scale_fill_manual(values=c("#1f78b4", "#a6cee3", "#fed9a6"))
	p <- p + theme(axis.title.x=element_blank(), axis.text=element_text(size=8))
	p <- p + labs(y="Observed Richness")
	pdf("richness_by_trophic_grouped.pdf")
	p
	dev.off()

	#save so plot can be done remotely
	save(list=c("ggdata", "ggdata2"), file="ggdata_trophic_grouped.rdata")

	# plot grouped, log10
	ggdata2 <- ggdata
	ggdata2$Trophic <- "All"
	ggdata2 <- rbind(ggdata, ggdata2)
	ggdata2$Trophic <- factor(ggdata2$Trophic, levels=c("Environmental", "PrimaryProducer", "Consumer", "All"))
	ggdata2$observed <- log10(ggdata2$observed)
	ddg <- 0.6
	p <- ggplot(ggdata2, aes(x=Trophic, y=observed, fill=Habitat))
	p <- p +  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
	p <- p + geom_violin(trim = TRUE, position = position_dodge(ddg))
	p <- p + geom_point(pch = ".", position = position_jitterdodge(dodge.width=ddg, jitter.width=0.20), size=0.7)
	p <- p + geom_boxplot(width = 0.15, position = position_dodge(ddg), outlier.shape=NA, coef = 0, alpha=0, color="white")
	p <- p + theme(legend.position = c(0.2, 0.15))
	p <- p + theme(legend.background=element_rect(fill=NA, color=NA))
	p <- p + theme(legend.title = element_blank())
	p <- p + theme(aspect.ratio=0.6)
	# order:                             marine     river      terr
	p <- p + scale_fill_manual(values=c("#1f78b4", "#a6cee3", "#fed9a6"))
	p <- p + theme(axis.title.x=element_blank(), axis.text=element_text(size=8))
	p <- p + labs(y="Log10 Observed Richness")
	pdf("log10richness_by_trophic_grouped.pdf")
	p
	dev.off()

