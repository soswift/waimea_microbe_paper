#!/usr/bin/env Rscript

# goal is to compare paired beta-div between TERRESTRIAL and RIVERINE
	# samples at each distance from shore


# read in data
	metadata_fp <- "../rawdata/waimea_full/combined_prep_sample_metadata_16s.csv"
	metadata <- read.delim(metadata_fp, sep=",", header=T, stringsAsFactors=FALSE)
	# Sean's NMDS
	nmds <- readRDS("all_sample_ordinations.rds")
	nmdspoints <- nmds$points

# organize plotty data
	# check compatibility
	all(rownames(nmdspoints) %in% metadata$sequencing_id)
	# OK
	# subset metadata then sort 2 match
	metadata <- metadata[metadata$sequencing_id %in% rownames(nmdspoints), ]
	data <- metadata[order(metadata$sequencing_id), ]
	nmdspoints <- nmdspoints[order(rownames(nmdspoints)), ]
	all(data$sequencing_id == rownames(nmdspoints))
	data$nmds1 <- nmdspoints[,1]
	data$nmds2 <- nmdspoints[,2]

# clean up habitat and trophic variables
	data$trophic[data$trophic %in% c("Omnivore" ,"Herbivore", "Detritivore", "Carnivore")] <- "Consumer"
	data$trophic[data$trophic == "PrimaryProducer"] <- "Primary producer"

# add pch column to data
	data$pch <- ""
	data$pch[data$empo_2 == "Animal"] <- "♥"
	data$pch[data$empo_2 == "Plant"] <- "☘"
	data$pch[data$empo_2 == "Fungus"] <- "☂"
	data$pch[data$empo_2 == "Saline"] <- "▲"
	data$pch[data$empo_2 == "Non-saline"] <- "●"

# add habitat color to data
	data$habcol <- ""
	data$habcol[data$habitat == "Marine"] <- "#1f78b4"
	data$habcol[data$habitat == "Riverine"] <- "#a6cee3"
	data$habcol[data$habitat == "Terrestrial"] <- "#fed9a6"

# add trophic color to data
	data$trocol <- ""
	data$trocol[data$trophic == "Primary producer"] <- "#66c2a5"
	data$trocol[data$trophic == "Consumer"] <- "#756bb1"
	data$trocol[data$trophic == "Environmental"] <- "#fc8d62"

# subsets
	data_mar <- data[data$habitat == "Marine", ]
	data_riv <- data[data$habitat == "Riverine", ]
	data_ter <- data[data$habitat == "Terrestrial", ]


# plot everything, base R
	cairo_pdf("plots_baseR_jld.pdf", onefile=TRUE)
	plot(data$nmds2, data$nmds1, type="n", main="All")
	text(x=data$nmds2, y=data$nmds1, labels=data$pch, col=data$habcol)

	# plot marine only
	plot(data_mar$nmds2, data_mar$nmds1, type="n")
	text(x=data_mar$nmds2, y=data_mar$nmds1, labels=data_mar$pch, col=data_mar$trocol)

	# plot riverine only
	plot(data_riv$nmds2, data_riv$nmds1, type="n")
	text(x=data_riv$nmds2, y=data_riv$nmds1, labels=data_riv$pch, col=data_riv$trocol)

	# plot terrestrial only
	plot(data_ter$nmds2, data_ter$nmds1, type="n")
	text(x=data_ter$nmds2, y=data_ter$nmds1, labels=data_ter$pch, col=data_ter$trocol)
	dev.off()

# plot everything, ggplot
	library(ggplot2)

	noaxis <- theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
		axis.ticks.x=element_blank(), axis.title.y=element_blank(), 
		axis.text.y=element_blank(), axis.ticks.y=element_blank())
	zerolines <- geom_hline(0) + geom_vline(0)

	p_all <- ggplot(data=data, aes(x=nmds2, y=nmds1, color=habitat, shape=empo_2)) + 
		#geom_hline(color="lightgray", yintercept=0) + 
		#geom_vline(color="lightgray", xintercept=0) + 
		geom_point(size=3) +
		scale_shape_manual(values=c("♥", "☂", "●", "☘", "▲" )) + 
		scale_color_manual(values=c("#1f78b4", "#a6cee3", "#fed9a6")) +
		theme_classic() + theme(legend.position = "none") + noaxis

	p_mar <- ggplot(data=data_mar, aes(x=nmds2, y=nmds1, color=trophic, shape=empo_2)) + 
		geom_point(size=3) +
		scale_shape_manual(values=c("♥", "☘", "▲" )) + 
		scale_color_manual(values=c("#756bb1", "#fc8d62", "#66c2a5")) +
		theme_classic() + theme(legend.position = "none") + noaxis

	p_riv <- ggplot(data=data_riv, aes(x=nmds2, y=nmds1, color=trophic, shape=empo_2)) + 
		geom_point(size=3) +
		scale_shape_manual(values=c("♥", "●", "☘")) + 
		scale_color_manual(values=c("#756bb1", "#fc8d62", "#66c2a5")) +
		theme_classic() + theme(legend.position = "none") + noaxis

	p_ter <- ggplot(data=data_ter, aes(x=nmds2, y=nmds1, color=trophic, shape=empo_2)) + 
		geom_point(size=3) +
		scale_shape_manual(values=c("♥", "☂", "●", "☘")) + 
		scale_color_manual(values=c("#756bb1", "#fc8d62", "#66c2a5")) +
		theme_classic() + theme(legend.position = "none") + noaxis

	p_blank <- ggplot() + theme_classic()



	library(gridExtra)



	grid.arrange(p_all, p_mar,p_riv,p_ter, ncol=2, nrow=2)


	cairo_pdf("nmds_plots_combined_gg_jld.pdf")
	lm <- matrix(c(1,1,1,1,2,2,2,3,2,2,2,4,2,2,2,5), nrow=4, byrow=T)
	p_final <- grid.arrange(p_blank, p_all, p_mar,p_riv,p_ter, 
		layout_matrix=lm)
	dev.off()
