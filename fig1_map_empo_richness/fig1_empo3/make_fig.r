# packages
	library(specificity)
	library(ggplot2)
	library(ggtree)
	library(ggstance)


# read stuff in
	metadata_fp <- "../rawdata/waimea_full/combined_prep_sample_metadata_16s.csv"
	metadata <- read.delim(metadata_fp, sep=",", header=T, stringsAsFactors=FALSE)

# drop stuff from metadata where empo or habitat is missing
	metadata <- metadata[metadata$empo_3 != "", ]
	metadata <- metadata[metadata$habitat != "", ]

# make empo3 tree
	empo3data <- data.frame(empo1=metadata$empo_1, empo2=metadata$empo_2, empo3=metadata$empo_3)
	empo3data <- apply(X=empo3data, MARGIN=2, FUN=function(x){gsub(pattern="\\(saline\\)", replacement="sal", x=x)})
	empo3data <- apply(X=empo3data, MARGIN=2, FUN=function(x){gsub(pattern="\\(non-saline\\)", replacement="nsal", x=x)})
	empo3data <- apply(X=empo3data, MARGIN=2, FUN=function(x){gsub(pattern=" ", replacement="_", x=x)})
	tr <- read.tree(text=onto2nwk(empo3data))
	# fix tree tips
	tr$tip.label <- sub(pattern="nsal", replacement="(ns)", x=tr$tip.label)
	tr$tip.label <- sub(pattern="sal", replacement="(s)", x=tr$tip.label)
	tr$tip.label <- gsub(pattern="_", replacement=" ", x=tr$tip.label)

# plain empo3 tree for node label reference later
	pdf("plain_empo3_tree.pdf")
	plot(tr, show.node.label=T)
	dev.off()

# change metadata empo3 to match tr
	metadata$empo_3 <- sub(pattern="\\(saline\\)", replacement="(s)", x=metadata$empo_3)
	metadata$empo_3 <- sub(pattern="\\(non-saline\\)", replacement="(ns)", x=metadata$empo_3)

# aggregate metadata into habitat counts by empo3
	metadata$one <- 1
	md_agg <- aggregate(metadata$one, by=list(metadata$habitat, metadata$empo_3), FUN=sum)
	colnames(md_agg) <- c("habitat", "id", "count")

# use ggtree to add stacked bars
	# DOESN'T WORK, DON'T HAVE TIME FOR THIS SHIT
	# p1 <- ggtree(tr) + geom_tiplab()
	# p2 <- facet_plot(p1, panel="stacked", data=metadata, 
	# 	geom=md_agg, aes(x=count, fill=as.factor(habitat)))

# aggregate to proportions
	md_totals <- aggregate(metadata$one, by=list(metadata$empo_3), FUN=sum)
	colnames(md_totals) <- c("empo_3", "total")
	md_totals$Terrestrial <- aggregate(metadata$habitat=="Terrestrial", by=list(metadata$empo_3), FUN=sum)$x
	md_totals$Riverine    <- aggregate(metadata$habitat=="Riverine",    by=list(metadata$empo_3), FUN=sum)$x
	md_totals$Marine      <- aggregate(metadata$habitat=="Marine",      by=list(metadata$empo_3), FUN=sum)$x
	md_totals$TerProp <- md_totals$Terrestrial / md_totals$total
	md_totals$RivProp <- md_totals$Riverine / md_totals$total
	md_totals$MarProp <- md_totals$Marine / md_totals$total

	# md_totals$empo_3 <- gsub(pattern="_", replacement=" ", x=md_totals$empo_3)

# use regular R to do it
	pdf("empo_sample_visualization.pdf")
	par(mfrow=c(1,3))
	# tree
		par(mar=c(4, 0, 0, 0))
		plot(tr, show.node.label=F, cex=1.5, edge.width=2 )
		mtext("EMP ontology", line=2, side=1)
	# proportions by habitat
		par(mar=c(4, 0, 0, 0))
		plot(x=0, type="n", ylim=c(1, length(tr$tip.label)), 
			xlim=c(0,1), xaxt="n", yaxt="n", bty="n", xlab="", ylab="")
		for(i in 1:nrow(md_totals)){
			empo3_i <- tr$tip.label[i]
			row_i <- which(md_totals$empo_3 == empo3_i)
			bp2 <- md_totals$TerProp[row_i]
			bp3 <- bp2 + md_totals$RivProp[row_i]
			# terrestrial
			rect(
				xleft=0, xright=bp2, ybottom=i-0.5,	ytop=i+0.5,
				col="#a6611a", lwd=2
			)
			# riverine
			rect(
				xleft=bp2, xright=bp3, ybottom=i-0.5, ytop=i+0.5,
				col="#a6cee3", lwd=2
			)
			# marine
			rect(
				xleft=bp3, xright=1, ybottom=i-0.5, ytop=i+0.5,
				col="#1f78b4", lwd=2
			)
		}
		axis(side=1, at=c(0, 0.33, 0.66, 1), labels=c("0", "1/3", "2/3", "1"),
			tick=FALSE, pos=0.75)
		mtext("Habitat proportion", line=2, side=1)
	# total number of samples
	log <- FALSE
	if(!log){
		par(mar=c(4, 0, 0, 0))
		plot(x=0, type="n", ylim=c(1, length(tr$tip.label)), 
			xlim=c(0,max(md_totals$total)), xaxt="n", yaxt="n", bty="n", 
			xlab="", ylab="")
		for(i in 1:nrow(md_totals)){
			empo3_i <- tr$tip.label[i]
			row_i <- which(md_totals$empo_3 == empo3_i)
			rect(xleft=0, xright=md_totals$total[row_i],
				ybottom=i-0.5, ytop=i+0.5, col="gray", lwd=2)
		}
		axis(side=1, at=c(0, 100, 200, 300, 400), tick=FALSE, pos=0.75)
		mtext("Number of samples", line=2, side=1)
	}else{
		par(mar=c(4, 0, 0, 0))
		md_totals$log10total <- log10(md_totals$total)
		plot(x=0, type="n", ylim=c(1, length(tr$tip.label)), 
			xlim=c(0,max(md_totals$log10total)), xaxt="n", yaxt="n", bty="n", 
			xlab="", ylab="")
		for(i in 1:nrow(md_totals)){
			empo3_i <- tr$tip.label[i]
			row_i <- which(md_totals$empo_3 == empo3_i)
			rect(xleft=0, xright=md_totals$log10total[row_i],
				ybottom=i-0.5, ytop=i+0.5, col="gray", lwd=2)
		}
		axis(side=1, at=c(0, 1, 2), labels=c(1, 10, 100), tick=FALSE, pos=0.75)
		mtext("Number of samples", line=2, side=1)

	}
	dev.off()