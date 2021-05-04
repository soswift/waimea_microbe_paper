# set up env
	library(ggplot2)
	source("plottreeab.r")
	library("ape")


	plotdata <- read.delim("data_jld_organized.txt.gz", sep="\t", stringsAsFactors=FALSE)

# add in type
	plotdata$type <- ""
	plotdata$type[grepl("waimea", plotdata$Seq)] <- "Waimea"
	plotdata$type[grepl("emp", plotdata$Seq)] <- "EMP"


# re-create mahdi figure
	col2trans <- function(col, alpha){
		rgb(
			red=col2rgb(col)[1]/255,
			green=col2rgb(col)[2]/255,
			blue=col2rgb(col)[3]/255,
			alpha
		)
	}
	pdf("dist_density_empvswai.pdf", useDingbats=FALSE)
	col_emp <- "orange"
	col_wai <- "blue"
	opacity <- 0.4
	trans_emp <- col2trans(col_emp, opacity)
	trans_wai <- col2trans(col_wai, opacity)

	emp_density <- density(plotdata$dist[plotdata$type == "EMP"])
	wai_density <- density(plotdata$dist[plotdata$type == "Waimea"])

	transblue <- rgb(0,0,1,alpha=opacity,maxColorValue=1)
	transred  <- rgb(1,0,0,alpha=opacity,maxColorValue=1)
	ymin <- 0
	ymax <- max(emp_density$y) + 0.02 * max(emp_density$y)
	ymid <- mean(c(ymin,ymax))
	plot(0, type="n", xlab="Phylogenetic distance to nearest EMP neighbor", ylab="",
		yaxt="n", xlim=c(0, 1), main="", ylim=c(ymin,ymax), xaxs="i", yaxs="i")
	polygon(x=emp_density$x, y=emp_density$y, col=trans_emp, lwd=3, border=NA)
	polygon(x=wai_density$x, y=wai_density$y, col=trans_wai,  lwd=3, border=NA)
	points(emp_density, col=col_emp, lwd=3, type="l")
	points(wai_density, col=col_wai, lwd=3, type="l")
	legendys <- c(ymid-0.05*(ymax-ymin), ymid+0.15*(ymax-ymin) )
	points(x=c(0.5, 0.5), y=legendys, pch=21, cex=6, lwd=3,
		bg=c(trans_wai, trans_emp), col=c(col_wai, col_emp))
	text(x=c(0.57, 0.57), y=legendys, labels=c("Waimea ASVs", "EMP ASVs"), 
		adj=c(0,0.5), cex=1)
	dev.off()



# make cool figure

	#  example data
	# tr <- rtree(30)
	# d <- data.frame(label=rep(tr$tip.label, each=100))
	# d$x <- rnorm(nrow(d))
	# d$g <- rep(c("a", "b"), length.out=nrow(d))
	# plottreeab(tr, d)

# load greengenes tree and format appropriately
	ggtree_fam <- read.tree("gg_13_8/gg_fam_tree.nwk")
	# check all fams are in tree
	all(unique(plotdata$Family) %in% ggtree_fam$tip.label)
	# they aren't ... just drop extras
	plotdata <- plotdata[plotdata$Family %in% ggtree_fam$tip.label, ]
	# remove blank
	plotdata <- plotdata[plotdata$Family != "", ]
	# also get rid of fams with less than 300 EMP observations and 300 Waimea observations
	familycountsdf <- data.frame(
		Family=as.vector(unique(plotdata$Family)), 
		nEMP=0, nWAI=0, stringsAsFactors=FALSE
	)
	for(i in 1:nrow(familycountsdf)){
		cl <- familycountsdf$Family[i]
		familycountsdf$nEMP[i] <- sum(plotdata$Family == cl & plotdata$type == "EMP")
		familycountsdf$nWAI[i] <- sum(plotdata$Family == cl & plotdata$type == "Waimea")
	}

	fams2keep <- familycountsdf$Family[familycountsdf$nEMP>300 & familycountsdf$nWAI>300]
	
	plotdata <- plotdata[plotdata$Family %in% fams2keep, ]


	# also prune extra tips
	ggtree_fam <- keep.tip(ggtree_fam, tip=unique(plotdata$Family))

	# rename data for plot function
	plotdata_less <- data.frame(
		label=plotdata$Family,
		x=plotdata$dist,
		g=plotdata$type
	)

	# find top effect sizes
	effect_size_df <- data.frame(label=unique(plotdata_less$label))
	effect_size_df$rawdiff <- sapply(X=effect_size_df$label, FUN=function(lab){
		wai_mean <- mean(plotdata_less$x[plotdata_less$label==lab & plotdata_less$g=="Waimea"])
		emp_mean <- mean(plotdata_less$x[plotdata_less$label==lab & plotdata_less$g=="EMP"])
		return(abs(wai_mean - emp_mean))
	})
	effect_size_df <- effect_size_df[order(effect_size_df$rawdiff, decreasing=T), ]

	# subset to top 50
	plotdata_less <- plotdata_less[plotdata_less$label %in% effect_size_df$label[1:50], ]
	plotdata_less$label <- as.character(plotdata_less$label)
	plotdata_less$g <- as.character(plotdata_less$g)
	ggtree_fam <- keep.tip(ggtree_fam, unique(plotdata_less$label))


	# make plots
	pdf("plot_fam.pdf", useDingbats = FALSE)

	plottreeab(ggtree_fam, plotdata_less, xlab="Phylogenetic distance to closest EMP ASV", 
		type="quartiles", colbg="gray")


	dev.off()

	
	# make plots
	pdf("plot_fam_bars.pdf", useDingbats = FALSE)

	plottreeab(ggtree_fam, plotdata_less, xlab="Phylogenetic distance to closest EMP ASV", 
		type="quartiles", colbg="gray")


	dev.off()


	# make plots
	pdf("plot_fam_dots.pdf", useDingbats = FALSE)

	plottreeab(ggtree_fam, plotdata_less, xlab="Phylogenetic distance to closest EMP ASV", 
		type="dots")


	dev.off()

