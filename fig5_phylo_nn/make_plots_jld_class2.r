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
	ggtree_class <- read.tree("gg_13_8/gg_class_tree.nwk")
	# check all classes are in tree
	all(unique(plotdata$Class) %in% ggtree_class$tip.label)
	# they aren't ... just drop extras
	plotdata <- plotdata[plotdata$Class %in% ggtree_class$tip.label, ]

# number of genera per class
	gen_count_df <- data.frame(Class=unique(plotdata$Class), stringsAsFactors=FALSE)
	gen_count_df$ngenera <- sapply(X=gen_count_df$Class, FUN=function(cls){
		length(unique(plotdata$Genus[plotdata$Class==cls]))
	})




	# also get rid of classes with less than 300 EMP observations and 300 Waimea observations
	classcountsdf <- data.frame(
		Class=as.vector(unique(plotdata$Class)), 
		nEMP=0, nWAI=0, stringsAsFactors=FALSE
	)
	for(i in 1:nrow(classcountsdf)){
		cl <- classcountsdf$Class[i]
		classcountsdf$nEMP[i] <- sum(plotdata$Class == cl & plotdata$type == "EMP")
		classcountsdf$nWAI[i] <- sum(plotdata$Class == cl & plotdata$type == "Waimea")
	}
	classes2keep <- classcountsdf$Class[classcountsdf$nEMP>300 & classcountsdf$nWAI>300]
	
	plotdata <- plotdata[plotdata$Class %in% classes2keep, ]


	# also prune extra tips
	ggtree_class <- keep.tip(ggtree_class, tip=unique(plotdata$Class))

	# rename data for plot function
	plotdata_less <- data.frame(
		label=plotdata$Class,
		x=plotdata$dist,
		g=plotdata$type
	)

	# read in P-values and format
	pvals_df <- read.delim("mahdi_p_values_fdr/class_fdr_corrected.csv", sep=",", header=T)
	pvals_df$class <- gsub(".__", "", pvals_df$class)
	pvals_df$class <- gsub("\\]", "", pvals_df$class)
	pvals_df$class <- gsub("\\[", "", pvals_df$class)

	# subset for significant classes
	plotdata_less <- plotdata_less[plotdata_less$label %in% pvals_df$class, ]
	plotdata_less$label <- as.character(plotdata_less$label)
	ggtree_class <- keep.tip(ggtree_class, unique(plotdata_less$label))

	# make plots
	source("plottreeab.r")


	pdf("treeplots.pdf", useDingbats = FALSE)

	plottreeab(ggtree_class, plotdata_less, 
		xlabel="Phylogenetic distance to closest EMP ASV", 
		type="quartiles", colbg="gray")

	plottreeab(ggtree_class, plotdata_less, 
		xlabel="Phylogenetic distance to closest EMP ASV", 
		type="quartiles2", colbg="gray")

	plottreeab(ggtree_class, plotdata_less, 
		xlabel="Phylogenetic distance to closest EMP ASV", 
		type="ovals")


	plottreeab(ggtree_class, plotdata_less, 
		xlabel="Phylogenetic distance to closest EMP ASV", 
		type="violins")


	plottreeab(ggtree_class, plotdata_less, 
		xlabel="Phylogenetic distance to closest EMP ASV", 
		type="dots", pointcex=0.7)

	plottreeab(ggtree_class, plotdata_less, 
		xlabel="Phylogenetic distance to closest EMP ASV", 
		type="overlap", pointcex=0.7, trans=0.10)

	plottreeab(ggtree_class, plotdata_less, 
		xlabel="Phylogenetic distance to closest EMP ASV", 
		type="joy", pointcex=0.7, trans=0.10, cola="red")


	dev.off()

