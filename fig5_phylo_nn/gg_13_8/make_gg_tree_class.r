	library(ape)

# read stuff in
	ggtree <- ape::read.tree("97_otus.tree.gz")
	ggtax <- read.delim("97_otu_taxonomy.txt.gz", sep="\t", head=FALSE, stringsAsFactors=FALSE)
	colnames(ggtax) <- c("tip", "tax")

# check something
	all(ggtax$tip %in% ggtree$tip.label)

# split taxonomy
 	ggtax$tax <- gsub(pattern="\\[", replacement="", x=ggtax$tax)
 	ggtax$tax <- gsub(pattern="\\]", replacement="", x=ggtax$tax)

	tax_split <- strsplit(ggtax$tax, split="; ")
	# are all splits length 7?
	all(lapply(X=tax_split, FUN=length) == 7)
		# TRUE
	tax_split <- do.call("rbind", tax_split)
	colnames(tax_split) <- c("k","p","c","o","f","g","s")
	tax_split <- sub(pattern=".__", replacement="", x=tax_split)

	ggtax <- data.frame(ggtax, tax_split, stringsAsFactors=FALSE)
	rm(tax_split)

# drop tax with no class info
	ggtax <- ggtax[ggtax$c != "", ]

# function to get a random representative of a given class
	get_random_rep <- function(cls){ sample(as.character(ggtax$tip[ggtax$c == cls]), size=1) }

# function to make a tree with one random member from each class
	random_class_tree <- function(classes=unique(ggtax$c)){
		reps <- sapply(X=classes, FUN=get_random_rep)
		out <- ape::keep.tip(ggtree, tip=as.character(reps))
		rep2class <- function(rep){ classes[reps==rep] }
		out$tip.label <- sapply(X=out$tip.label, FUN=rep2class)
		return(out)
	}

# make a bunch of random class trees and consensus them together
	rclasstrees <- replicate(500, random_class_tree(), simplify=FALSE)
	cons <- ape::consensus(rclasstrees, p=0.5, check.labels=TRUE)

# write out
	write.tree(cons, file="gg_class_tree.nwk")