#!/usr/bin/env Rscript


# greet
  message("it's started!")

# read in data
  otutable_fp <- "../../data/topicdata/100topics_sample_topics.csv"
  metadata_fp <- "../../data/raw/sequence/waimea_16S_v1/combined_prep_sample_metadata_16s.csv"
  taxonomy_fp <- "../../data/raw/sequence/waimea_16S_v1/annotations_100.taxonomy"
  otutable <- read.delim(otutable_fp, sep=",", header=F, stringsAsFactors=FALSE)
  metadata <- read.delim(metadata_fp, sep=",", header=T, stringsAsFactors=FALSE)
  taxonomy <- read.delim(taxonomy_fp, sep="\t", header=T, stringsAsFactors=FALSE)

# format otutable
  rownames(otutable) <- otutable[,1]
  otutable <- otutable[,2:ncol(otutable)]


# get rid of stuff
  metadata <- metadata[metadata$habitat %in% c("Marine", "Riverine", "Terrestrial"), ]
  # remove samples from otutable not in metadata
  otutable <- otutable[rownames(otutable) %in% metadata$sequencing_id,]
  # remove otus that are not present in any samples
  otutable <- otutable[,colSums(otutable) > 0]


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

# make column of combined habitat + empo_3 categories
  metadata$habempo3 <- paste(metadata$habitat, metadata$empo_3)

# collapse otutable by habempo3 categories
  agg_by_habempo3 <- function(x){
    ag <- aggregate(x, by=list(metadata$habempo3), FUN=sum)
    out <- ag$x
    names(out) <- ag[,1]
    return(out)
  }
  # it's in transpose orientation... good! 
  habempo3table <- apply(X=otutable, MAR=1, FUN=agg_by_habempo3)
  message("otutable aggregation complete")

# format habempo3table to make binary
  habempo3table <- habempo3table > 0

# load in gen alg function
  source("draft.r")

# run algo for each n empo3s
  message("running gen alg")
  results_list <- list()
  for(i in 1:nrow(habempo3table)){
    results_list[[i]] <- genOptRich(
      otutable = habempo3table,
      popsize = 300,
      best2keep = 20,
      regs2keep = 20,
      gens2stop = 40,
      maxgens = 1000,
      mutrate = 3,
      noisy=TRUE,
      n=i
    )
    message(paste0(i, "/", nrow(habempo3table)))
  }

  save(list="results_list", file="results_list.rdata")


# create accumulation curve
ns <- sapply(X=results_list, FUN=function(x){sum(x$bestState)})
rs <- sapply(X=results_list, FUN=function(x){max(x$richnessByGen)})

# make plot
pdf("plot1_accumulation.pdf")
plot(rs ~ ns, type="o", pch=20, cex=2, xlab="Number of categories", ylab="Maximum richness")
dev.off()

# get data for plot 2

plotmat <- sapply(X=results_list, FUN=function(x){x$bestState})
rownames(plotmat) <- rownames(habempo3table)

# order plotmat nicely
plotmat <- plotmat[order(rowSums(plotmat)), ]

pdf("plot2_mat.pdf")
par(mar=c(4,9,1,1))
plot(0, type="n", 
  xlim=c(0, nrow(plotmat)),
  ylim=c(0, nrow(plotmat)),
  xlab="Number of categories",
  ylab="", yaxt="n", xaxt="n",
  xaxs="i", yaxs="i"
)
plotcol <- function(st){
  n <- sum(st)
  for(i in 1:length(st)){
    if(st[i] == TRUE){
      rect(
        xleft=n-1,
        xright=n,
        ybottom=i-1,
        ytop=i,
        col="black"
      )
    }
  }
}
shit <- apply(X=plotmat, MAR=2, FUN=plotcol)
n <- ncol(plotmat)
axis(side=1, at = (0:(n-1))+0.5, labels=1:n, cex.axis=0.7)
axis(side=2, at=(0:(n-1))+0.5, 
  labels=rownames(plotmat), 
  las=2, cex.axis=0.7)
abline(v=0:n, col="gray")
abline(h=0:n, col="gray")

dev.off()







pdf("combined_plot.pdf")

par(mfrow=c(2,1))

par(mar=c(0.5,9,3,1))

# make plot
nsp <- ns + 0.5
plot(rs ~ nsp, type="o", pch=20, cex=4, xlab="", 
  ylab="", xaxt="n", yaxt="n",  xaxs="i", xlim=c(1, nrow(plotmat)+1))
abline(v=2:nrow(plotmat), col="gray")
points(rs ~ nsp, type="l")
mtext(side=2, text="Maximum topic richness", line=5)
yaxvals <- (2:8)*100
yaxlabs <- as.character(yaxvals)

axis(side=2, at=yaxvals, labels=parse(text=yaxlabs), cex.axis=0.7, las=2)
# axis(side=3, at=nsp, labels=rs, cex.axis=0.7, las=2)
# fudge <- dist(range(rs))[1] * 0.025
# text(x=nsp, y=rs+fudge, labels=rs, adj=c(0.5, 0.5), cex=0.7)
text(x=nsp, y=rs, labels=rs, adj=c(0.5, 0.5), cex=0.7, col="white")

par(mar=c(4,9,0,1))


plot(0, type="n", 
  xlim=c(0, nrow(plotmat)),
  ylim=c(0, nrow(plotmat)),
  xlab="Number of categories",
  ylab="", yaxt="n", xaxt="n",
  xaxs="i", yaxs="i"
)
plotcol <- function(st){
  n <- sum(st)
  for(i in 1:length(st)){
    if(st[i] == TRUE){
      rect(
        xleft=n-1,
        xright=n,
        ybottom=i-1,
        ytop=i,
        col="black"
      )
    }
  }
}
shit <- apply(X=plotmat, MAR=2, FUN=plotcol)
n <- ncol(plotmat)
axis(side=1, at = (0:(n-1))+0.5, labels=1:n, cex.axis=0.7)
axis(side=2, at=(0:(n-1))+0.5, 
  labels=rownames(plotmat), 
  las=2, cex.axis=0.7)
abline(v=0:n, col="gray")
abline(h=0:n, col="gray")

dev.off()




