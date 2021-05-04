#!/usr/bin/env Rscript


# greet
  message("it's started!")

# load stuff
  library(data.table)
  library(parallel)
  # for rarefaction:
  source("rare_funs_endos.r")
  # for gen alg:
  source("genOptRich.r")



# read in data
  otutable_fp <- "../../data/waimea_full/abundance_table_100.shared"
  metadata_fp <- "../../data/waimea_full/combined_prep_sample_metadata_16s.csv"
  taxonomy_fp <- "../../data/waimea_full/annotations_100.taxonomy"
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
  # remove samples with ultra low sequencing depth (<6000 reads)
  otutable <- otutable[rowSums(otutable) >= 6000, ]
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


# make column of combined habitat + empo_3 categories
  metadata$habempo3 <- paste(metadata$habitat, metadata$empo_3)

# check numbers of habempo3s
  table(metadata$habempo3)
  # min number is 11. seems fine. go down to 10 for a nice even number.

# collapse otutable by habempo3 categories
  agg_by_cat <- function(x, cats=metadata$habempo3){
    cols <- lapply(
      X=unique(cats),
      FUN=function(cat){
        rowSums(x[,cats==cat,drop=F])
      }
    )
    out <- do.call("cbind", cols)
    colnames(out) <- unique(cats)
    return(out)
  }

# wrapper function to
  # a) get 10 samples from each habempo3
  # b) rarefy to 60seqs/sample
  # c) aggregate by habempo3
  # d) run gen alg
  gen_alg_rare_wrap <- function(seed, otutable2=otutable, categories=metadata$habempo3,
   nsam=10, raredep=6000){
    # a) get 10 samples from each habempo3
    set.seed(seed)
    sampstoget <- as.vector(sapply(X=unique(categories), 
      FUN=function(x){sample(which(categories==x), size=nsam)}))
    otb <- otutable2[ , sampstoget]
    otb_cats <- categories[sampstoget]
    # b) rarefy to 60seqs/sample
    otb <- rarefy_otu_table_cpp(otb, raredep, seed)
    otb <- otb[rowSums(otb)>0,]
    # c) aggregate by habempo3
    otb <- agg_by_cat(otb, cats=otb_cats)

    results_list <- list()
    for(i in 1:ncol(otb)){
      results_list[[i]] <- genOptRich(
        otutable = t(otb),
        popsize = 300,
        best2keep = 20,
        regs2keep = 20,
        gens2stop = 10,
        maxgens = 300,
        mutrate = 3,
        noisy=F,
        n=i
      )
    }
    message("one iter done")
    richs <- colSums(otb > 0)
    catnames <- colnames(otb)
    return(list(res_list=results_list, richs=richs, cats=catnames))
  }

# run gen_alg_rare_wrap() 500 times in parallel
  seeds <- 12345 + 1:500
  perm_results <- mclapply(X=seeds, FUN=gen_alg_rare_wrap, mc.cores=20)

# save
  save(list="perm_results", file="perm_results.rdata")
  # load(file="perm_results.rdata")

# check that cats are in same order for all perms (should be!)
  cats2check <- sapply(X=perm_results, FUN=function(x){x$cats})
  # next line should be TRUE:
  all(apply(X=cats2check, MAR=1, FUN=function(x){all(x==x[1])}))

# prepare data for plotting
  # sum results into matrix of proportions by column
  get_mat_from_rep <- function(a){
    out <- sapply(X=a$res_list, FUN=function(x){x$bestState})
    rownames(out) <- a$cats
    return(out)
  }
  mats <- lapply(X=perm_results, FUN=get_mat_from_rep)
  plotmat <- Reduce("+", mats)/length(mats)

  # create accumulation curve
  # each column of rs corresponds to ns
  ns <- sapply(X=perm_results[[1]]$res_list, FUN=function(x){sum(x$bestState)})
  rs <- sapply(X=perm_results, FUN=function(a){sapply(X=a$res_list, FUN=function(x){max(x$richnessByGen)})})
  rs_mean <- apply(X=rs, MAR=1, FUN=mean)
  rs_sd <- apply(X=rs, MAR=1, FUN=sd)

# order plotmat
  plotmat <- plotmat[order(rowSums(plotmat)), ]

# change plotmat labels to be prettier
  rownames(plotmat) <- sub(pattern="non-saline", replacement="ns", x=rownames(plotmat))
  rownames(plotmat) <- sub(pattern="saline", replacement="s", x=rownames(plotmat))

# function to get trophic from empo3
  colors_table <- read.table("../../colorstable.txt", comment.char="", stringsAsFactors=FALSE, sep="\t", header=T)
  empo3_habtro_map <- read.table("empo32trophic.txt", stringsAsFactors=FALSE, header=T, sep="\t")
  empo32col <- function(empo3, type="black"){
    if(type %in% c("trophic", "tro")){
      tro <- empo3_habtro_map$trophic[empo3_habtro_map$empo3 == empo3]
      return(colors_table$trocols[colors_table$tros == tro])
    }else if(type %in% c("habitat", "hab")){
      hab <- empo3_habtro_map$habitat[empo3_habtro_map$empo3 == empo3]
      return(colors_table$habcols[colors_table$habs == hab])
    }else{
      return(type)
     }

  }

# do plotting
  pdf("combined_plot_habitat.pdf")
  
  coloring <- "habitat" # can also be "black" or "red" or "habitat" or any other color

  par(mfrow=c(2,1))

  par(mar=c(0.5,9,3,1))

  # make empty plot (top)
  plot(0, type="n", ylab="", xaxt="n", yaxt="n", xaxs="i", 
    xlim=c(1, nrow(plotmat)+1), ylim=range(rs), xlab="")
  # bin lines
  abline(v=2:nrow(plotmat), col="gray")
  # grey lines
  apply(X=rs, MAR=2, FUN=function(y){points(ns+0.5, y, col="gray", type="l")})
  # apply(X=rs[,1:10], MAR=2, FUN=function(y){points(ns+0.5, y, col="gray", type="l")})
  # mean
  points(ns+0.5, rs_mean, pch=20, cex=3, type="b", lwd=2)
  # labels
  mtext(side=2, text="Total maximized ASV richness", line=5)
  yaxvals <- 2:8 * 10000
  yaxlabs <- as.character(yaxvals)
  axis(side=2, at=yaxvals, labels=parse(text=yaxlabs), cex.axis=0.7, las=2)

  # make empty plot (bottom)
  par(mar=c(4,9,0,1))
  plot(0, type="n", 
    xlim=c(0, nrow(plotmat)),
    ylim=c(0, nrow(plotmat)),
    xlab="Number of categories",
    ylab="", yaxt="n", xaxt="n",
    xaxs="i", yaxs="i"
  )
  # rowcolors
  rowcolors <- sapply(X=rownames(plotmat), FUN=empo32col, type=coloring)
  # quick function to make a color transparent
  col2trans <- function(col, alpha){
    rgb(
      red=col2rgb(col)[1]/255,
      green=col2rgb(col)[2]/255,
      blue=col2rgb(col)[3]/255,
      alpha
    )
  }

  plotcol <- function(st, n){
    for(i in 1:length(st)){
      col_i <- col2trans(rowcolors[i], st[i])
      # col_i <- rgb(0,0,0,st[i])
      rect(
        xleft=n-1,
        xright=n,
        ybottom=i-1,
        ytop=i,
        col=col_i,
        border=NA
      )
    }
  }
  for(i in 1:ncol(plotmat)){plotcol(plotmat[,i], i)}

  n <- ncol(plotmat)
  axis(side=1, at = (0:(n-1))+0.5, labels=1:n, cex.axis=0.7)
  axis(side=2, at=(0:(n-1))+0.5, 
    labels=rownames(plotmat), 
    las=2, cex.axis=0.7)
  abline(v=0:n, col="gray")
  abline(h=0:n, col="gray")

dev.off()


# plot of richness by empo3
  # get richnesses
  richmat <- sapply(X=perm_results, FUN=function(x){x$richs})
  # change labels as per above
  rownames(richmat) <- sub(pattern="non-saline", replacement="ns", x=rownames(richmat))
  rownames(richmat) <- sub(pattern="saline", replacement="s", x=rownames(richmat))
  # sort richmat to match plotmat
  richmat <- richmat[ order(match(rownames(richmat), rownames(plotmat))) , ]
  # check
  all(rownames(richmat) == rownames(plotmat))
  # divide by 10 (its richness of 10 samples) to get per-sample richness
  richmat <- richmat / 10
  # aggregate into quartiles
  rich_quarts <- as.data.frame(t(apply(X=richmat, MARGIN=1, FUN=quantile)))
  colnames(rich_quarts) <- c("min", "lq", "median", "uq", "max")

  # just make plot same as before
  pdf("richness_by_category.pdf", useDingbats = FALSE)
  par(mfrow=c(2,1))

  par(mar=c(0.5,0,3,25))
  plot(0, type="n", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")

  par(mar=c(4,1,0,25))

  plot(0, type="n", yaxt="n", ylab="",
    xlim=range(c(rich_quarts$min, rich_quarts$max)), ylim=c(0, nrow(richmat)),
    xlab="Sample ASV richness", xaxs="i", yaxs="i", bty="n", cex.axis=0.7)

  # draw whiskers
  segments(x0=rich_quarts$min, x1=rich_quarts$max, 
    y0=(1:nrow(rich_quarts))-0.5, y1=(1:nrow(rich_quarts))-0.5,
    col="black")

  # draw rectangles
  rect(xleft=rich_quarts$lq, xright = rich_quarts$uq, 
    ybottom=(0:(nrow(rich_quarts)- 1)), 
    ytop=(1:nrow(rich_quarts)), border="white", col="black")

  # draw median lines
  # segments(x0=rich_quarts$median, x1=rich_quarts$median, 
  #   y0=(0:(nrow(rich_quarts)- 1)), y1=(1:nrow(rich_quarts)), 
  #   col="white", lwd=1)


  points(x=rich_quarts$median, y=(1:nrow(rich_quarts))-0.5,
    pch=18, col="white", cex=0.7)

  # gridlines
  abline(h=1:nrow(rich_quarts), col="gray")

  # write out
  dev.off()




# partitioning alpha diversity - nested vs unique
  r_max <- sapply(X=perm_results, FUN=function(a){sapply(X=a$res_list, FUN=function(x){max(x$richnessByGen)})})
  r_cum <- matrix(0, nrow=nrow(r_max), ncol=ncol(r_max))
  for(i in 1:nrow(r_cum)){
    for(j in 1:ncol(r_cum)){
      a <- perm_results[[j]]
      r_cum[i,j] <- sum( a$richs[ a$res_list[[i]]$bestState ] )
    }
    
  }



plot(r_max, r_cum, pch=20, cex=0.1)
abline(a=0, b=1)

nsampmat <- matrix(rep(1:nrow(rs), ncol(rs)), ncol=ncol(rs), byrow=F)

plot(r_cum ~ nsampmat)
points(r_max ~ nsampmat, col="red")