  load(file="perm_results.rdata")

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

# create accumulation curve, split by total richness and nestedness component,
  # i.e. "unique to 1 category", "shared between 2 categories", etc
  # done by taking average of all the "nest" matrices:
    # rows: nestedness/sharedness. 0 means "number of otus found in 0 categories",
    # 1 means "number of otus found in only 1 category", et cetera
    # cols: same as plotmat. number of categories algo was told to choose
  nests_list <- lapply(X=perm_results, FUN=function(x){x$nest})
  nests_means <- Reduce("+", nests_list)/length(nests_list)
  # get rid of zeroes
  nests_means <- nests_means[rownames(nests_means) > 0, ]
  nests_means_cum <- apply(X=nests_means, FUN=function(x){rev(cumsum(rev(x)))}, MARGIN=2)

# create accumulation curve totals
  # each column of total_rs corresponds to ns
  ns <- sapply(X=perm_results[[1]]$res_list, FUN=function(x){sum(x$bestState)})
  total_rs <- sapply(X=perm_results, FUN=function(a){sapply(X=a$res_list, FUN=function(x){max(x$richnessByGen)})})
  total_rs_mean <- apply(X=total_rs, MAR=1, FUN=mean)
  total_rs_sd <- apply(X=total_rs, MAR=1, FUN=sd)

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
  
  coloring <- "trophic" # can also be "black" or "red" or "habitat" or any other color

  par(mfrow=c(2,1))

  par(mar=c(0.5,9,3,1))

  # make empty plot (top)
  plot(0, type="n", ylab="", xaxt="n", yaxt="n", xaxs="i", yaxs="i", 
    xlim=c(1, nrow(plotmat)+1), ylim=c(0, max(total_rs)), xlab="")

  # polygons for nestedness levels
  # polycols <- c( '#fa9fb5','#f768a1','#dd3497','#ae017e' ) # ugly flesh
  # polycols <- c("#fed98e", "#fe9929", "#d95f0e", "#993404") # aggro orange
  # polycols <- (c("#b2e2e2", "#66c2a4", "#2ca25f", "#006d2c")) # mint, too similar to PP
   polycols <- rev(c("#bdc9e1", "#74a9cf", "#2b8cbe", "#045a8d"))
  # polycols <- c(
  # 	rgb(0.9, 0.9, 0.9, maxColorValue=1),
  # 	rgb(0.8, 0.8, 0.8, maxColorValue=1),
  # 	rgb(0.6, 0.6, 0.6, maxColorValue=1),
  # 	rgb(0.4, 0.4, 0.4, maxColorValue=1)
  # )
  for(n in 1:4){
    polygon(
      x=c(0.5, ns+0.5,              max(ns)+1,                max(ns)+1 ),
      y=c(0,   nests_means_cum[n,], max(nests_means_cum[n,]), 0         ),
      col=polycols[n], border="black", lwd=2)
  }


  # bin lines
  abline(v=2:nrow(plotmat), col="gray")
  # grey lines
  apply(X=total_rs, MAR=2, FUN=function(y){points(ns+0.5, y, col="darkgray", type="l")})
  # apply(X=total_rs[,1:10], MAR=2, FUN=function(y){points(ns+0.5, y, col="gray", type="l")})
  # mean
  points(ns+0.5, total_rs_mean, pch=20, cex=3, type="b", lwd=2)
  # labels
  mtext(side=2, text="Total maximized ASV richness", line=5)
  yaxvals <- 0:8 * 10000
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
    ytop=(1:nrow(rich_quarts)), border="black", col=rowcolors)

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






