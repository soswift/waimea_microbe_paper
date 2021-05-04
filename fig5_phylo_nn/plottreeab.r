

library(ape)

# quick function to make a color transparent
col2trans <- function(col, alpha){
	rgb(
		red=col2rgb(col)[1]/255,
		green=col2rgb(col)[2]/255,
		blue=col2rgb(col)[3]/255,
		alpha
	)
}


# data:
#  column "label" which corresponds to tip labels on tree
#  column "x" for the variable to be plotted
#  column "g" for groups to be plotted; must have 2 levels
#  type:
#    quartiles = simple boxplot just showing median and quartiles
#      uses args: cola, colb, colbg
#    dots = dots for medians, connected by line
#      uses args: cola, colb, colbg, linecol, linelwd, pointcex
#    overlap = quantiles but overlapping and transparent
#      uses args: cola, colb, trans, pointcex

plottreeab <- function(tree, data, xlabel="", cola="red", colb="blue", 
	colbg="gray98", linecol="gray", linelwd=2, pointcex=2, trans=0.20, 
	type="quartiles", xlab="distance"){

	# data <- plotdata_less
	# tree <- ggtree_class
	# cola <- "orange"
	# colb <- "blue"
	# colbg <- "gray98"
	# linecol="gray"
	# linelwd=2
	# pointcex=2
	# trans=0.20
	# type="violins"
	# xlabel="distance"
	plotags_violins <- function(tip){
		# tip <- tree$tip.label[1]
		groupa_name <- as.character(sort(unique(data$g))[1])
		groupb_name <- as.character(sort(unique(data$g))[2])
		# y = height for plot, middle
		y <- which(tree$tip.label == tip)
		# calculate densities
		dens_a <- density(data$x[data$label==tip & data$g == groupa_name])
		y <- which(tree$tip.label == tip)
		dens_b <- density(data$x[data$label==tip & data$g == groupb_name])
		y <- which(tree$tip.label == tip)
		# rescale so max dens is 0.5
		dens_a$y <- (dens_a$y / max(dens_a$y))/2
		dens_b$y <- (dens_a$y / max(dens_a$y))/2
		# draw polygons
		polygon(
			x=c(dens_a$x, rev(dens_a$x)),
			y=c(y+dens_a$y, y-rev(dens_a$y)),
			col=col2trans(cola, trans), border=NA, lwd=linelwd
		)
		polygon(
			x=c(dens_b$x, rev(dens_b$x)),
			y=c(y+dens_b$y, y-rev(dens_b$y)),
			col=col2trans(colb, trans), border=NA, lwd=linelwd
		)
		polygon(
			x=c(dens_a$x, rev(dens_a$x)),
			y=c(y+dens_a$y, y-rev(dens_a$y)),
			col=NA, lwd=linelwd
		)
		polygon(
			x=c(dens_b$x, rev(dens_b$x)),
			y=c(y+dens_b$y, y-rev(dens_b$y)),
			col=NA, lwd=linelwd
		)
	}

	plotags_joy <- function(tip){
		# tip <- tree$tip.label[1]
		groupa_name <- as.character(sort(unique(data$g))[1])
		groupb_name <- as.character(sort(unique(data$g))[2])
		# y = height for plot, middle
		y <- which(tree$tip.label == tip)
		# calculate densities
		dens_a <- density(data$x[data$label==tip & data$g == groupa_name])
		y <- which(tree$tip.label == tip)
		dens_b <- density(data$x[data$label==tip & data$g == groupb_name])
		y <- which(tree$tip.label == tip)
		# rescale so max dens is 1.5 (0.5 overlap)
		dens_a$y <- (dens_a$y / max(dens_a$y)) * 1.25
		dens_b$y <- (dens_a$y / max(dens_a$y)) * 1.25

		# draw polygons
		polygon(
			x=dens_a$x,
			y=y+dens_a$y-0.5,
			col="white", border=NA, lwd=linelwd
		)
		polygon(
			x=dens_b$x,
			y=y+dens_b$y-0.5,
			col="white", border=NA, lwd=linelwd
		)

		points(x=dens_a$x, y=y+dens_a$y-0.5, type="l", col=cola, lwd=linelwd)
		points(x=dens_b$x, y=y+dens_b$y-0.5, type="l", col=colb, lwd=linelwd)

		segments(x0 = min(data$x), x1=max(data$x), y0=y-0.5, col="black", lwd=linelwd)

	}


	plotags_quartiles <- function(tip, ret="plot"){
		# tip <- tree$tip.label[1]
		y <- which(tree$tip.label == tip)

		groupa_name <- as.character(sort(unique(data$g))[1])
		groupb_name <- as.character(sort(unique(data$g))[2])

		aq <- quantile( data$x[data$label==tip & data$g == groupa_name], probs=c(0.05, 0.25, 0.50, 0.75, 0.975))
		bq <- quantile( data$x[data$label==tip & data$g == groupb_name], probs=c(0.05, 0.25, 0.50, 0.75, 0.975))
		# am <- mean ( data$x[data$label==tip & data$g == groupa_name] )
		# bm <- mean ( data$x[data$label==tip & data$g == groupb_name] )

		if(ret=="lims"){
			alims <- c(min(aq[1], bq[1]), max(aq[5], bq[5]))
		}else{
			segments(x0=aq[1], x1=aq[5], y0=y+0.15, col=cola, lwd=linelwd)
			rect(
				xleft=aq[2], xright=aq[4],
				ybottom=y, ytop=y+0.3, col="white", border=NA
			)
			rect(
				xleft=aq[2], xright=aq[4],
				ybottom=y, ytop=y+0.3, col=col2trans(cola, trans), border=cola
			)
			segments(x0=aq[3], y0=y, y1=y+0.30, col=cola, lwd=linelwd)
			#points(x=aq[3], y=y+0.15, pch=20, col="white")
			# points(x=am, y=y+0.15, pch=20, col="white")

			segments(x0=bq[1], x1=bq[5], y0=y-0.15, col=colb, lwd=linelwd)
			rect(
				xleft=bq[2], xright=bq[4],
				ybottom=y-0.3, ytop=y, col="white", border=NA
			)
			rect(
				xleft=bq[2], xright=bq[4],
				ybottom=y-0.3, ytop=y, col=col2trans(colb, trans), border=colb
			)
			segments(x0=bq[3], y0=y, y1=y+-0.30, col=colb, lwd=linelwd)
			# points(x=bq[3], y=y-0.15, pch=20, col="white")
			# points(x=bm, y=y+0.15, pch=20, col="white")
		}
	}

	plotags_quartiles2 <- function(tip, ret="plot"){
		# tip <- tree$tip.label[1]
		y <- which(tree$tip.label == tip)

		groupa_name <- as.character(sort(unique(data$g))[1])
		groupb_name <- as.character(sort(unique(data$g))[2])

		aq <- quantile( data$x[data$label==tip & data$g == groupa_name], probs=c(0.05, 0.25, 0.50, 0.75, 0.975))
		bq <- quantile( data$x[data$label==tip & data$g == groupb_name], probs=c(0.05, 0.25, 0.50, 0.75, 0.975))
		# am <- mean ( data$x[data$label==tip & data$g == groupa_name] )
		# bm <- mean ( data$x[data$label==tip & data$g == groupb_name] )

		if(ret=="lims"){
			alims <- c(min(aq[2], bq[2]), max(aq[4], bq[4]))
		}else{
			#segments(x0=aq[1], x1=aq[5], y0=y+0.15, col=cola, lwd=linelwd)
			rect(
				xleft=aq[2], xright=aq[4],
				ybottom=y, ytop=y+0.3, col="white", border=NA
			)
			rect(
				xleft=aq[2], xright=aq[4],
				ybottom=y, ytop=y+0.3, col=col2trans(cola, trans), border=cola
			)
			segments(x0=aq[3], y0=y, y1=y+0.30, col=cola, lwd=linelwd)
			#points(x=aq[3], y=y+0.15, pch=20, col="white")
			# points(x=am, y=y+0.15, pch=20, col="white")

			#segments(x0=bq[1], x1=bq[5], y0=y-0.15, col=colb, lwd=linelwd)
			rect(
				xleft=bq[2], xright=bq[4],
				ybottom=y-0.3, ytop=y, col="white", border=NA
			)
			rect(
				xleft=bq[2], xright=bq[4],
				ybottom=y-0.3, ytop=y, col=col2trans(colb, trans), border=colb
			)
			segments(x0=bq[3], y0=y, y1=y+-0.30, col=colb, lwd=linelwd)
			# points(x=bq[3], y=y-0.15, pch=20, col="white")
			# points(x=bm, y=y+0.15, pch=20, col="white")
		}


	}

	plotags_ellipses <- function(tip){
		require(plotrix)
		# tip <- tree$tip.label[1]
		y <- which(tree$tip.label == tip)

		groupa_name <- as.character(sort(unique(data$g))[1])
		groupb_name <- as.character(sort(unique(data$g))[2])
		ax <- data$x[data$label==tip & data$g == groupa_name]
		bx <- data$x[data$label==tip & data$g == groupb_name]
		am <- mean ( ax )
		bm <- mean ( bx )
		asd <- sd(ax)
		bsd <- sd(bx)

		draw.ellipse(x=am, y=y, a=asd, b=0.5, border=cola, lwd=linelwd)
		draw.ellipse(x=bm, y=y, a=bsd, b=0.5, border=colb, lwd=linelwd)
		points(x=am, y=y, col=cola, pch=15)
		points(x=bm, y=y, col=colb, pch=16)
	}


	plotags_quartilesold <- function(tip){
		y <- which(tree$tip.label == tip)
		# box for "a"
		rect(
			xleft=agdata$lq[agdata$label==tip & agdata$groupab=="a"],
			xright=agdata$uq[agdata$label==tip & agdata$groupab=="a"],
			ybottom=y, ytop=y+0.25, col=cola, border=NA
		)
		rect(
			xleft=agdata$lq[agdata$label==tip & agdata$groupab=="b"],
			xright=agdata$uq[agdata$label==tip & agdata$groupab=="b"],
			ybottom=y, ytop=y-0.25, col=colb, border=NA
		)
		segments(
			x0=agdata$x[agdata$label==tip & agdata$groupab=="a"],
			y0=y, y1=y+0.25, lwd=3, col="black"
		)
		segments(
			x0=agdata$x[agdata$label==tip & agdata$groupab=="b"],
			y0=y, y1=y-0.25, lwd=3, col="black"
		)
	}


	plotags_dots <- function(tip){
		y <- which(tree$tip.label == tip)
		# line from a to b
		segments(
			x0=agdata$x[agdata$label==tip & agdata$groupab=="a"],
			x1=agdata$x[agdata$label==tip & agdata$groupab=="b"],
			y0=y, y1=y, lwd=linelwd, col=linecol
		)
		points(
			x=agdata$x[agdata$label==tip & agdata$groupab=="a"],
			y=y, pch=15, col=cola, cex=pointcex
		)
		points(
			x=agdata$x[agdata$label==tip & agdata$groupab=="b"],
			y=y, pch=16, col=colb, cex=pointcex
		)
	}

	plotags_overlap <- function(tip){
		y <- which(tree$tip.label == tip)
		# box for "a"
		rect(
			xleft=agdata$lq[agdata$label==tip & agdata$groupab=="a"],
			xright=agdata$uq[agdata$label==tip & agdata$groupab=="a"],
			ybottom=y-0.5, ytop=y+0.5, col=col2trans(cola, trans), border=NA
		)
		rect(
			xleft=agdata$lq[agdata$label==tip & agdata$groupab=="b"],
			xright=agdata$uq[agdata$label==tip & agdata$groupab=="b"],
			ybottom=y-0.5, ytop=y+0.5, col=col2trans(colb, trans), border=NA
		)
		points(
			x=agdata$x[agdata$label==tip & agdata$groupab=="a"],
			y=y, pch=15, col=cola, cex=pointcex
		)
		points(
			x=agdata$x[agdata$label==tip & agdata$groupab=="b"],
			y=y, pch=16, col=colb, cex=pointcex
		)
	}

	ntip <- length(tree$tip.label)
	par(mfrow=c(1,2))
	# plot tree
	par(mar=c(4,0,0,0))
	plot(tree, adj=1, align.tip.label=TRUE, edge.width=2)

	# scale bar
	add.scale.bar(x=2, y=length(tree$tip.label))


	# aggregate data, calculate xlims, and plot
	agdata <- aggregate(x=data$x, by=list(label=data$label, group=data$g), FUN=median)
	agdata$lq <- aggregate(x=data$x, by=list(label=data$label, group=data$g), FUN=function(x){quantile(x)[2]})$x
	agdata$uq <- aggregate(x=data$x, by=list(label=data$label, group=data$g), FUN=function(x){quantile(x)[4]})$x
	
	# change group to a and b
	lvls <- sort(unique(agdata$group))
	if(length(lvls) != 2){stop("only 2 group levels are supported at this time")}
	agdata$groupab <- "a"
	agdata$groupab[agdata$group == lvls[2]] <- "b"

	# plot more stuff
	par(mar=c(4,0,0,0))
	if(type %in% c("violins", "joy")) {
		xlims <- range(data$x)
	}else if(type == "quartiles"){
		xlimmat <- sapply(X=tree$tip.label, FUN=plotags_quartiles, ret="lims")
		xlims <- c(
			min(xlimmat[1,]),
			max(xlimmat[2,])
		)
	}else if(type == "quartiles2"){
		xlimmat <- sapply(X=tree$tip.label, FUN=plotags_quartiles2, ret="lims")
		xlims <- c(
			min(xlimmat[1,]),
			max(xlimmat[2,])
		)
	}else{
		xlims <- range(c(agdata$lq, agdata$uq))
	}
	plot(0, type="n", ylim=c(1,ntip), xlim=xlims, yaxt="n", ylab="", bty="n", xlab=xlabel) 

	# add background shading for odd tips
	plotshade <- function(ymid){rect(xleft=xlims[1], xright=xlims[2],
		ybottom=ymid-0.5, ytop=ymid+0.5, col=colbg, border=NA)}
	trash <- sapply(X=((1:ntip)[(1:ntip %% 2 == 1)]), FUN=plotshade)



	if(type=="quartiles"){
		trash <- sapply(X=tree$tip.label, FUN=plotags_quartiles)
	}else if(type=="quartiles2"){
		trash <- sapply(X=tree$tip.label, FUN=plotags_quartiles2)
	}else if(type == "dots"){
		trash <- sapply(X=tree$tip.label, FUN=plotags_dots)
	}else if(type=="overlap"){
		trash <- sapply(X=tree$tip.label, FUN=plotags_overlap)
	}else if(type=="violins"){
		trash <- sapply(X=tree$tip.label, FUN=plotags_violins)
	}else if(type=="joy"){
		trash <- sapply(X=rev(tree$tip.label), FUN=plotags_joy)
	}else if(type=="ovals"){
		trash <- sapply(X=rev(tree$tip.label), FUN=plotags_ellipses)
	}else{
		stop(paste0("plot type \"", type, "\" not defined."))
	}


	# for(i in 1:length(tree$tip.label)){
	# 	print(i)
	# 	plotags(tree$tip.label[i])
	# }



}