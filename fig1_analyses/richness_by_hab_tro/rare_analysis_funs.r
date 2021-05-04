require(parallel)



# function to permute column order
	# only useful for "sample" style rarefaction
	perm_col_order <- function(otb){
		otb[,sample(1:ncol(otb))]
	}

# function to make alpha collectors curve by samples in table otb
	# note these samples could be real samples OR rarefactions
	rare_perm_curve_samp <- function(otb){
		# permute column order then accumulate table
		otb <- otb[, sample(1:ncol(otb))]
		# do cumulative sum
		otb_cs <- t(apply(X=otb, FUN=cumsum, MAR=1))
		# calculate richness
		rich <- apply(X=otb_cs, MAR=2, FUN=function(x){sum(x>0)})
		nobs <- apply(X=otb_cs, MAR=2, FUN=sum)
		nsam <- 1:ncol(otb_cs)
		return(data.frame(rich, nobs, nsam, row.names=NULL))
	}

# function to make perm curves
	make_perm_curves <- function(otb, nperm=500, ncores=20){
		mclapply(
			X=1:nperm,
			FUN=function(x){rare_perm_curve_samp(perm_col_order(otb))},
			mc.cores=ncores
		)
	}

# function to make rbrtabs
	# rbrtabs = "rbind rtabs"
	# required because gamlss can only look in global namespace for inputs for
	# some IDIOTIC reason.
	make_rbrtabs <- function(rtabs){
		# transform rtabs into a giant dataframe via rbind
		rbrtabs <- do.call("rbind", rtabs)
		rbrtabs$weight <- log(rbrtabs$nsam)
		rbrtabs$rep <- factor(rep(1:length(rtabs), sapply(X=rtabs, FUN=nrow)))

		# log variables for gamlss
		rbrtabs$lnrich <- log(rbrtabs$rich)
		rbrtabs$lnnsam <- log(rbrtabs$nsam)
	
		# remove first sample since its bogus
		rbrtabs <- rbrtabs[rbrtabs$nsam > 1,]

		return(rbrtabs)
	}

# function to get 95% confidence interval log models
	# doing it just with the 95CIs from lm()
	model_curves_gamlss <- function(rbrtabs){
		require(gamlss)

        # set controll settings for gamlss
        i_control_sett <- glim.control()
        i_control_sett$"c.cyc" <- 300
        i_control_sett$"bf.cyc" <- 300
        g_control_sett <- gamlss.control()
        g_control_sett$"n.cyc" <- 300

		m1_gam <- gamlss(
			formula=lnrich~lnnsam, #+random(rep),
			sigma.formula=~lnnsam,
			weights=weight,
			family="NO",
			data=rbrtabs,
			i.control=i_control_sett,
			control=g_control_sett
		)

		coefs <- coef(m1_gam)
		cis <- confint(m1_gam)

		return(list(
			LCI_025=list(logc=as.numeric(cis[1,1]), z=as.numeric(cis[2,2])),
			HCI_975=list(logc=as.numeric(cis[1,2]), z=as.numeric(cis[2,1])),
			mean=list(logc=as.numeric(coefs[1]), z=as.numeric(coefs[2]))
		))

	}


# function to get richness from a model at some point x
	get_logmodel_rich <- function(mod, x){
		# points(x=x, y=exp( logc+z*log(x) ), type="l")
		exp( mod$logc+mod$z*log(x))
	}

# function to get slope from a model at some point x
	get_logmodel_slope <- function(mod, x){
		# points(x=x, y=exp( logc+z*log(x) ), type="l")
		exp(mod$logc) * mod$z * x^(mod$z - 1)
	}

# function to find nsamp where slope = m
	get_nsamp_from_slope <- function(mod, m){
		a <- mod$logc
		z <- mod$z
		((exp(-1*a)*m)/z)^(1/(z-1))
	}

# plot model
	plotmodel <- function(mod, xmax, xmin=1, nx=100, ...){
		xs <- seq(from=xmin, to=xmax, length.out=nx)
		ys <- sapply(X=xs, FUN=function(x){get_logmodel_rich(mod, x)})
		points(ys ~ xs, type="l", ...)
	}


# function to plot model fitting
	# rtabs - output from make_perm_curves()
	# mods - output from model_curves_gamlss()
	plot_model_fitting <- function(rbrtabs, mods){

		par(mfrow=c(2,1))

		# linear scale
		plot(rbrtabs$lnrich ~ rbrtabs$lnnsam, pch=20, cex=0.5, col="gray",
			xlab="ln Number of Samples", ylab="ln Richness")
		abline(a=mods$mean$logc, b=mods$mean$z, col="blue", lwd=2)
		abline(a=mods$LCI_025$logc, b=mods$LCI_025$z, col="blue", lty=2)
		abline(a=mods$HCI_975$logc, b=mods$HCI_975$z, col="blue", lty=2)


		# log scale
		plot(rbrtabs$rich ~ rbrtabs$nsam, pch=20, cex=0.5, col="gray", 
			xlab="Number of Samples", ylab="Richness")
		unique_sams <- sort(unique(rbrtabs$nsam))
		pred_y_linear_mean <- sapply(X=unique_sams, FUN=function(x){get_logmodel_rich(mods$mean, x)})
		pred_y_linear_lci <- sapply(X=unique_sams, FUN=function(x){get_logmodel_rich(mods$LCI_025, x)})
		pred_y_linear_hci <- sapply(X=unique_sams, FUN=function(x){get_logmodel_rich(mods$HCI_975, x)})
		points(x=unique_sams, y=pred_y_linear_mean, type="l", col="blue", lwd=2)
		points(x=unique_sams, y=pred_y_linear_lci, type="l", col="blue", lty=2)
		points(x=unique_sams, y=pred_y_linear_hci, type="l", col="blue", lty=2)
	}



# function to make curves for a given hab+tro
	add_info_to_model <- function(rbrtabs, logmods, m){
		
		nsamp_slopes_m <- sapply(X=logmods, FUN=get_nsamp_from_slope, m=m)

		exp_rich_m <- sapply(X=1:length(logmods), FUN=function(x){
			get_logmodel_rich(logmods[[x]], nsamp_slopes_m[x])
		})

		obs_max_rich <- max(rbrtabs$rich)

		obs_nsamp=max(rbrtabs$nsam)
		pred_rich_1samp=get_logmodel_rich(mod=logmods$mean, x=1)

		return(list(
			coefs=logmods,
			nsamps_m=nsamp_slopes_m,
			exp_rich_m=exp_rich_m,
			obs_max_rich=obs_max_rich,
			obs_nsamp=obs_nsamp,
			pred_rich_1samp=pred_rich_1samp
		))
	}


	# fl - a list of outputs from add_into_to_model
	plot_fits_list <- function(fl, cols, pchs, xlim=NULL, ylim=NULL, pmain=NULL){
		# make empty plot
		ymin <- min(sapply(X=fl, FUN=function(x){x$pred_rich_1samp}))
		ymax <- max(sapply(X=fl, FUN=function(x){max(x$exp_rich_m)}))
		xmax <- max(sapply(X=fl, FUN=function(x){max(x$nsamps_m)}))
		xmin <- 1
		if(is.null(xlim)){xlim <- c(xmin, xmax)}
		if(is.null(ylim)){ylim <- c(ymin, ymax)}
		plot(0, xlab="Number of Samples", ylab="Richness", type="n",
			xlim=xlim, ylim=ylim, main=pmain)

		# function to plot a fit in fl
		plot_fit <- function(fit, fcol, fpch){
			# fit <- fl[[1]]
			# fcol <- "red"
			# fpch <- 20
			xmax <- round(max(fit$nsamps_m))
			# curve for mean
			plotmodel(fit$coefs[[3]], xmax, col=fcol, lty=1)
			# curves for LCI+HCI
			plotmodel(fit$coefs[[1]], xmax, col=fcol, lty=2)
			plotmodel(fit$coefs[[2]], xmax, col=fcol, lty=2)
			# point for where actual sampling ends
			pcy <- fit$obs_max_rich
			pcx <- fit$obs_nsamp
			points(pcy~pcx, pch=fpch, col=fcol, cex=2)
		}

		for(i in 1:length(fl)){
			plot_fit(fl[[i]], cols[i], pchs[i])
		}

	}


