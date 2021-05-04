#!/usr/bin/env Rscript


# This script uses the Vazquez null model for bipartite networks to make null matrices
# from the habitat-site combination matrices already subset in ../split_data. 


# load packages
	library(Rcpp)
	sourceCpp("nodf.cpp")
	sourceCpp("vazfill.cpp")
	source("vaznull_faster.r")
	library(parallel)

	message("env loaded and compiled.")

# global parameters
	g_ncores <- 20
	g_nsim <- 300

# files to use
	habsite_files <- paste0("../split_data/", list.files("../split_data/"))

# quick function for H2. simplified version of bipartite::H2fun
	plogp <- function(p){
		output <- rep(0, length(p))
		output[p > 0] <- p[p>0] * log(p[p>0])
		return(output)
	}
	H2simple <- function(web){
		tot <- sum(web)
		rs <- rowSums(web)
		cs <- colSums(web)
		return( -sum(plogp(web/tot)) )
	}

	message("H2 funs defined.")


# do nestedness for each file. store results as list
	outputs_list <- list()
	for(i in 1:length(habsite_files)){

		output_i <- list()
		# provenance
		output_i$filename <- habsite_files[i]
		output_i$habitat <- sub("-.*", "", sub(".*site_", "", habsite_files[i]))
		output_i$site <- sub("_empo.*", "", sub(".*-", "", habsite_files[i]))

		message(paste0(i, "/",length(habsite_files),"; ", 
			output_i$habitat, "; ", output_i$site))

		# format matrix
		file_i <- habsite_files[i]
		mat_i <- read.table(file_i, header=T, sep=",")
		mat_i <- as.matrix(mat_i[,colnames(mat_i) != "OTU_ID"])
		message("...mat setup")


		# empirical values
		output_i$H2_emp <- H2simple(mat_i)
		output_i$NODF_emp <- NODFcpp(mat_i>0)
		message("...emp H2+NODF calculated")


		# make g_nsim nulls
		nulls_i <- mclapply(
			X=1:g_nsim, 
			FUN=function(x){vaznull_faster(mat_i, 200, noisy=FALSE)}, 
			mc.cores=g_ncores
		)
		message("...made sims")

		# calculate H2sim for each one
		output_i$H2_sim <- simplify2array(mclapply(
			X=nulls_i,
			FUN=H2simple,
			mc.cores=g_ncores
		))
		message("...sim H2 calculated")

		# calculate NODFsim for each one
		output_i$NODF_sim <- simplify2array(mclapply(
			X=nulls_i,
			FUN=function(mat){NODFcpp(mat>0)},
			mc.cores=g_ncores
		))
		message("...sim NODF calculated")


		outputs_list[[i]] <- output_i

	}

	saveRDS(outputs_list, file="H2NODF.RDS")
	message("Saved (all done)")
