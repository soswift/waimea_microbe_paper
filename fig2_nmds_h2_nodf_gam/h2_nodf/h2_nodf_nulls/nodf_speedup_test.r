# Testing how much faster NODFcpp is

# load packages
	library(Rcpp)
	sourceCpp("nodf.cpp")
	sourceCpp("vazfill.cpp")
	source("vaznull_faster.r")
	library(RInSp)


# files to use
	habsite_files <- paste0("../split_data/", list.files("../split_data/"))


	i <- 1
	file_i <- habsite_files[i]
	mat_i <- read.table(file_i, header=T, sep=",")
	mat_i <- as.matrix(mat_i[,colnames(mat_i) != "OTU_ID"])
	mat_i <- mat_i > 0
	mat_i_RInSp <- import.RInSp(mat_i)

	# test RInSp::NODF
		ptm <- proc.time()
		test <- RInSp::NODF(mat_i_RInSp)
		proc.time() - ptm

	# test my cpp version
		ptm <- proc.time()
		test <- NODFcpp(mat_i)
		proc.time() - ptm


	# test my vaznull implementation
		ptm <- proc.time()
		test <- vaznull_faster(mat_i, 300)
		proc.time() - ptm

		ptm <- proc.time()
		test <- vaznull(1, mat_i)
		proc.time() - ptm


