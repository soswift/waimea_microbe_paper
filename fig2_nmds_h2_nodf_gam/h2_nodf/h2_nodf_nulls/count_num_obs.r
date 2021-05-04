

	outputlist <- list()
	habsite_files <- paste0("../split_data/", list.files("../split_data/"))
	for(i in 1:length(habsite_files)){
		output_i <- list()
		# provenance
		output_i$filename <- habsite_files[i]
		output_i$habitat <- sub("-.*", "", sub(".*site_", "", habsite_files[i]))
		output_i$site <- sub("_empo.*", "", sub(".*-", "", habsite_files[i]))

		file_i <- habsite_files[i]
		mat_i <- read.table(file_i, header=T, sep=",")
		mat_i <- as.matrix(mat_i[,colnames(mat_i) != "OTU_ID"])
		output_i$nObs=sum(mat_i)
		output_i$nOTUs=nrow(mat_i)
		output_i$nCats=ncol(mat_i)

		outputlist[[i]] <- output_i

	}

	outdf <- do.call("rbind", outputlist)

	write.csv(x=outdf, file="matrix_dims_ns.csv", row.names=FALSE, quote=FALSE)