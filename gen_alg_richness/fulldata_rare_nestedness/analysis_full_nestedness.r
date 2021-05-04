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
  otutable_fp <- "../../data/raw/sequence/waimea_16S_v1/abundance_table_100.shared"
  metadata_fp <- "../../data/raw/sequence/waimea_16S_v1/combined_prep_sample_metadata_16s.csv"
  taxonomy_fp <- "../../data/raw/sequence/waimea_16S_v1/annotations_100.taxonomy"
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

    # d) nestedness stuff
    # note: cols correspond to columns of otb, rows are counts of species
    # shared with n number of samples, where n is rowname
    otu_nestedness <- sapply(X=results_list, FUN=function(x){
      table(factor( rowSums( otb[, x$bestState, drop=F] > 0 ), levels=0:ncol(otb)))
    })


    message("one iter done")
    richs <- colSums(otb > 0)
    catnames <- colnames(otb)
    return(list(res_list=results_list, richs=richs, cats=catnames, nest=otu_nestedness))
  }

# run gen_alg_rare_wrap() 500 times in parallel
  seeds <- 12345 + 1:500
  perm_results <- mclapply(X=seeds, FUN=gen_alg_rare_wrap, mc.cores=20)

# save
  save(list="perm_results", file="perm_results.rdata")
  # load(file="perm_results.rdata")

