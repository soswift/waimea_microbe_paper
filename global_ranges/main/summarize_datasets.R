# R script summarizing the HPC dataset and the EMP dataset(s)
library(data.table)
library(parallel)
set.seed(2020)

options(mc.cores = 8)

# Define Functions ---------------------

## read_hpc_abund reads in abundance file from 16S metaflowmics pipeline.
# Cleans up columns, transposes so OTUs are rows and samples are columns.

read_hpc_abund <- function(abundance_file){
  
  # read in
  abundance_table <- fread(abundance_file)
  
  # drop extra columns
  abundance_table[, c('label','numOtus') :=NULL]
  
  # transpose, keep OTU names
  abundance_table <- transpose(abundance_table, keep.names = "OTU_ID", make.names = "Group")
  
  setkey(abundance_table,"OTU_ID")
  
  return(abundance_table)
}

# Define inputs ---------------------

# define input file names

# all EMP
emp_sums_csv <- "sample-frequency-detail.csv"

# subset EMP
emp_abund_tsv <- "/home/seanos/cmaiki-lts/seanos/waimea/biom/subset_global16S/global16S_subset.tsv"
emp_meta_tsv  <- "/home/seanos/cmaiki-lts/seanos/waimea/biom/subset_global16S/clean_map_global16S.tsv"

# all UH HPC
hpc_abund_tsv <- "/home/seanos/cmaiki-lts/seanos/waimea/pipeline/waimea_16S_v3/abundance_table_100.shared"
hpc_meta_csv  <- "/home/seanos/cmaiki-lts/seanos/waimea/pipeline/waimea_16S_v3/combined_prep_sample_metadata_16s.csv"

# reduced UH HPC
reduced_hpc_csv <- "mm_16s_hiseqs_abundance_table.csv"

# HPC ------------------------
# Read in and subset data from the Metaflowmics pipeline running on the UH-HPC
# These data contain only samples from Waimea


# 1. This is the full waimea data set (i.e. all ASVs retained)

# read in abundance file
hpc_abund <- read_hpc_abund(hpc_abund_tsv)

# read in metadata file
hpc_meta <- fread(hpc_meta_csv)

# subset metadata to samples that sequenced succesfully
hpc_meta <- hpc_meta[sequencing_id %in% colnames(hpc_abund)]

# subset to project samples
hpc_meta <- hpc_meta[project_name %in% c("Keck","JoshMS")]
samples_to_keep <- as.character(hpc_meta$sequencing_id)


hpc_abund <- hpc_abund[ ,..samples_to_keep]

sink(file = "hpc_all_summary.txt")
 print("hpc_all")
 print(paste0("mean seq depth = ", mean(colSums(hpc_abund))))
 print(paste0("sd of seq depth = ", sd(colSums(hpc_abund))))
 print(paste0("total samples = ", ncol(hpc_abund)))
 print(paste0("total ASVs = ", nrow(hpc_abund)))
sink()

rm(hpc_abund)

# 2. This is the reduced dataset.
# cleaned using the clean_16S function, removed taxa based on rarity and abundance
# these data were used for ordinations and heatmaps
# analyses based on distance matrices did not benefit from rare taxa

# read in
reduced_hpc <- fread(reduced_hpc_csv, header = T)
hpc_meta <- fread(hpc_meta_csv)

# subset to project samples
hpc_meta <- hpc_meta[project_name %in% c("Keck","JoshMS")]
samples_to_keep <- as.character(hpc_meta$sequencing_id)

reduced_hpc <- reduced_hpc[ Group %in% samples_to_keep]

# drop rownames
reduced_hpc[ , Group:= NULL]

reduced_hpc <- reduced_hpc[rowSums(reduced_hpc)>0 ,]

sink(file = "hpc_reduced_summary.txt")
print("hpc_reduced")
print(paste0("mean seq depth = ", mean(rowSums(reduced_hpc))))
print(paste0("sd of seq depth = ", sd(rowSums(reduced_hpc))))
print(paste0("total samples = ", nrow(reduced_hpc)))
print(paste0("total ASVs = ", ncol(reduced_hpc)))
sink()


# EMP ------------------------------------------

# 1. All Waimea and EMP
# File is too big to read in
# Happily, Qiita provides read sums by sample
all_emp_sums <-fread(emp_sums_csv)

sink("emp_all_summary.txt")
print(paste0("mean seq depth = ", mean(all_emp_sums$V2)))
print(paste0("sd of seq depth = ", sd(all_emp_sums$V2)))
print(paste0("total samples = ", length(all_emp_sums$V2)))
print("total ASVs = 1911880") # see Qiita portal for this information
sink()

# 2. EMP global16S subset
# these samples contain all EMP and waimea samples, limited to ASVs that occured in waimea
# read in abundance file
emp_abund <- fread(emp_abund_tsv)
colnames(emp_abund)[1] <- "OTU_ID"
emp_abund[ , OTU_ID := NULL]
# read in metadata file
emp_meta <- fread(emp_meta_tsv)
colnames(emp_meta)[1] <- "SampleID"

emp_meta <- emp_meta[SampleID %in% colnames(emp_abund),]

# efficiently sum reads by sample
samp_sums <- emp_abund[ ,lapply(.SD, sum)]
samp_sums <- transpose(samp_sums)

# subset to non-Waimea samples
waimea_samples <- emp_meta[ geo_loc_name == "Waimea Valley", SampleID]
emp_abund[ , c(waimea_samples) := NULL]
asv_presence <- apply(emp_abund, 1, function(x) any(x > 0))



sink("emp_subset_summary.txt")
print(paste0("mean seq depth = ", mean(samp_sums$V1)))
print(paste0("sd of seq depth = ", sd(samp_sums$V1)))
print(paste0("total samples = ", ncol(emp_abund)))
print(paste0("total ASVs = ", nrow(emp_abund)))
print(paste0("Waimea ASVs in EMP = ", length(asv_presence[asv_presence == TRUE])))
