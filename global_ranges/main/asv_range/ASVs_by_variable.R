library(data.table)
### sub_by_samples(): subsets an OTU table to match a vector of sample names.
# Assumes that OTU table columns are samples, rows are otus.
# If a column of OTU names should kept, indicate the column name. 

sub_by_samples <- function(samples, otu_table, OTU_name_column = NA) {
  
  if(is.na(OTU_name_column)){
    # No OTU name column present
    
    # subset columns by selected sample names
    new_otu  <- otu_table[ , ..samples]
    
    # remove OTUs with rowSums that = 0
    new_otu[, Sums := rowSums(.SD)]
    
  } else {
    # OTU name column present, add it to the vector of sample names
    samples <- c(OTU_name_column, samples)
    
    # subset columns by selected sample names, preserving OTU names
    new_otu  <- otu_table[ , ..samples]
    
    # remove OTUs with rowSums that = 0, exlcuding OTU name column
    new_otu[, Sums := rowSums(.SD), .SDcols = !OTU_name_column]
  }
  
  new_otu <- new_otu[Sums != 0, .SD, .SDcols = !"Sums"]
  return(new_otu)
}

# In all other functions, OTUs and metadata will always be linked using column standardID
set_standardID <- function(metadata, id_column){
  ids <- metadata[[id_column]]
  metadata[, standardID := ids]
}



# Identify data files

otu_file  <- "../../biom/subset_global16S/global16S_subset.tsv"
meta_file <- "../../biom/subset_global16S/clean_map_global16S.tsv"

# Read in ----------

print("reading in data")

tm <- proc.time()

all_otus <- fread(otu_file)
colnames(all_otus)[1] <- "OTU_ID"

all_meta <- fread(meta_file)
colnames(all_meta)[1] <- "SampleID"
keep_samples <- all_meta$SampleID %in% colnames(all_otus)
all_meta <- all_meta[keep_samples]

proc.time() - tm

# set standard ID for meta
all_meta <- set_standardID(all_meta, "SampleID")

# Subset waimea samples
wai_meta <- all_meta[site_type == "Core Waimea"]

wai_otu <- sub_by_samples(samples = wai_meta$standardID,  otu_table = all_otus, OTU_name_column = "OTU_ID")

# list ASVs in each habitat
ASV_list <- list()

print("making lists")
tm <- proc.time()

for( hab in unique(wai_meta$habitat)){
 hab_otu <- sub_by_samples(samples = wai_meta[habitat %in% hab, SampleID],
				 otu_table = wai_otu,
				 OTU_name_column = "OTU_ID")

 ASV_list[[hab]] <- hab_otu$OTU_ID
 print(hab)
}

proc.time() - tm
print("getting ASVs by site")
tm <- proc.time()
for( site in unique(wai_meta$site_name)){
 site_otu <- sub_by_samples(samples = wai_meta[site_name %in% site, SampleID],
				 otu_table = wai_otu,
				 OTU_name_column = "OTU_ID")
 ASV_list[[site]] <- site_otu$OTU_ID
 print(site)
}
proc.time() - tm

print("getting ASVs by trophic")
tm <- proc.time()

wai_meta[!(trophic %in% c("PrimaryProducer","Environmental")), trophic := "Consumer"]
for( guild in unique(wai_meta$trophic)){
  guild_otu <- sub_by_samples(samples = wai_meta[trophic %in% guild, SampleID],
				 otu_table = wai_otu,
				 OTU_name_column = "OTU_ID")
 ASV_list[[guild]] <- guild_otu$OTU_ID
print(guild)
}

print("ASVs by trophic and habitat")
wai_meta[, habi_troph := paste(habitat, trophic)]
for(habitroph in unique(wai_meta$habi_troph)){
  habitroph_otu <- sub_by_samples(samples = wai_meta[habi_troph %in% habitroph, SampleID],
				 otu_table = wai_otu,
				 OTU_name_column = "OTU_ID")
 ASV_list[[habitroph]] <- habitroph_otu$OTU_ID
print(habitroph)
}

# write out
# capture.output(print(ASV_list), file= "ASVS_by_variable_group.txt")
saveRDS(ASV_list, file= "ASVs_by_variable.rds")

