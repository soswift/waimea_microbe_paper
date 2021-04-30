
library(data.table)

# Get Occupancy and Abundance For ASVs ------------------

# occ_abund() gets the occupancy (n samples) and abundance (sum of sequences) for each OTU/ASV in a table
# assumes names are stored in column "OTU_ID"
occ_abund <- function(otu_table){
  otu_names <- otu_table$OTU_ID
  abund <- rowSums(otu_table[, .SD, .SDcols = !"OTU_ID"])
  occ   <- rowSums(otu_table[ , lapply(.SD, as.logical), .SDcols = !"OTU_ID"])
  
  return(data.frame(otu = otu_names, abundance = abund, occupancy = occ))
}  

# Files ------------------------------------------------
# Run occ_abund() on reasonably sized chunks of a large file
otu_file  <- "../../biom/subset_global16S/global16S_subset.tsv"

print("reading in")

# how many total lines in the file?
otu_row_n <- as.numeric(sub(" (.+)","", system(paste("wc -l ",otu_file), intern = TRUE)))

print(paste0("total rows = ",otu_row_n))

# split total rows into 100 chunks (skip first two rows, which have header info)
chunks <- floor(seq(from = 2, to = otu_row_n, length.out = 100))

all_occ_abund <- data.frame(otu = character(), abundance = integer(), occupancy = integer())

for(i in 1:length(chunks)) {
  print(paste("Reading rows",chunks[i]+1, " - ",chunks[i+1]))
  # read chunk
  chunk_otus <- fread( cmd=  paste0("sed -n '", chunks[i]+1 ,",", chunks[i+1],"p' ", otu_file),
                       header = F)
  colnames(chunk_otus)[1] <- "OTU_ID"
  chunk_occ_abund <- occ_abund(chunk_otus)

  all_occ_abund <- rbind(all_occ_abund, chunk_occ_abund)
  print(paste("Finished", i, "of", length(chunks), "chunks"))
}

write.csv(all_occ_abund, "all_ASV_occupancy_abundance.csv", row.names = F)


