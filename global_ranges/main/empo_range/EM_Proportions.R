# Document Summary:
# This script corellates the number of samples an OTU is found in with it's global geographic range
# Written as a script to run on HPC.

library(data.table)

# Define Functions ----------
# sum_sample_reads() gets the total reads for each sample in an otu table and adds it to the metadata
# Calculate proportion of reads that represent each EMPO3


# Define Inputs ----------

# Toy data version
 # otu_file  <- "../../toy_data/toy_data/toy_emp_data.csv"
 # meta_file <- "../../toy_data/toy_data/toy_emp_metadata.csv"

# Read data version
 otu_file  <- "../../biom/subset_global16S/global16S_subset.tsv"
 meta_file <- "../../biom/subset_global16S/clean_map_global16S.tsv"

# Read in data in chunks ----------
# reading all at once = memory disaster

print("reading in data")

tm <- proc.time()

print("reading in")
# pull out header names
sample_columns <- names(fread(otu_file, nrows = 3))[-1]

# how many total lines are in the file?
otu_row_n <- as.numeric(sub(" (.+)","", system(paste("wc -l ",otu_file), intern = TRUE)))
print(paste0("total rows = ",otu_row_n))

# split total rows into 100 chunks (skip first two rows, which have header info)
chunks <- floor(seq(from = 2, to = otu_row_n, length.out = 100))

# define vector  of zeroes to hold sums
sample_sums <-rep.int(0, times =  length(sample_columns))
names(sample_sums) <- sample_columns

# read in chunks sequentially
i <- 1
while(i < length(chunks)) {
  
  print(paste("Reading rows",chunks[i]+1, " - ",chunks[i+1]))
  
  chunk_otus <- fread( cmd=  paste0("sed -n '", chunks[i]+1 ,",", chunks[i+1],"p' ", otu_file) , header = F)
  chunk_otus[ , V1 := NULL]
  sample_sums <- sample_sums + colSums(chunk_otus)
  
  print(paste("Finished", i, "of", length(chunks), "chunks"))
  i <- i+1
}

proc.time() - tm

rm(chunk_otus)
gc()
# Read in meta 
all_meta <- fread(meta_file)
colnames(all_meta)[1] <- "SampleID"
keep_samples <- all_meta$SampleID %in% sample_columns
all_meta <- all_meta[keep_samples]

# set standard ID for meta
all_meta <- set_standardID(all_meta, "SampleID")

# add 'Dataset' column to identify EMP from Waimea
all_meta[ , Dataset := ifelse(site_type == "Core Waimea", "Waimea", "EMP")]

## Calculate proportions for EMPO 3 --------------------------------------

# Calculate the proportion of samples that represent each EMPO3
prop.samp <- all_meta[ , .(sample_count = .N), by = .(Dataset, empo_3) ]

# add 'read_count' column for each sample. 
all_meta[ , read_count := sample_sums[SampleID]]

# Calculate proportion of sequences that represent each EMPO3
prop.seq <- all_meta[ , .(read_sum = sum(read_count)), by = .(Dataset, empo_3)]

# Write out the sums by empo_3 and the raw read counts by sample (just in case)
all_props <-merge(prop.samp, prop.seq, by = c("Dataset", "empo_3"))
fwrite(all_props, "all_empo3_proportions.csv")
fwrite(all_meta[ , .(SampleID, read_count)], "read_counts_by_sample.csv")
