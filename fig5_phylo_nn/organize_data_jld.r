# read in dists
	empdists <- read.delim("emp_nearest_emps.csv.gz", header=TRUE, sep=",", stringsAsFactors=FALSE)
	waidists <- read.delim("waimea_nearest_emps.csv.gz", header=TRUE, sep=",", stringsAsFactors=FALSE)
	emptax <- read.delim("emp_seq_to_annoation.txt.gz", sep="\t", header=FALSE, stringsAsFactors=FALSE)
	waitax <- read.delim("waimea_seq_to_annoation.txt.gz", sep="\t", header=FALSE, stringsAsFactors=FALSE)
	colnames(waitax) <- c("Seq", "Tax", "NA")
	colnames(emptax) <- c("Seq", "Tax", "NA")

# drop negative dists (WTF???)
	waidists <- waidists[waidists$dist > 0, ]
	empdists <- empdists[empdists$dist > 0, ]

# merge waitax and waidists
	waiall <- merge(waitax, waidists, by="Seq")
	rm(waidists, waitax)

# merge emptax and empdists - much more difficult! 
	# not all empdists$Seq are in emptax, but those taxs ARE in waitax!
	getemptax <- function(empseq){
		emptax_matches <- emptax$Seq == empseq
		if(any(emptax_matches)){
			return(emptax$Tax[ which(emptax_matches)[1] ])
		}else{
			waiall_matches <- waiall$Neighbor == empseq
			if(any(waiall_matches)){
				return(waiall$Tax[ which(waiall_matches)[1] ])
			}else{
				return("")
			}
		}
	}
	empall <- data.frame(
		Seq=empdists$Seq,
		Tax="", # empty, will be filled in loop below
		"NA"=NA, # just to preserve matchup with waiall
		Neighbor=empdists$Neighbor,
		dist=empdists$dist,
		stringsAsFactors=FALSE
	)
	pb <- txtProgressBar(min=0, max=length(unique(empall$Seq)), style=3)
	for(s in unique(empall$Seq)){
		empall$Tax[empall$Seq == s] <- getemptax(s)
		setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
	}
	close(pb)
	colnames(empall)[3] <- "NA"

# merge data into one big df
	plotdata <- rbind(waiall, empall)
	
# deal with taxonomy
	# throw out unknown taxonomy
	plotdata <- plotdata[plotdata$Tax != "",]
	# get rid of garbage characters
	plotdata$Tax <- gsub(x=plotdata$Tax, pattern="\'", replacement="")
	plotdata$Tax <- gsub(x=plotdata$Tax, pattern="\\[", replacement="")
	plotdata$Tax <- gsub(x=plotdata$Tax, pattern="\\]", replacement="")
	plotdata$Tax <- gsub(x=plotdata$Tax, pattern=" ", replacement="")
	# split, add to df
	taxmat <- do.call("rbind",  strsplit(plotdata$Tax, split=","))
	# add column names
	colnames(taxmat) <- c("kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
	# remove gg "k__" stuff
	taxmat <- sub(pattern=".__", replacement="", x=taxmat)
	# recombine
	plotdata <- data.frame(plotdata, taxmat)

# save
	write.table(plotdata, file="data_jld_organized.txt", sep="\t", row.names=F, quote=F)
	system("gzip data_jld_organized.txt")




