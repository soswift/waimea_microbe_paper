# make fake data
	# rows are empo3 categories ("samples")
	# cols asvs/otus/species/whatever
	# otutable <- matrix(
	# 	sample(c(0,1), prob=c(90, 1), size=19000, replace=T),
	# 	nrow=19
	# )
	# in real life this would have to be aggregated from
	# a full contingency table. not too difficult.

# define genetic algorithm function.
	# goal is to choose identities of n empo3 cats
	# that maximize diversity
	# takes otutable and n as arguments
	# otutable - see above
	# n - number of empo3 categories to get
	# popsize - how many random states to startwith
	# best2keep - how many of the best states to keep between generations
	# regs2keep - how many of the *NOT* best states to keep between generations
	# gens2stop - how many generations must pass with no improvement to trigger stopping
	# mutrate - how many "mutations" do do per iter.
	genOptRich <- function(otutable, n, popsize=1000,
		best2keep=50, regs2keep=25, gens2stop=0,
		maxgens=10000, mutrate=2, noisy=F){
		ncats <- nrow(otutable)

		# function to calculate richness of a state
		# (state = a member of pop)
		# x (state) is a boolean vector with length
		# equal to nrow(otutable)
		getRich <- function(x){
			rs<- colSums(otutable[x, , drop=F])
			return(sum(rs > 0))
		}

		# check for edge cases - n==ncats, n==1, n>ncats, n<0:
		if(n==ncats){
			# don't need to run genetic algorithm.
			return(list(
				bestState=rep(T, ncats),
				richnessByGen=getRich(rep(T, ncats))
			))

		}else if(n==1){
			richs <- rowSums(otutable >0)
			bs <- rep(F, ncats)
			bs[which.max(richs)] <- T
			return(list(
				bestState=bs,
				richnessByGen=getRich(bs)
			))

		}else if(n>ncats){
			stop("n must be between 1 and ncol(otutable)")
		}else if(n<0){
			stop("n must be between 1 and ncol(otutable)")
		}

		# create starting population, randomly
		getStartingVec <- function(n, nc=ncats){
			out <- rep(F, nc)
			out[sample(1:nc, size=n, replace=F)] <- T
			return(out)
		}
		pop <- lapply(X=rep(n, popsize), FUN=getStartingVec)
		# defind mutation function
		# just finds a false (pT) and a true (pF) and
		# NOTs their values. This way, mutation happens,
		# but then number of trues and falses stay the same.
		mutate <- function(x, rate){
			for(i in 1:rate){
				pT <- sample(which(x), size=1)
				pF <- sample(which(!x), size=1)
				x[pT] <- !x[pT]
				x[pF] <- !x[pF]
			}
			return(x)
		}

		# function to calculate which positions of a
		# vector have the n highest values
		whichMaxN <- function(x, n){ order(x, decreasing=T)[1:n] }

		# function to do a generation
		generation <- function(pop, richs){
			# choose some to keep cuz they're good
			tokeep <- whichMaxN(richs, best2keep)
			# choose some to keep cuz they're lucky
			tokeep <- c(tokeep, sample((1:popsize)[! (1:popsize) %in% tokeep], size=regs2keep))

			# make new generation - start with tokeep
			newgen <- pop[tokeep]
			# do mutation to fill in new generation
			while(length(newgen) < popsize){
				newgen <- c(newgen, lapply(X=newgen, FUN=mutate, rate=mutrate))
			}
			# trim extras
			newgen <- newgen[1:popsize]

			return(newgen)
		}

		# variables for "how many gens with no improvement"
		# and "what was the last maximum"
		gensNoImprovment <- 0
		lastMax <- 0
		currentGen <- 0
		richPerGen <- rep(0, maxgens)

		# do loop because recursion... dangerous writing while(TRUE)!
		while(TRUE){
			# increment currentGen
			currentGen <- currentGen + 1

			# get richnesses
			richs <- sapply(X=pop, FUN=getRich)

			# store best richness
			richPerGen[currentGen] <- max(richs)

			if(gensNoImprovment <= gens2stop && currentGen <= maxgens){
				# make next gen
				pop <- generation(pop, richs)
			}else{
				break()
			}

			# figure out how many gens its been with no improvenment
			if(currentGen >= 2){
				if(richPerGen[currentGen] == richPerGen[currentGen - 1]){
					gensNoImprovment <- gensNoImprovment + 1
				}else{
					gensNoImprovment <- 0
				}
			}else{
				gensNoImprovment <- 0
			}

			if(noisy){
				message(paste0("(n=", n, ") gen=", currentGen,
					" rich=", richPerGen[currentGen]))
			}

		}

		# all done - make output and return it
		output <- list(
			bestState=pop[[which.max(richs)]],
			richnessByGen=richPerGen[1:currentGen]
		)
		return(output)

	}


