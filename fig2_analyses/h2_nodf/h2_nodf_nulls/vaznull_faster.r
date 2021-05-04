#' vaznull_faster
#'
#' A faster way to do the Vazquez null model for bipartite networks. Compared
#' to bipartite::vaznull, this is roughly 60 times faster. 
#'
#' @author John L. Darcy
#'
#' @param web matrix of integer values representing interactions.
#' @param samplesize integer, number of cells to fill simultaneously. Recommended
#'   to choose a value that's roughly 1/1000th the size of web. Higher values will
#'   give better speedup, but will reduce accuracy vs 1, which is the actual canonical
#'   Vazquez null model. 
#' @param noisy bool. should progress reports be printed?
#' @param noisefreq integer. how often (number of iterations) to report progress?
#' 
#' @return a matrix just like web but all null-model-y. 
#'
#' @examples
#'   # none yet written
#'   
#' 
#' @export
vaznull_faster <- function(web, samplesize, noisy=TRUE, noisefreq=100) {
	rs.p <- rowSums(web)/sum(web)
	cs.p <- colSums(web)/sum(web)
	P1 <- tcrossprod(rs.p, cs.p)
	if(min(P1)==0){P1 <- P1+min(P1[P1>0])}

	finalmat <- vazfill(Pvec=as.vector(P1), nRows=nrow(P1), nCols=ncol(P1),
		nSel=samplesize, noisy=noisy, iters2report=noisefreq)

	conn.remain <- sum(web > 0) - sum(finalmat > 0)
	if (conn.remain > 0) {
		if (length(which(finalmat == 0)) == 1) {
			add <- which(finalmat == 0)
		}else{
			add <- sample(which(finalmat == 0), conn.remain, prob = P1[finalmat == 0])
		}
		finalmat[add] <- 1
	}
	int.remain <- sum(web) - sum(finalmat)
	if (int.remain > 0) {
		add <- sample(which(finalmat > 0), int.remain, prob = P1[which(finalmat > 
			0)], replace = TRUE)
		finalmat[as.numeric(names(table(add)))] <- finalmat[as.numeric(names(table(add)))] + 
			table(add)
	}
	return(finalmat)
}

