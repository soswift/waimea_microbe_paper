#include <Rcpp.h>
using namespace Rcpp;




inline IntegerVector getSelRows(IntegerVector sel, int nrow){
	IntegerVector output(sel.size());
	for(int i=0; i < sel.size(); i++){
		output[i] =  sel[i] % nrow;
	}
	return(output);
}

inline IntegerVector getSelCols(IntegerVector sel, double nrow){
	IntegerVector output(sel.size());
	for(int i=0; i < sel.size(); i++){
		output[i] =  floor(sel[i] / nrow);
	}
	return(output);
}

// test code for above (read into Rcpp without 'inline')
// matrix(getSelRows(0:24, 5), nrow=5)
// matrix(getSelCols(0:24, 5), nrow=5)




/* 
  Pvec: as.vector(P) of some rxc matrix of nonzero probabilities P
    (probabilities as in numeric values relative to each other, see ?sample)
  Nrows: nrow(P)
  Ncols: ncol(P)
  nSel: Number of probabilities to select per iteration. 
  noisy: display progress T/F
  iters2report: how many iterations pass before reporting progress if noisy.
*/
// [[Rcpp::export]]
IntegerMatrix vazfill(NumericVector Pvec, int nRows, int nCols, 
	int nSel=5, bool noisy=true, int iters2report=100){

	// clone Pvec and nSel so as to not modify them
	NumericVector Pvec2 = clone(Pvec);

	// check that nRows * nCols == Pvec2.size()
	if((nRows * nCols) != Pvec2.size()){
		stop("Error: Nrows * Ncols unequal to length(Pvec).");
	}

	int nFilledRows = 0;
	int nFilledCols = 0;
	IntegerVector filledCols(nCols);
	IntegerVector filledRows(nRows);
	int dimSum = nRows + nCols;

	IntegerVector positions = seq_len(Pvec2.size()) - 1;

	IntegerVector sel(nSel);
	IntegerVector selRows(nSel);
	IntegerVector selCols(nSel);

	IntegerMatrix outputMat(nRows, nCols);

	int n_unused= nRows * nCols;
	// Rcout << "n_unused = " << n_unused <<"\n";

	int noisecounter = 0;

	while((nFilledRows + nFilledCols) < dimSum){
		// sample cels to populate using Pvec2, find their row and col indices
		// if not enough values left to sample, sample fewer!
		if(n_unused < nSel){
			stop("Ran out of values to sample. nSel too high.");
		}else{
			sel = sample(positions, nSel, false, Pvec2);
		}
		selRows = getSelRows(sel, nRows);
		selCols = getSelCols(sel, nRows);

		// for each selected cel, check if row OR col is empty. 
		for(int i=0; i<nSel; i++){
			if(outputMat(selRows[i], selCols[i]) > 0){
				Rcout << "REPEAT FILL ERROR" <<"\n";
				Rcout << "P = " << Pvec2[sel[i]]<<"\n";
			}
			if(filledCols[selCols[i]] == 0 || filledRows[selRows[i]] == 0){
				// put a 1 in the matrix
				outputMat(selRows[i],selCols[i]) = 1;
				// put a 1 in filledCols
				filledCols[selCols[i]] = 1;
				// put a 1 in filledRows
				filledRows[selRows[i]] = 1;
			}
			Pvec2[sel[i]] = 0;
			n_unused = n_unused - 1;
		}
		// update nFilledRows and nFilledCols
		nFilledRows = sum(filledRows);
		nFilledCols = sum(filledCols);
		// report on progress maybe
		if(noisy){
			noisecounter ++;
			if(noisecounter == iters2report){
				// report on progress
				Rcout << nFilledRows+nFilledCols <<" / " << dimSum <<"\n";
				noisecounter = 0;
			}
		}
	}

	return(outputMat);
}

/*
	R code for testing

	web <- matrix(sample(c(0,0,0,1,1,2,3), 5000, replace=TRUE), nrow=50, ncol=100)
	rs.p <- rowSums(web)/sum(web)
	cs.p <- colSums(web)/sum(web)
	P <- P1 <- tcrossprod(rs.p, cs.p)
	# P <- P * 1/min(P)

	range(P)

	sourceCpp("vazfill.cpp")
	test <- vazfill(as.vector(P), nrow(P), ncol(P))

	
*/