#include <Rcpp.h>
using namespace Rcpp;

// wrapper around R's RNG such that we get a uniform distribution over
// [0,n) as required by the STL algorithm
inline int randWrapper(const int n) { return floor(unif_rand()*n); }

// function to make collectors curve out of otutable otb with samps as cols
// minval is the minimum value to count a positive hit
// [[Rcpp::export]]
IntegerVector collect_samps_otb(const NumericMatrix& otb, const float minval) {
	int nc = otb.ncol();
	int nr = otb.nrow();
	IntegerVector output(nc);
	int ncm1 = nc - 1;
	// colidx is a vector of column indices for otb, so that otb doesn't need to
	// be manipulated (copied) in this function.
	IntegerVector colidx = seq(0, ncm1);
	// permute colidx
	std::random_shuffle(colidx.begin(), colidx.end(), randWrapper );
	// vector to store cumulative row sums
	NumericVector cumsum(nr);
	// loop through output (cols)
	int k=0;
	for(int j=0; j<nc; j++){
		// add each value otb[i,k] to cumsum[j]
		// where k is colidx[j]
		for(int i=0; i<nr; i++){
			k=colidx[j];
			cumsum[i] += otb(i,k);
			if(cumsum[i] > minval){
				output[j] += 1;
			}
		}
	}
	return(output);
}

