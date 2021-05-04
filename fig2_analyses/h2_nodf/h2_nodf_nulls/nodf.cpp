#include <Rcpp.h>
using namespace Rcpp;

inline int count_if(LogicalVector x) {
    int counter = 0;
    for(int i = 0; i < x.size(); i++) {
        if(x[i] == TRUE) {
            counter++;
        }
    }
    return counter;
}


// [[Rcpp::export]]
double NODFcpp(const LogicalMatrix& web){
	// tally row and col numbers
	int nrows = web.nrow();
	int ncols = web.ncol();
	// do rorws
	IntegerVector MTrows = rowSums(web); // marginal totals for rows
	double NpRows=0; // npaired for rows
	for(int i=0; i < nrows; i++){
		for(int j=0; j < nrows; j++){
			if(MTrows[i] > MTrows[j]){
				NpRows += count_if(web(i,_) * web(j,_)) / (double)MTrows[j];
			}
		}
	}

	IntegerVector MTcols = colSums(web);   // marginal totals for columns
	double NpCols=0;   // npaired for rows
	for(int i=0; i < ncols; i++){
		for(int j=0; j < ncols; j++){
			if(MTcols[i] > MTcols[j]){
				NpCols += count_if(web(_,i) * web(_,j)) / (double)MTcols[j];
			}
		}
	}

	double numerator = (NpRows + NpCols) * 200;
	double denominator = ((nrows * (nrows - 1)) + (ncols * (ncols - 1)));

	return(numerator / denominator);
}

