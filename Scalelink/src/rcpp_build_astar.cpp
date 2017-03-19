
#include <Rcpp.h>
using namespace Rcpp;

//' @title
//' buildAstar
//' @description
//' Builds the A* matrix
//' 
//' @param foinew numeric \code{\link[base]{matrix}} representing the file of interest
//' 
//' @param ldfnew numeric \code{\link[base]{matrix}} representing the linking data file
//'
//' @param nreport integer determining the frequency of progress reporing
//'
//' @debug Boolean indicating whether to output additional debugging information
//' 
//' @details
//' \code{buildAstar} takes a matrix representing the file of interest and 
//' a matrix representing the linking data file and creates a matrix that 
//' can then be used to generating linking scores. Reporting frequency as this
//' occurs can be specified via the nreport option. This is implemented in C++
//' to provide a speed increase over implementing it directly in the R equivalent.
//' @export
// [[Rcpp::export]]
NumericMatrix buildAstar(NumericMatrix foinew, NumericMatrix ldfnew, int nreport, bool debug) {
	unsigned long long nrecfoi = foinew.nrow();
	unsigned long long nrecldf = ldfnew.nrow();
	unsigned nident = foinew.ncol();
	unsigned long long total = nrecfoi * nrecldf;

	unsigned step = nrecfoi / 100 * nreport;
	unsigned nextstep = step;
	unsigned cent = nreport;

	NumericMatrix astar((nident*2) + 2, (nident*2) + 2);

	std::vector<std::vector<unsigned long long> > accum(nident*2, std::vector<unsigned long long>(nident*2));
	std::vector<int> foiexp(nident * 2);

	Rcpp::Rcout << "Progress: |0";
	std::flush(Rcpp::Rcout);
	for (unsigned i = 0u; i < nrecfoi; i++) {
		for (unsigned j = 0u; j < nrecldf; j++) {
			for (unsigned k = 0u; k < nident; k++) {
				if (foinew(i, k) - ldfnew(j, k) != 0) {
					foiexp[k*2] = 1;
					foiexp[(k*2) + 1] = 0;
				} else {
					foiexp[(k*2) + 1] = 1;
					foiexp[k*2] = 0;
				}
			}

			for (unsigned x = 0u; x < foiexp.size(); x++) {
				if (foiexp[x] == 1) {
					for (unsigned y = 0u; y <= x; y++) {
						accum[x][y] += foiexp[y];
					}
				}
			}
		}
		if (i > nextstep && cent < 100) {
			Rcpp::Rcout << ".." << cent;
			std::flush(Rcpp::Rcout);
			cent += nreport;
			nextstep += step;
		}
		//if (i % nreport == 0) {
		//	Rcpp::Rcout << "FOI record number: " << i << std::endl;
		//}
	}

	if (debug == true) {
		Rcpp::Rcout << "..100|" << std::endl;
		Rcpp::Rcout << "Accumulations: " << std::endl;
		for (unsigned i = 0u; i < accum.size(); i++) {
			for (unsigned j = 0u; j < accum.size(); j++) {
				Rcpp::Rcout << accum[i][j] << " ";
			}
			Rcpp::Rcout << std::endl;
		}
	}

	// Test that sums add up
	bool sumerr = false;
	for (unsigned i = 0u; i < nident * 2; i += 2) {
		for (unsigned j = 0u; j < i; j += 2) {
			unsigned long long test = total;
			test -= accum[i][j];
			test -= accum[i+1][j];
			test -= accum[i][j+1];
			test -= accum[i+1][j+1];
			if (test != 0) {
				sumerr = true;
			}
		}
	}

	for (unsigned i = 0u; i < foiexp.size(); i++) {
		astar(i, i) = accum[i][i] * (nident - 1.0) / (nident * nident) * 2.0;
		for (unsigned j = 0u; j < i; j++) {
			astar(i, j) -= accum[i][j] / (nident * nident) * 2.0;
		}
	}

	// Fill in other half of symmetric matrix
	for (unsigned i = 0u; i < foiexp.size(); i++) {
		for (unsigned j = i + 1; j < foiexp.size(); j++) {
			astar(i, j) = astar(j, i);
		}
	}

	// Add q and r

	// -q
	for (unsigned i = 0u; i < nident; i++) {
		astar(i * 2, foiexp.size()) -= 1.0 / nident;
		astar(i * 2, foiexp.size()) *= total;
	}

	// -r
	for (unsigned i = 0u; i < nident; i++) {
		astar((i * 2) + 1, foiexp.size() + 1) -= 1.0 / nident;
		astar((i * 2) + 1, foiexp.size() + 1) *= total;
	}

	// q
	for (unsigned i = 0u; i < nident; i++) {
		astar(foiexp.size(), i * 2) = 1.0 / nident;
		astar(foiexp.size(), i * 2) *= total;
	}

	// r
	for (unsigned i = 0u; i < nident; i++) {
		astar(foiexp.size() + 1, (i * 2) + 1) = 1.0 / nident;
		astar(foiexp.size() + 1, (i * 2) + 1) *= total;
	}
	if (debug == false) {
		Rcpp::Rcout << "..100|" << std::endl;
	}

	if (sumerr == true) {
		Rcpp::Rcout << "Warning: unexpected sum (should be " << total << "), possible overflow" << std::endl;
	}

	return(astar);
}