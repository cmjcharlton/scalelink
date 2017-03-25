// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

struct Accumulate : public Worker
{
	// source data
	const RMatrix<double> foi;
	const RMatrix<double> ldf;

	// accumulated value
	std::vector<std::vector<unsigned long long> > accum;

	// temporary storage
	std::vector<int> foiexp;

	// constructors
	Accumulate(const NumericMatrix foi, const NumericMatrix ldf) : foi(foi), ldf(ldf) {
		unsigned nident = foi.ncol();
		foiexp.resize(nident);
		accum.resize(nident * 2);
		for (unsigned i = 0u; i < nident * 2; ++i) {
			accum[i].resize(nident * 2);
		}
	}
	Accumulate(const Accumulate& acc, Split) : foi(acc.foi), ldf(acc.ldf) {
		unsigned nident = foi.ncol();
		foiexp.resize(nident);
		accum.resize(nident * 2);
		for (unsigned i = 0u; i < nident * 2; ++i) {
			accum[i].resize(nident * 2);
		}
	}
	// accumulate just the element of the range I've been asked to
	void operator()(std::size_t begin, std::size_t end) {
		for (std::size_t i = begin; i < end; ++i) {
			for (std::size_t j = 0u; j < ldf.nrow(); ++j) {

				for (std::size_t x = 0u; x < foi.ncol(); ++x) {
					if (foi(i, x) - ldf(j, x) != 0) {
						foiexp[x] = 1;
					} else {
						foiexp[x] = 0;
					}
				}

				for (std::size_t x = 0u; x < foiexp.size(); ++x) {
					for (std::size_t y = 0u; y <= x; ++y) {
						if (foiexp[x] & foiexp[y]) {
							accum[x*2][y*2]++;
						}
						if (!foiexp[x] & foiexp[y]) {
							accum[(x*2)+1][y*2]++;
						}
						if (foiexp[x] & !foiexp[y]) {
							accum[x*2][(y*2)+1]++;
						}
						if (!foiexp[x] & !foiexp[y]) {
							accum[(x*2)+1][(y*2)+1]++;
						}
					}
				}

/*
				for (std::size_t x = 0u; x < foiexp.size(); ++x) {
					for (std::size_t y = 0u; y <= x; ++y) {
						accum[x*2][y*2] += foiexp[x] & foiexp[y];
						accum[(x*2)+1][y*2] += !foiexp[x] & foiexp[y];
						accum[x*2][(y*2)+1] += foiexp[x] & !foiexp[y];
						accum[(x*2)+1][(y*2)+1] += !foiexp[x] & !foiexp[y];
					}
				}
*/



/*
				for (std::size_t x = 0u; x < foi.ncol(); ++x) {
					for (std::size_t y = 0u; y <= x; ++y) {
						accum[x*2][y*2] += (foi(i, x) == ldf(j, x)) & (foi(i, y) == ldf(j, y));
						accum[(x*2)+1][y*2] += (foi(i, x) != ldf(j, x)) & (foi(i, y) == ldf(j, y));
						accum[x*2][(y*2)+1] += (foi(i, x) == ldf(j, x)) & (foi(i, y) != ldf(j, y));
						accum[(x*2)+1][(y*2)+1] += (foi(i, x) != ldf(j, x)) & (foi(i, y) != ldf(j, y));
					}
				}
*/
			}
		}
		// Rcpp::Rcout << "Processed FOI records " << begin << " to " << end << std::endl;
	}
     
	// join my value with that of another Accumulate
	void join(const Accumulate& rhs) {
		for (std::size_t i = 0u; i < accum.size(); ++i) {
			for (std::size_t j = 0u; j < accum[i].size(); ++j) {
				accum[i][j] += rhs.accum[i][j];
			}
		}
	}
};

//' @title
//' buildAstar
//' @description
//' Builds the A* matrix
//' 
//' @param foinew numeric \code{\link[base]{matrix}} representing the file of interest
//' 
//' @param ldfnew numeric \code{\link[base]{matrix}} representing the linking data file
//'
//' @param grainsize integer determining minimum grain size for parallisation
//'
//' @param debug Boolean indicating whether to output additional debugging information
//' 
//' @details
//' \code{buildAstar} takes a matrix representing the file of interest and 
//' a matrix representing the linking data file and creates a matrix that 
//' can then be used to generating linking scores. Reporting frequency as this
//' occurs can be specified via the nreport option. This is implemented in C++
//' to provide a speed increase over implementing it directly in the R equivalent.
//' @export
// [[Rcpp::export]]
NumericMatrix buildAstar(NumericMatrix foinew, NumericMatrix ldfnew, int grainsize, bool debug) {
	unsigned long long nrecfoi = foinew.nrow();
	unsigned long long nrecldf = ldfnew.nrow();
	unsigned nident = foinew.ncol();
	unsigned long long total = nrecfoi * nrecldf;

	NumericMatrix astar((nident*2) + 2, (nident*2) + 2);

	// declare the Accumulate instance 
	Accumulate acc(foinew, ldfnew);
   
	if (grainsize < nrecfoi) {
		// call parallel_reduce to start the work
		parallelReduce(0, nrecfoi, acc, grainsize);
	} else {
		// run the calculations in series
		unsigned nreport = 10;
		unsigned step = nrecfoi / 100 * nreport;
		unsigned nextstep = step;
		unsigned cent = nreport;

		// std::vector<std::vector<unsigned long long> > accum(nident*2, std::vector<unsigned long long>(nident*2));
		std::vector<int> foiexp(nident * 2);

		Rcpp::Rcout << "Progress: |0";
		std::flush(Rcpp::Rcout);
		for (unsigned i = 0u; i < nrecfoi; ++i) {
			for (unsigned j = 0u; j < nrecldf; ++j) {
				for (unsigned k = 0u; k < nident; ++k) {
					if (foinew(i, k) - ldfnew(j, k) != 0) {
						foiexp[k*2] = 1;
						foiexp[(k*2) + 1] = 0;
					} else {
						foiexp[(k*2) + 1] = 1;
						foiexp[k*2] = 0;
					}
				}

				for (unsigned x = 0u; x < foiexp.size(); ++x) {
					if (foiexp[x] == 1) {
						for (unsigned y = 0u; y <= x; ++y) {
							acc.accum[x][y] += foiexp[y];
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
		Rcpp::Rcout << "..100|" << std::endl;
	}

	if (debug == true) {
		Rcpp::Rcout << "Accumulations: " << std::endl;
		for (unsigned i = 0u; i < acc.accum.size(); i++) {
			for (unsigned j = 0u; j < acc.accum.size(); j++) {
				Rcpp::Rcout << acc.accum[i][j] << " ";
			}
			Rcpp::Rcout << std::endl;
		}
	}

	// Test that sums add up
	bool sumerr = false;
	for (unsigned i = 0u; i < nident * 2; i += 2) {
		for (unsigned j = 0u; j < i; j += 2) {
			unsigned long long test = total;
			test -= acc.accum[i][j];
			test -= acc.accum[i+1][j];
			test -= acc.accum[i][j+1];
			test -= acc.accum[i+1][j+1];
			if (test != 0) {
				sumerr = true;
			}
		}
	}

	for (unsigned i = 0u; i < nident * 2; ++i) {
		astar(i, i) = acc.accum[i][i] * (nident - 1.0) / (nident * nident) * 2.0;
		for (unsigned j = 0u; j < i; j++) {
			astar(i, j) -= acc.accum[i][j] / (nident * nident) * 2.0;
		}
	}

	// Fill in other half of symmetric matrix
	for (unsigned i = 0u; i < nident * 2; ++i) {
		for (unsigned j = i + 1; j < nident * 2; j++) {
			astar(i, j) = astar(j, i);
		}
	}

	// Add q and r

	// -q
	for (unsigned i = 0u; i < nident; ++i) {
		astar(i * 2, nident * 2) -= 1.0 / nident;
		astar(i * 2, nident * 2) *= total;
	}

	// -r
	for (unsigned i = 0u; i < nident; ++i) {
		astar((i * 2) + 1, nident * 2 + 1) -= 1.0 / nident;
		astar((i * 2) + 1, nident * 2 + 1) *= total;
	}

	// q
	for (unsigned i = 0u; i < nident; ++i) {
		astar(nident * 2, i * 2) = 1.0 / nident;
		astar(nident * 2, i * 2) *= total;
	}

	// r
	for (unsigned i = 0u; i < nident; ++i) {
		astar(nident * 2 + 1, (i * 2) + 1) = 1.0 / nident;
		astar(nident * 2 + 1, (i * 2) + 1) *= total;
	}

	if (sumerr == true) {
		Rcpp::Rcout << "Warning: unexpected sum (should be " << total << "), possible overflow" << std::endl;
	}

	return(astar);
}