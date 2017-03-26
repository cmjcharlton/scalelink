// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

struct LongLongMatrix {
	std::size_t nrow;
	std::size_t ncol;
	std::vector<std::array<uint64_t, 4> > data;

	LongLongMatrix() : nrow(0), ncol(0) {
	}

	LongLongMatrix(std::size_t nrows, std::size_t ncols) : nrow(nrows), ncol(ncols), data(nrows * ncols) {
	}

	LongLongMatrix(const LongLongMatrix &mat) : nrow(mat.nrow), ncol(mat.ncol) {
		data.assign(mat.data.begin(), mat.data.end());
	}

	uint64_t& operator()(std::size_t i, std::size_t j, std::size_t k) {
		return data[i * ncol + j][k];
	}

	uint64_t operator()(std::size_t i, std::size_t j, std::size_t k) const {
		return data[i * ncol + j][k];
	}
};

LongLongMatrix& operator += (LongLongMatrix & result, const LongLongMatrix & rhs) {
	for (size_t i = 0u; i < result.data.size(); ++i) {
		result.data[i][0] += rhs.data[i][0];
		result.data[i][1] += rhs.data[i][1];
		result.data[i][2] += rhs.data[i][2];
		result.data[i][3] += rhs.data[i][3];
	}
	//std::transform(result.data.begin(), result.data.end(), rhs.data.begin(), result.data.begin(), std::plus<uint64_t>());
	return result;
}

LongLongMatrix operator + (LongLongMatrix lhs, const LongLongMatrix & rhs) {
	return lhs += rhs;
}

struct Accumulate : public Worker
{
	// source data
	const RMatrix<double> foi;
	const RMatrix<double> ldf;

	// accumulated value
	LongLongMatrix accum;

	// temporary storage
	std::vector<int> foiexp;

	// constructors
	Accumulate(const NumericMatrix foi, const NumericMatrix ldf) : foi(foi), ldf(ldf) {
		std::size_t nident = foi.ncol();
		foiexp.resize(nident);
		accum = LongLongMatrix(nident, nident);
	}
	Accumulate(const Accumulate& acc, Split) : foi(acc.foi), ldf(acc.ldf) {
		std::size_t nident = foi.ncol();
		foiexp.resize(nident);
		accum = LongLongMatrix(nident, nident);
	}
	// accumulate just the element of the range I've been asked to
	void operator()(std::size_t begin, std::size_t end) {
		for (std::size_t i = begin; i < end; ++i) {
			for (std::size_t j = 0u; j < ldf.nrow(); ++j) {
				for (std::size_t x = 0u; x < foi.ncol(); ++x) {
					foiexp[x] = (foi(i, x) != ldf(j, x));
				}
				for (std::size_t x = 0u; x < foiexp.size(); ++x) {
					for (std::size_t y = 0u; y <= x; ++y) {
						accum(x, y, (foiexp[x]*2) + foiexp[y])++;
					}
				}
			}
		}
		// Rcpp::Rcout << "Processed FOI records " << begin << " to " << end << std::endl;
	}
     
	// join my value with that of another Accumulate
	void join(const Accumulate& rhs) {
		accum += rhs.accum;
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
	uint64_t nrecfoi = foinew.nrow();
	uint64_t nrecldf = ldfnew.nrow();
	std::size_t nident = foinew.ncol();
	uint64_t total = nrecfoi * nrecldf;

	NumericMatrix astar((nident*2) + 2, (nident*2) + 2);

	// declare the Accumulate instance 
	Accumulate acc(foinew, ldfnew);
   
	if (grainsize < nrecfoi) {
		// call parallel_reduce to start the work
		parallelReduce(0, nrecfoi, acc, grainsize);
	} else {
		// run the calculations in series
		std::size_t nreport = 10;
		std::size_t step = nrecfoi / 100 * nreport;
		std::size_t nextstep = step;
		std::size_t cent = nreport;

		// std::vector<std::vector<unsigned long long> > accum(nident*2, std::vector<unsigned long long>(nident*2));
		std::vector<int> foiexp(nident);

		Rcpp::Rcout << "Progress: |0";
		std::flush(Rcpp::Rcout);
		for (std::size_t i = 0u; i < nrecfoi; ++i) {
			for (std::size_t j = 0u; j < nrecldf; ++j) {
				for (std::size_t x = 0u; x < foiexp.size(); ++x) {
					foiexp[x] = (foinew(i, x) != ldfnew(j, x));
				}
				for (std::size_t x = 0u; x < foiexp.size(); ++x) {
					for (std::size_t y = 0u; y <= x; ++y) {
						acc.accum(x, y, (foiexp[x]*2) + foiexp[y])++;
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
		for (std::size_t i = 0u; i < acc.accum.nrow; i++) {
			for (std::size_t j = 0u; j < acc.accum.ncol; j++) {
				for (std::size_t k = 0u; k < 2; k++) {
					Rcpp::Rcout << acc.accum(i, j, k) << " ";
				}
			}
			Rcpp::Rcout << std::endl;
			for (std::size_t j = 0u; j < acc.accum.ncol; j++) {
				for (std::size_t k = 2u; k < 4; k++) {
					Rcpp::Rcout << acc.accum(i, j, k) << " ";
				}
			}
			Rcpp::Rcout << std::endl;
		}
	}

	// Test that sums add up
	bool sumerr = false;
	for (std::size_t i = 0u; i < nident; ++i) {
		for (std::size_t j = 0u; j < i; ++j) {
			uint64_t test = total;
			test -= acc.accum(i, j, 0);
			test -= acc.accum(i, j, 1);
			test -= acc.accum(i, j, 2);
			test -= acc.accum(i, j, 3);
			if (test != 0) {
				sumerr = true;
			}
		}
	}

	for (std::size_t i = 0u; i < nident; ++i) {
		astar(i*2, i*2) = acc.accum(i, i, 0) * (nident - 1.0) / (nident * nident) * 2.0;
		astar((i*2)+1, (i*2)+1) = acc.accum(i, i, 3) * (nident - 1.0) / (nident * nident) * 2.0;
		for (std::size_t j = 0u; j < i; j++) {
			astar(i*2, j*2) -= acc.accum(i, j, 0) / (nident * nident) * 2.0;
			astar(i*2, (j*2)+1) -= acc.accum(i, j, 1) / (nident * nident) * 2.0;
			astar((i*2)+1, j*2) -= acc.accum(i, j, 2) / (nident * nident) * 2.0;
			astar((i*2)+1, (j*2)+1) -= acc.accum(i, j, 3) / (nident * nident) * 2.0;
		}
	}

	// Fill in other half of symmetric matrix
	for (std::size_t i = 0u; i < nident * 2; ++i) {
		for (std::size_t j = i + 1; j < nident * 2; j++) {
			astar(i, j) = astar(j, i);
		}
	}

	// Add q and r

	// -q
	for (std::size_t i = 0u; i < nident; ++i) {
		astar(i * 2, nident * 2) -= 1.0 / nident;
		astar(i * 2, nident * 2) *= total;
	}

	// -r
	for (std::size_t i = 0u; i < nident; ++i) {
		astar((i * 2) + 1, nident * 2 + 1) -= 1.0 / nident;
		astar((i * 2) + 1, nident * 2 + 1) *= total;
	}

	// q
	for (std::size_t i = 0u; i < nident; ++i) {
		astar(nident * 2, i * 2) = 1.0 / nident;
		astar(nident * 2, i * 2) *= total;
	}

	// r
	for (std::size_t i = 0u; i < nident; ++i) {
		astar(nident * 2 + 1, (i * 2) + 1) = 1.0 / nident;
		astar(nident * 2 + 1, (i * 2) + 1) *= total;
	}

	if (sumerr == true) {
		Rcpp::Rcout << "Warning: unexpected sum (should be " << total << "), possible overflow" << std::endl;
	}

	return(astar);
}