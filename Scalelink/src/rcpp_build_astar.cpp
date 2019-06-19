// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

struct LongLongMatrix {
	std::size_t nrow;
	std::size_t ncol;
	std::vector<uint64_t> data;

	LongLongMatrix() : nrow(0), ncol(0) {
	}

	LongLongMatrix(std::size_t nrows, std::size_t ncols) : nrow(nrows), ncol(ncols), data(nrows * ncols) {
	}

	LongLongMatrix(const LongLongMatrix &mat) : nrow(mat.nrow), ncol(mat.ncol) {
		data.assign(mat.data.begin(), mat.data.end());
	}

	uint64_t& operator()(std::size_t i, std::size_t j) {
		return data[i * ncol + j];
	}

	uint64_t operator()(std::size_t i, std::size_t j) const {
		return data[i * ncol + j];
	}
};

LongLongMatrix& operator += (LongLongMatrix & result, const LongLongMatrix & rhs) {
	//for (size_t i = 0u; i < result.data.size(); ++i) {
	//	result.data[i] += rhs.data[i];
	//}
	std::transform(result.data.begin(), result.data.end(), rhs.data.begin(), result.data.begin(), std::plus<uint64_t>());
	return result;
}

LongLongMatrix operator + (LongLongMatrix lhs, const LongLongMatrix & rhs) {
	return lhs += rhs;
}

tbb::mutex m;
bool dbg = false;

struct Accumulate : public Worker
{
	// type of each column
	std::vector<char> coltype;

	// source data
	std::vector<RVector<int> > foicols_int;
	std::vector<RVector<int> > ldfcols_int;
	std::vector<RVector<double> > foicols_num;
	std::vector<RVector<double> > ldfcols_num;
	std::vector<std::vector<bool> > foicols_lgl; // Note, we cannot use RVector<bool> here as this does not wrap LogicalVector
	std::vector<std::vector<bool> > ldfcols_lgl;
	std::vector<std::vector<std::string> > foicols_str;
	std::vector<std::vector<std::string> > ldfcols_str;

	// comparison functions
	std::vector<std::function<int(const int, const int)> > intfuncs;
	std::vector<std::function<int(const double, const double)> > realfuncs;
	std::vector<std::function<int(const std::string&, const std::string&)> > strfuncs;

	// accumulated value
	LongLongMatrix accum;

	// dimensions
	const size_t nvar;
	const size_t nfoi;
	const size_t nldf;

	// start of each identifier
	std::vector<uint64_t> idbase;

	// total categories
	const int totcats;

	// temporary storage
	std::vector<int> foiexp;

	// debugging
	std::vector<std::pair<int, int> > blockrange;

	// constructors
	Accumulate(const DataFrame foi, const DataFrame ldf, const std::vector<char> &coltype, std::vector<std::function<int(const int, const int)> > intfuncs, std::vector<std::function<int(const double, const double)> > realfuncs, std::vector<std::function<int(const std::string&, const std::string&)> > strfuncs, std::vector<uint64_t> &idbase, const int totcats) :
		coltype(coltype), intfuncs(intfuncs), realfuncs(realfuncs), strfuncs(strfuncs), nvar(coltype.size()), nfoi(foi.nrows()), nldf(ldf.nrows()), idbase(idbase), totcats(totcats) {
		for (std::size_t x = 0u; x < nvar; ++x) {
			switch (coltype[x]) {
			case 'I':
				{
					IntegerVector foivec = foi[x];
					foicols_int.emplace_back(RVector<int>(foivec));

					IntegerVector ldfvec = ldf[x];
					ldfcols_int.emplace_back(RVector<int>(ldfvec));
				}
				break;
			case 'R':
				{
					NumericVector foivec = foi[x];
					foicols_num.emplace_back(RVector<double>(foivec));

					NumericVector ldfvec = ldf[x];
					ldfcols_num.emplace_back(RVector<double>(ldfvec));
				}
				break;
			case 'L':
				{
					LogicalVector foivec = foi[x];
					std::vector<bool> foitmp(nfoi);
					std::copy(foivec.cbegin(), foivec.cend(), foitmp.begin()); // Copy LogicalVector elements to bool
					foicols_lgl.emplace_back(foitmp);

					LogicalVector ldfvec = ldf[x];
					std::vector<bool> ldftmp(nldf);
					std::copy(ldfvec.cbegin(), ldfvec.cend(), ldftmp.begin()); // Copy LogicalVector elements to bool
					ldfcols_lgl.emplace_back(ldftmp);
				}
				break;
			case 'S':
				{
					StringVector foivec = foi[x];
					std::vector<std::string> foitmp(nfoi);
					std::copy(foivec.cbegin(), foivec.cend(), foitmp.begin()); // Copy Rcpp::String elements to std::string
					foicols_str.emplace_back(foitmp);
				
					StringVector ldfvec = ldf[x];
					std::vector<std::string> ldftmp(nldf);
					std::copy(ldfvec.cbegin(), ldfvec.cend(), ldftmp.begin()); // Copy Rcpp::String elements to std::string
					ldfcols_str.emplace_back(ldftmp);
				}
			break;
			}
		}

		foiexp.resize(nvar);
		accum = LongLongMatrix(totcats, totcats);
	}
	Accumulate(const Accumulate& acc, Split) : coltype(acc.coltype), foicols_int(acc.foicols_int), ldfcols_int(acc.ldfcols_int), foicols_num(acc.foicols_num), ldfcols_num(acc.ldfcols_num), foicols_lgl(acc.foicols_lgl), ldfcols_lgl(acc.ldfcols_lgl), foicols_str(acc.foicols_str), ldfcols_str(acc.ldfcols_str), intfuncs(acc.intfuncs), realfuncs(acc.realfuncs), strfuncs(acc.strfuncs), nvar(acc.nvar), nfoi(acc.nfoi), nldf(acc.nldf), idbase(acc.idbase), totcats(acc.totcats) {
		foiexp.resize(nvar);
		accum = LongLongMatrix(totcats, totcats);
	}
	// accumulate just the element of the range I've been asked to
	void operator()(std::size_t begin, std::size_t end) {
		if (dbg == true) {
			blockrange.push_back(std::pair<int, int>(begin, end));
		}
		for (std::size_t i = begin; i < end; ++i) {
			for (std::size_t j = 0u; j < nldf; ++j) {
				int int_count = 0;
				int num_count = 0;
				int str_count = 0;
				int lgl_count = 0;
				for (std::size_t x = 0u; x < nvar; ++x) {
					switch (coltype[x]) {
					case 'I':
						foiexp[x] = intfuncs[int_count](foicols_int[int_count][i], ldfcols_int[int_count][j]);
						++int_count;
						break;
					case 'R':
						foiexp[x] = realfuncs[num_count](foicols_num[num_count][i], ldfcols_num[num_count][j]);
						++num_count;
						break;
					case 'L':
						foiexp[x] = foicols_lgl[lgl_count][i] != ldfcols_lgl[lgl_count][j];
						++lgl_count;
						break;
					case 'S':
						foiexp[x] = strfuncs[str_count](foicols_str[str_count][i], ldfcols_str[str_count][j]);
						++str_count;
						break;
					}
				}
				for (std::size_t x = 0u; x < nvar; ++x) {
					for (std::size_t y = 0u; y <= x; ++y) {
						accum(idbase[x] + foiexp[x], idbase[y] + foiexp[y])++;
					}
				}
			}
		}
		// Rcpp::Rcout << "Processed FOI records " << begin << " to " << end << std::endl;
	}
     
	// join my value with that of another Accumulate
	void join(const Accumulate& rhs) {
		if (dbg == true) {
			blockrange.insert(blockrange.end(), rhs.blockrange.begin(), rhs.blockrange.end());
		}
		accum += rhs.accum;
	}
};

//' @title
//' buildAstar
//' @description
//' Builds the A* matrix
//' 
//' @param foinew \code{\link[base]{data.frame}} representing the file of interest
//' 
//' @param ldfnew \code{\link[base]{data.frame}} representing the linking data file
//'
//' @param ncats numeric \code{\link[base]{vector}} indicating the number of categories for each variable
//'
//' @param cutpoints numeric \code{\link[base]{list}} of \code{\link[base]{vector}} indicating the cutpoints to use for for each variable
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
NumericMatrix buildAstar(DataFrame foinew, DataFrame ldfnew, IntegerVector ncats, List cutpoints, int grainsize, bool debug) {
	dbg = debug;
	uint64_t nrecfoi = foinew.nrow();
	uint64_t nrecldf = ldfnew.nrow();
	std::size_t nident = ncats.size();
	uint64_t total = nrecfoi * nrecldf;

	std::vector<uint64_t> idbase(nident); // Start of each variable
	for (std::size_t i = 1; i < nident; ++i) {
		idbase[i] = idbase[i - 1] + ncats[i - 1];
	}

	int totcats = 0;
	for (std::size_t i = 0; i < nident; ++i) {
		totcats += ncats[i];
	}

	NumericMatrix astar(totcats + 2, totcats + 2);

	std::vector<char> coltype;
	for (std::size_t x = 0u; x < nident; ++x) {
		if (TYPEOF(foinew[x]) != TYPEOF(ldfnew[x])) {
			Rcpp::Rcout << "Variable types do not match between LDF and FOI data frames" << std::endl;
			return(astar);
		}
		switch (TYPEOF(foinew[x])) {
		case INTSXP:
			coltype.emplace_back('I');
			if (debug == true) {
				Rcpp::Rcout << "Variable " << x << " is integer" << std::endl;
			}
			break;
		case REALSXP:
			coltype.emplace_back('R');
			if (debug == true) {
				Rcpp::Rcout << "Variable " << x << " is real" << std::endl;
			}
			break;
		case LGLSXP:
			coltype.emplace_back('L');
			if (debug == true) {
				Rcpp::Rcout << "Variable " << x << " is logical" << std::endl;
			}
			break;
		case STRSXP:
			coltype.emplace_back('S');
			if (debug == true) {
				Rcpp::Rcout << "Variable " << x << " is string" << std::endl;
			}
			break;
		default:
			Rcpp::Rcout << "Unsupported variable type" << std::endl;
			return(astar);
		}
	}

	// Define comparison functions
	std::vector<std::function<int(const int, const int)> > intfuncs;
	std::vector<std::function<int(const double, const double)> > realfuncs;
	std::vector<std::function<int(const std::string&, const std::string&)> > strfuncs;
	for (std::size_t x = 0u; x < nident; ++x) {
		NumericVector cutstmp = cutpoints[x];
		std::vector<double> cuts(cutstmp.begin(), cutstmp.end());
		std::size_t numcats = ncats[x];
		switch (coltype[x]) {
		case 'I':
			if (numcats == 2) {
				intfuncs.emplace_back([](const int foi, const int ldf)->int { return ldf != foi; });
			}
			else {
				intfuncs.emplace_back([cuts](const int foi, const int ldf)->int {
					size_t cat = cuts.size();
					for (std::size_t i = 0; i < cuts.size(); ++i) {
						if (std::abs(foi - ldf) < cuts[i]) {
							break;
						}
						--cat;
					}
					return cat;
				});
			}
			break;
		case 'R': 
			if (numcats == 2) {
				realfuncs.emplace_back([](const double foi, const double ldf)->int { return ldf != foi; });
			}
			else {
				realfuncs.emplace_back([cuts](const double foi, const double ldf)->int {
					size_t cat = cuts.size();
					for (std::size_t i = 0; i < cuts.size(); ++i) {
						if (std::abs(foi - ldf) < cuts[i]) {
							break;
						}
						--cat;
					}
					return cat;
					return 0;
				});
			}
			break;
		case 'L':
			break;
		case 'S':
			if (numcats == 2) {
				strfuncs.emplace_back([](const std::string &foi, const std::string &ldf)->int { return ldf != foi; });
			}
			else {
				strfuncs.emplace_back([cuts](const std::string &foi, const std::string &ldf)->int {
					std::set<std::string> bigrams1;
					std::set<std::string> bigrams2;
					if (foi.length() != 0 && ldf.length() != 0) {
						for (std::size_t i = 0; i < foi.length() - 1; ++i) {
							bigrams1.insert(foi.substr(i, 2));
						}
						for (std::size_t i = 0; i < ldf.length() - 1; ++i) {
							bigrams2.insert(ldf.substr(i, 2));
						}
						std::size_t common = 0;
						for (const auto& val : bigrams2) {
							common += bigrams1.count(val);
						}
						std::size_t total = bigrams1.size() + bigrams2.size();
						double dice = static_cast<double>(common * 2) / static_cast<double>(total);
						
						size_t cat = cuts.size();
						for (std::size_t i = 0; i < cuts.size(); ++i) {
							if (dice > cuts[cat - 1]) {
								break;
							}
							--cat;
						}
						return cat;
					}
					return 0;
				});
			}
			break;
		}
	}

	LongLongMatrix accum(totcats, totcats);
	if (grainsize < nrecfoi) {
		// declare the Accumulate instance 
		Accumulate acc(foinew, ldfnew, coltype, intfuncs, realfuncs, strfuncs, idbase, totcats);
		// call parallel_reduce to start the work
		parallelReduce(0, nrecfoi, acc, grainsize);
		if (debug == true) {
			Rcpp::Rcout << "Blocks processed:" << std::endl;
			for (std::size_t i = 0u; i < acc.blockrange.size(); ++i) {
				Rcpp::Rcout << "Range: " << acc.blockrange[i].first << " to " << acc.blockrange[i].second << std::endl;
			}
		}
		accum += acc.accum;
	} else {
		// run the calculations in series
		std::size_t nreport = 10;
		std::size_t step = nrecfoi / 100 * nreport;
		std::size_t nextstep = step;
		std::size_t cent = nreport;

		std::vector<int> foiexp(nident);

		std::vector<IntegerVector> foicols_int;
		std::vector<IntegerVector> ldfcols_int;
		std::vector<NumericVector> foicols_num;
		std::vector<NumericVector> ldfcols_num;
		std::vector<LogicalVector> foicols_lgl;
		std::vector<LogicalVector> ldfcols_lgl;
		//std::vector<StringVector> foicols_str;
		//std::vector<StringVector> ldfcols_str;
		std::vector<std::vector<std::string> > foicols_str;
		std::vector<std::vector<std::string> > ldfcols_str;
		
		for (std::size_t x = 0u; x < nident; ++x) {
			switch (coltype[x]) {
			case 'I':
				foicols_int.emplace_back(foinew[x]);
				ldfcols_int.emplace_back(ldfnew[x]);
				break;
			case 'R': 
				foicols_num.emplace_back(foinew[x]);
				ldfcols_num.emplace_back(ldfnew[x]);
				break;
			case 'L':
				foicols_lgl.emplace_back(foinew[x]);
				ldfcols_lgl.emplace_back(ldfnew[x]);
				break;
			case 'S':
				//foicols_str.emplace_back(foinew[x]);
				//ldfcols_str.emplace_back(ldfnew[x]);

				// Copy contents of Rcpp::StringVector to std::vector<std::string>, as this saves converting individual elements later when comparing
				// NOTE: This will approximately double the storage space required
				{
					StringVector foivec = foinew[x];
					std::vector<std::string> foitmp(nrecfoi);
					std::copy(foivec.cbegin(), foivec.cend(), foitmp.begin()); // Copy Rcpp::String elements to std::string
					foicols_str.emplace_back(foitmp);

					StringVector ldfvec = ldfnew[x];
					std::vector<std::string> ldftmp(nrecldf);
					std::copy(ldfvec.cbegin(), ldfvec.cend(), ldftmp.begin()); // Copy Rcpp::String elements to std::string
					ldfcols_str.emplace_back(ldftmp);
				}
				break;
			}
		}

		Rcpp::Rcout << "Progress (%): |0";
		std::flush(Rcpp::Rcout);
		for (std::size_t i = 0u; i < nrecfoi; ++i) {
			for (std::size_t j = 0u; j < nrecldf; ++j) {
				int int_count = 0;
				int num_count = 0;
				int str_count = 0;
				int lgl_count = 0;
				for (std::size_t x = 0u; x < foiexp.size(); ++x) {
					switch (coltype[x]) {
					case 'I':
						foiexp[x] = intfuncs[int_count](foicols_int[int_count][i], ldfcols_int[int_count][j]);
						++int_count;
						break;
					case 'R':
						foiexp[x] = realfuncs[num_count](foicols_num[num_count][i], ldfcols_num[num_count][j]);
						++num_count;
						break;
					case 'L':
						foiexp[x] = foicols_lgl[lgl_count][i] != ldfcols_lgl[lgl_count][j];
						++lgl_count;
						break;
					case 'S':
						//foiexp[x] = strfuncs[str_count](Rcpp::as<std::string>(foicols_str[str_count][i]), Rcpp::as<std::string>(ldfcols_str[str_count][j]));
						foiexp[x] = strfuncs[str_count](foicols_str[str_count][i], ldfcols_str[str_count][j]);
						++str_count;
						break;
					}
					
				}
				for (std::size_t x = 0u; x < foiexp.size(); ++x) {
					for (std::size_t y = 0u; y <= x; ++y) {
						accum(idbase[x] + foiexp[x], idbase[y] + foiexp[y])++;
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
		for (std::size_t i = 0u; i < accum.nrow; i++) {
			for (std::size_t j = 0u; j < accum.ncol; j++) {
				Rcpp::Rcout << accum(i, j) << " ";
			}
			Rcpp::Rcout << std::endl;
		}
	}

	// Test that sums add up
	bool sumerr = false;
	for (std::size_t i = 0u; i < nident; ++i) {
		for (std::size_t j = 0u; j < i; ++j) {
			uint64_t test = total;
			for (int x = 0; x < ncats[i]; ++x) {
				for (int y = 0; y < ncats[j]; ++y) {
					test -= accum(idbase[i] + x, idbase[j] + y);
				}
			}
			if (test != 0) {
				sumerr = true;
			}
		}
	}

	for (int i = 0; i < totcats; ++i) {
		astar(i, i) = accum(i, i) * (nident - 1.0) / static_cast<double>(nident * nident) * 2.0;
		for (int j = 0; j < i; j++) {
			astar(i, j) -= accum(i, j) / static_cast<double>(nident * nident) * 2.0;
		}
	}

	// Fill in other half of symmetric matrix
	for (int i = 0; i < totcats; ++i) {
		for (int j = i + 1; j < totcats; j++) {
			astar(i, j) = astar(j, i);
		}
	}

	// Add q and r

	// -q
	for (std::size_t i = 0u; i < nident; ++i) {
		astar(idbase[i], totcats) -= 1.0 / nident;
		astar(idbase[i], totcats) *= total;
	}

	// -r
	for (std::size_t i = 0u; i < nident; ++i) {
		astar(idbase[i] + ncats[i] - 1, totcats + 1) -= 1.0 / nident;
		astar(idbase[i] + ncats[i] - 1, totcats + 1) *= total;
	}

	// q
	for (std::size_t i = 0u; i < nident; ++i) {
		astar(totcats, idbase[i]) = 1.0 / nident;
		astar(totcats, idbase[i]) *= total;
	}

	// r
	for (std::size_t i = 0u; i < nident; ++i) {
		astar(totcats + 1, idbase[i] + ncats[i] - 1) = 1.0 / nident;
		astar(totcats + 1, idbase[i] + ncats[i] - 1) *= total;
	}

	if (sumerr == true) {
		Rcpp::Rcout << "Warning: unexpected sum (should be " << total << "), possible overflow" << std::endl;
	}

	return(astar);
}