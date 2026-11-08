#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>

using namespace Rcpp;

static std::string get_string_value(const CharacterMatrix& x, int i, int j) {
  Rcpp::String s = x(i, j);
  if (s == NA_STRING) {
    Rcpp::stop("Missing values are not supported by these distance functions.");
  }
  return std::string(s);
}

static int get_int_value(const IntegerMatrix& x, int i, int j) {
  int v = x(i, j);
  if (v == NA_INTEGER) {
    Rcpp::stop("Missing values are not supported by these distance functions.");
  }
  return v;
}

static std::unordered_map<std::string, int> make_category_map(const CharacterMatrix& x) {
  int n = x.nrow();
  int p = x.ncol();
  std::unordered_map<std::string, int> cat_map;
  int next_id = 0;

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < p; ++j) {
      std::string val = get_string_value(x, i, j);
      if (!cat_map.count(val)) {
        cat_map[val] = next_id++;
      }
    }
  }
  return cat_map;
}

static std::vector< std::vector<double> > make_frequency_matrix(
    const CharacterMatrix& x,
    const std::unordered_map<std::string, int>& cat_map) {

  int n = x.nrow();
  int p = x.ncol();
  int m = cat_map.size();
  std::vector< std::vector<double> > freq(n, std::vector<double>(m, 0.0));

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < p; ++j) {
      std::string val = get_string_value(x, i, j);
      int k = cat_map.at(val);
      freq[i][k] += 1.0;
    }
  }
  return freq;
}

static std::vector<int> sorted_row(const IntegerMatrix& x, int row) {
  int p = x.ncol();
  std::vector<int> out(p);
  for (int j = 0; j < p; ++j) out[j] = get_int_value(x, row, j);
  std::sort(out.begin(), out.end());
  return out;
}

static double ecdf_at_sorted_row(const std::vector<int>& sorted_values, int value) {
  return static_cast<double>(std::upper_bound(sorted_values.begin(), sorted_values.end(), value) - sorted_values.begin()) /
    static_cast<double>(sorted_values.size());
}

static std::vector<int> union_levels(const std::vector<int>& a, const std::vector<int>& b) {
  std::vector<int> levels;
  levels.reserve(a.size() + b.size());
  levels.insert(levels.end(), a.begin(), a.end());
  levels.insert(levels.end(), b.begin(), b.end());
  std::sort(levels.begin(), levels.end());
  levels.erase(std::unique(levels.begin(), levels.end()), levels.end());
  return levels;
}

static std::vector<double> sorted_average_ranks(const IntegerMatrix& x, int row) {
  int p = x.ncol();
  std::vector<int> vals(p);
  for (int j = 0; j < p; ++j) vals[j] = get_int_value(x, row, j);
  std::sort(vals.begin(), vals.end());

  std::vector<double> ranks;
  ranks.reserve(p);
  int start = 0;
  while (start < p) {
    int end = start;
    while (end + 1 < p && vals[end + 1] == vals[start]) ++end;

    // R's default rank() uses average ranks with 1-based positions.
    double avg_rank = (static_cast<double>(start + 1) + static_cast<double>(end + 1)) / 2.0;
    for (int k = start; k <= end; ++k) ranks.push_back(avg_rank);
    start = end + 1;
  }
  return ranks;
}

//' Hamming distance
//'
//' @name HammingDist_rcpp
//' @title Hamming distance via Rcpp
//' @description Computes the adjusted Hamming distance between rows: mismatches divided by the number of columns.
//' @param x n x p character matrix.
//' @return Numeric n x n distance matrix.
//' @export
// [[Rcpp::export]]
NumericMatrix HammingDist_rcpp(CharacterMatrix x) {
  int n = x.nrow();
  int p = x.ncol();
  NumericMatrix D(n, n);

  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      int mismatches = 0;
      for (int q = 0; q < p; ++q) {
        if (get_string_value(x, i, q) != get_string_value(x, j, q)) ++mismatches;
      }
      double d = static_cast<double>(mismatches) / static_cast<double>(p);
      D(i, j) = d;
      D(j, i) = d;
    }
  }
  return D;
}

//' Bhattacharyya distance
//'
//' @name BhattacharyyaDist_rcpp
//' @title Bhattacharyya distance via Rcpp
//' @description Computes pairwise Bhattacharyya distances between row-wise category frequency profiles.
//' @param x n x p character matrix.
//' @param offset small positive offset added before log; use 0.01 for strict equivalence with the original R fallback.
//' @return Numeric n x n distance matrix.
//' @export
// [[Rcpp::export]]
NumericMatrix BhattacharyyaDist_rcpp(CharacterMatrix x, double offset = 0.01) {
  int n = x.nrow();
  int p = x.ncol();
  std::unordered_map<std::string, int> cat_map = make_category_map(x);
  std::vector< std::vector<double> > freq = make_frequency_matrix(x, cat_map);
  int m = cat_map.size();

  // Original R fallback uses frequencies divided by the number of columns.
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < m; ++k) freq[i][k] /= static_cast<double>(p);
  }

  NumericMatrix D(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      double bc = 0.0;
      for (int k = 0; k < m; ++k) {
        bc += std::sqrt(freq[i][k] * freq[j][k]);
      }
      double d = -std::log(bc + offset);
      if (d < 0.0) d = 0.0;
      D(i, j) = d;
      D(j, i) = d;
    }
  }
  return D;
}

//' Chi-square distance
//'
//' @name ChisqDist_rcpp
//' @title Chi-square distance via Rcpp
//' @description Computes pairwise chi-square distances between row-wise category counts, divided by the number of columns.
//' @param x n x p character matrix.
//' @return Numeric n x n distance matrix.
//' @export
// [[Rcpp::export]]
NumericMatrix ChisqDist_rcpp(CharacterMatrix x) {
  int n = x.nrow();
  int p = x.ncol();
  std::unordered_map<std::string, int> cat_map = make_category_map(x);
  std::vector< std::vector<double> > freq = make_frequency_matrix(x, cat_map);
  int m = cat_map.size();

  NumericMatrix D(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      double dist = 0.0;
      for (int g = 0; g < m; ++g) {
        double fx = freq[i][g];
        double fy = freq[j][g];
        double denom = fx + fy;
        if (denom > 0.0) {
          double diff = fx - fy;
          dist += (diff * diff) / denom;
        }
      }
      dist /= static_cast<double>(p);  // matches .ChisqDist_R(); do not multiply by 2.
      D(i, j) = dist;
      D(j, i) = dist;
    }
  }
  return D;
}

//' Hellinger distance
//'
//' @name HellingerDist_rcpp
//' @title Hellinger distance via Rcpp
//' @description Computes pairwise Hellinger distances between row-wise category frequency profiles.
//' @param x n x p character matrix.
//' @return Numeric n x n distance matrix.
//' @export
// [[Rcpp::export]]
NumericMatrix HellingerDist_rcpp(CharacterMatrix x) {
  int n = x.nrow();
  int p = x.ncol();
  std::unordered_map<std::string, int> cat_map = make_category_map(x);
  std::vector< std::vector<double> > freq = make_frequency_matrix(x, cat_map);
  int m = cat_map.size();

  NumericMatrix D(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      double ss = 0.0;
      for (int g = 0; g < m; ++g) {
        double diff = std::sqrt(freq[i][g]) - std::sqrt(freq[j][g]);
        ss += diff * diff;
      }
      double d = std::sqrt(ss / (2.0 * static_cast<double>(p)));
      D(i, j) = d;
      D(j, i) = d;
    }
  }
  return D;
}

//' Cramer's V distance
//'
//' @name CramerVDist_rcpp
//' @title Cramer's V distance via Rcpp
//' @description Computes 1 - Cramer's V between each pair of rows using the same contingency-table formula as stats::chisq.test(x, y, correct = FALSE).
//' @param x n x p character matrix.
//' @return Numeric n x n distance matrix.
//' @export
// [[Rcpp::export]]
NumericMatrix CramerVDist_rcpp(CharacterMatrix x) {
  int n = x.nrow();
  int p = x.ncol();
  NumericMatrix D(n, n);

  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      bool identical_rows = true;
      std::unordered_map<std::string, int> row_map;
      std::unordered_map<std::string, int> col_map;
      int r_next = 0;
      int c_next = 0;

      for (int q = 0; q < p; ++q) {
        std::string a = get_string_value(x, i, q);
        std::string b = get_string_value(x, j, q);
        if (a != b) identical_rows = false;
        if (!row_map.count(a)) row_map[a] = r_next++;
        if (!col_map.count(b)) col_map[b] = c_next++;
      }

      double d = 0.0;
      if (identical_rows) {
        d = 0.0;  // V = 1, distance = 0
      } else if (r_next == 1 || c_next == 1) {
        d = 1.0;  // original fallback returns V = 0, distance = 1
      } else {
        std::vector< std::vector<double> > obs(r_next, std::vector<double>(c_next, 0.0));
        std::vector<double> row_sum(r_next, 0.0), col_sum(c_next, 0.0);

        for (int q = 0; q < p; ++q) {
          std::string a = get_string_value(x, i, q);
          std::string b = get_string_value(x, j, q);
          int rr = row_map[a];
          int cc = col_map[b];
          obs[rr][cc] += 1.0;
          row_sum[rr] += 1.0;
          col_sum[cc] += 1.0;
        }

        double chi2 = 0.0;
        double N = static_cast<double>(p);
        for (int rr = 0; rr < r_next; ++rr) {
          for (int cc = 0; cc < c_next; ++cc) {
            double expected = row_sum[rr] * col_sum[cc] / N;
            if (expected > 0.0) {
              double diff = obs[rr][cc] - expected;
              chi2 += (diff * diff) / expected;
            }
          }
        }

        int k = std::min(r_next, c_next);
        double V = std::sqrt(chi2 / (N * static_cast<double>(k - 1)));
        if (V < 0.0) V = 0.0;
        if (V > 1.0) V = 1.0;
        d = 1.0 - V;
      }

      D(i, j) = d;
      D(j, i) = d;
    }
  }
  return D;
}

//' Kolmogorov-Smirnov ordinal distance
//'
//' @name KSDist_rcpp
//' @title Kolmogorov-Smirnov ordinal distance via Rcpp
//' @description Computes row-wise KS distance by evaluating empirical CDFs on the union of observed row levels.
//' @param x n x p integer matrix.
//' @return Numeric n x n distance matrix.
//' @export
// [[Rcpp::export]]
NumericMatrix KSDist_rcpp(IntegerMatrix x) {
  int n = x.nrow();
  NumericMatrix D(n, n);
  std::vector< std::vector<int> > rows(n);
  for (int i = 0; i < n; ++i) rows[i] = sorted_row(x, i);

  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      std::vector<int> levels = union_levels(rows[i], rows[j]);
      double maxdiff = 0.0;
      for (std::size_t g = 0; g < levels.size(); ++g) {
        double Fx = ecdf_at_sorted_row(rows[i], levels[g]);
        double Fy = ecdf_at_sorted_row(rows[j], levels[g]);
        double diff = std::fabs(Fx - Fy);
        if (diff > maxdiff) maxdiff = diff;
      }
      D(i, j) = maxdiff;
      D(j, i) = maxdiff;
    }
  }
  return D;
}

//' Wasserstein DS ordinal distance
//'
//' @name WassersteinDS_rcpp
//' @title Data-support Wasserstein ordinal distance via Rcpp
//' @description Computes the row-wise ECDF distance on the union of levels observed in each row pair.
//' @param x n x p integer matrix.
//' @param power positive power. The current package default is 1.
//' @return Numeric n x n distance matrix.
//' @export
// [[Rcpp::export]]
NumericMatrix WassersteinDS_rcpp(IntegerMatrix x, double power = 1.0) {
  if (power <= 0.0) Rcpp::stop("'power' must be positive.");
  int n = x.nrow();
  NumericMatrix D(n, n);
  std::vector< std::vector<int> > rows(n);
  for (int i = 0; i < n; ++i) rows[i] = sorted_row(x, i);

  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      std::vector<int> levels = union_levels(rows[i], rows[j]);
      double acc = 0.0;
      for (std::size_t g = 0; g < levels.size(); ++g) {
        double Fx = ecdf_at_sorted_row(rows[i], levels[g]);
        double Fy = ecdf_at_sorted_row(rows[j], levels[g]);
        acc += std::pow(std::fabs(Fx - Fy), power);
      }
      double d = std::pow(acc / static_cast<double>(levels.size()), 1.0 / power);
      D(i, j) = d;
      D(j, i) = d;
    }
  }
  return D;
}

//' Wasserstein KS ordinal distance
//'
//' @name WassersteinKS_rcpp
//' @title Known-support Wasserstein ordinal distance via Rcpp
//' @description Computes the row-wise ECDF distance on the fixed support 1:maxM.
//' @param x n x p integer matrix.
//' @param maxM maximum integer support value.
//' @param power positive power. The current package default is 1.
//' @return Numeric n x n distance matrix.
//' @export
// [[Rcpp::export]]
NumericMatrix WassersteinKS_rcpp(IntegerMatrix x, int maxM, double power = 1.0) {
  if (power <= 0.0) Rcpp::stop("'power' must be positive.");
  if (maxM < 1) Rcpp::stop("'maxM' must be a positive integer.");

  int n = x.nrow();
  NumericMatrix D(n, n);
  std::vector< std::vector<int> > rows(n);
  for (int i = 0; i < n; ++i) rows[i] = sorted_row(x, i);

  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      double acc = 0.0;
      for (int level = 1; level <= maxM; ++level) {
        double Fx = ecdf_at_sorted_row(rows[i], level);
        double Fy = ecdf_at_sorted_row(rows[j], level);
        acc += std::pow(std::fabs(Fx - Fy), power);
      }
      double d = std::pow(acc / static_cast<double>(maxM), 1.0 / power);
      D(i, j) = d;
      D(j, i) = d;
    }
  }
  return D;
}

//' Wasserstein rank ordinal distance
//'
//' @name WassersteinRank_rcpp
//' @title Rank Wasserstein ordinal distance via Rcpp
//' @description Computes the row-wise distance between sorted average-rank vectors, matching sort(rank(.)) in R.
//' @param x n x p integer matrix.
//' @param power positive power. The current package default is 1.
//' @return Numeric n x n distance matrix.
//' @export
// [[Rcpp::export]]
NumericMatrix WassersteinRank_rcpp(IntegerMatrix x, double power = 1.0) {
  if (power <= 0.0) Rcpp::stop("'power' must be positive.");

  int n = x.nrow();
  int p = x.ncol();
  NumericMatrix D(n, n);
  std::vector< std::vector<double> > ranks(n);
  for (int i = 0; i < n; ++i) ranks[i] = sorted_average_ranks(x, i);

  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      double acc = 0.0;
      for (int q = 0; q < p; ++q) {
        acc += std::pow(std::fabs(ranks[i][q] - ranks[j][q]), power);
      }
      double d = std::pow(acc / static_cast<double>(p), 1.0 / power);
      D(i, j) = d;
      D(j, i) = d;
    }
  }
  return D;
}

// Podani's (1999) ordinal dissimilarity for integer-coded categories.
// Integer codes determine category order, not numerical spacing.
//
// For two different observed levels a < b within a column:
//
//   partial distance = (1 + count(a < X < b)) /
//                      (n_observed - count(min) - count(max) + 1).
//
// This is Podani's midrank/tie-adjusted expression, simplified using
// cumulative frequencies. Equal levels contribute zero. Constant columns
// contribute zero and remain in the average. Missing values are omitted
// pairwise; pairs with no shared observed columns receive NA.
//
// The returned R 'dist' vector uses lower-triangular, column-wise order.

//' Ordinal Podani distance via Rcpp
//'
//' @name PodaniDist_rcpp
//' @title Ordinal Podani distance via Rcpp
//' @description Computes Podani's rank-and-tie-adjusted dissimilarity
//'   for each pair of rows of an integer-coded ordinal matrix. Numeric
//'   spacing between category codes is ignored; only ordering matters.
//' @param x An n x p integer matrix, with rows as observations and
//'   columns as ordinal variables. Missing values (NA) are permitted.
//' @return A compact object of class \code{dist}, with equal weights over
//'   available variables. All-missing pairs have distance NA.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector PodaniDist_rcpp(Rcpp::IntegerMatrix x) {
    const int n = x.nrow();
    const int p = x.ncol();

    if (n < 2 || p < 1) {
        Rcpp::stop("x must have at least two rows and one column.");
    }

    // Column-major storage: category code for row i, column j at j*n + i.
    // Code -1 designates a missing observation.
    std::vector<int> codes(static_cast<size_t>(n) * p, -1);
    std::vector<std::vector<int> > cumulative(p);
    std::vector<double> denominators(p, 0.0);

    for (int j = 0; j < p; ++j) {
        std::vector<int> observed;
        observed.reserve(n);

        for (int i = 0; i < n; ++i) {
            const int value = x(i, j);
            if (value != NA_INTEGER) {
                observed.push_back(value);
            }
        }

        if (observed.empty()) {
            continue;
        }

        std::sort(observed.begin(), observed.end());
        std::vector<int> levels = observed;
        levels.erase(std::unique(levels.begin(), levels.end()), levels.end());

        const int n_levels = static_cast<int>(levels.size());
        std::vector<int> counts(n_levels, 0);

        for (size_t k = 0; k < observed.size(); ++k) {
            const int level = static_cast<int>(
                std::lower_bound(levels.begin(), levels.end(), observed[k]) -
                levels.begin()
            );
            ++counts[level];
        }

        cumulative[j].assign(n_levels + 1, 0);
        for (int k = 0; k < n_levels; ++k) {
            cumulative[j][k + 1] = cumulative[j][k] + counts[k];
        }

        if (n_levels > 1) {
            denominators[j] = static_cast<double>(
                observed.size() - counts.front() - counts.back() + 1
            );
        }

        for (int i = 0; i < n; ++i) {
            const int value = x(i, j);
            if (value == NA_INTEGER) {
                continue;
            }

            codes[static_cast<size_t>(j) * n + i] = static_cast<int>(
                std::lower_bound(levels.begin(), levels.end(), value) -
                levels.begin()
            );
        }
    }

    // R 'dist' ordering: (2,1), (3,1), ..., (n,1), (3,2), ...
    const R_xlen_t pairs = static_cast<R_xlen_t>(n) * (n - 1) / 2;
    Rcpp::NumericVector result(pairs);
    R_xlen_t position = 0;

    for (int i = 0; i < n - 1; ++i) {
        if (i % 256 == 0) {
            Rcpp::checkUserInterrupt();
        }

        for (int k = i + 1; k < n; ++k) {
            double sum = 0.0;
            int available = 0;

            for (int j = 0; j < p; ++j) {
                const int code_i = codes[static_cast<size_t>(j) * n + i];
                const int code_k = codes[static_cast<size_t>(j) * n + k];

                if (code_i < 0 || code_k < 0) {
                    continue;
                }

                ++available;
                if (code_i == code_k) {
                    continue;
                }

                const int lower = std::min(code_i, code_k);
                const int upper = std::max(code_i, code_k);
                const int between = cumulative[j][upper] - cumulative[j][lower + 1];

                sum += (1.0 + static_cast<double>(between)) / denominators[j];
            }

            result[position++] = available > 0 ? sum / available : NA_REAL;
        }
    }

    result.attr("Size") = n;
    result.attr("Diag") = false;
    result.attr("Upper") = false;
    result.attr("method") = "podani";
    result.attr("class") = "dist";

    if (x.hasAttribute("dimnames")) {
        Rcpp::List names = x.attr("dimnames");
        if (!Rf_isNull(names[0])) {
            result.attr("Labels") = names[0];
        }
    }

    return result;
}
