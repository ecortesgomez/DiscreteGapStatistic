#include <Rcpp.h>
#include <unordered_map>
#include <cmath>
using namespace Rcpp;

//' Bhattacharyya Distance
//'
//' @name BhattacharyyaDist_rcpp
//' @title Bhattacharyya Distance via Rcpp
//' @description Computes pairwise Bhattacharyya distance between rows.
//' @param x n X p character matrix.
//' @param offset small offset for log(0*0) cases.
//' @return Distance matrix between rows.
//' @export

// [[Rcpp::export]]
NumericMatrix BhattacharyyaDist_rcpp(CharacterMatrix x, double offset = 1e-8) {
   int n = x.nrow();
   int p = x.ncol();

   // --- Step 1: Identify unique categories ---
   std::unordered_map<std::string, int> cat_index;
   int cat_counter = 0;
   for (int i = 0; i < n * p; ++i) {
      std::string s = Rcpp::as<std::string>(x[i]);
      if (cat_index.find(s) == cat_index.end()) {
         cat_index[s] = cat_counter++;
      }
   }
   int m = cat_counter; // number of unique categories

   // --- Step 2: Build frequency matrix ---
   NumericMatrix freq(n, m);
   for (int i = 0; i < n; ++i) {
      for (int j = 0; j < p; ++j) {
         std::string s = Rcpp::as<std::string>(x(i, j));
         freq(i, cat_index[s]) += 1.0;
      }
   }

   // --- Step 3: Normalize to probabilities ---
   for (int i = 0; i < n; ++i) {
      double sum_row = 0.0;
      for (int k = 0; k < m; ++k) sum_row += freq(i, k);
      if (sum_row > 0) {
         for (int k = 0; k < m; ++k)
            freq(i, k) /= sum_row;
      }
   }

   // --- Step 4: Compute pairwise Bhattacharyya distances ---
   NumericMatrix dist(n, n);
   for (int i = 0; i < n; ++i) {
      dist(i, i) = 0.0;
      for (int j = i + 1; j < n; ++j) {
         double bc = 0.0;
         for (int k = 0; k < m; ++k)
            bc += std::sqrt(freq(i, k) * freq(j, k));

         // Add offset *after* summation
         double d = -std::log(bc + offset);

         dist(i, j) = d;
         dist(j, i) = d;
      }
   }

   return dist;
}

//' Chi-square Distance
//'
//' @name ChisqDist_rcpp
//' @title Chi-square Distance via Rcpp
//' @description Computes pairwise Chi-square distance between rows.
//' @param x n X p character matrix.
//' @return Distance matrix between rows.
//' @export

// [[Rcpp::export]]
NumericMatrix ChisqDist_rcpp(CharacterMatrix x) {
   int n = x.nrow();
   int p = x.ncol();

   // --- Determine unique categories across entire matrix ---
   std::unordered_map<std::string, int> cat_map;
   int cat_index = 0;
   for (int i = 0; i < n; ++i) {
      for (int j = 0; j < p; ++j) {
         std::string val = Rcpp::as<std::string>(x(i, j));
         if (!cat_map.count(val)) {
            cat_map[val] = cat_index++;
         }
      }
   }
   int m = cat_index;  // number of unique categories

   // --- Build frequency matrix (n x m) ---
   NumericMatrix freq(n, m);
   for (int i = 0; i < n; ++i) {
      for (int j = 0; j < p; ++j) {
         std::string val = Rcpp::as<std::string>(x(i, j));
         int k = cat_map[val];
         freq(i, k) += 1.0;
      }
   }

   // --- Compute pairwise distances ---
   NumericMatrix D(n, n);
   double two_over_p = 2.0 / p;

   for (int i = 0; i < n; ++i) {
      for (int j = i + 1; j < n; ++j) {
         double dist = 0.0;
         for (int g = 0; g < m; ++g) {
            double fx = freq(i, g);
            double fy = freq(j, g);
            double denom = fx + fy;
            if (denom > 0.0) {
               double num = fx - fy;
               dist += (num * num) / denom;
            }
         }
         dist *= two_over_p;
         D(i, j) = dist;
         D(j, i) = dist;
      }
   }

   return D;
}
