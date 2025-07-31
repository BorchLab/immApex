// ── src/fastEditEdges.cpp ───────────────────────────────────────────────────
#include <Rcpp.h>
#include <numeric>
#include <cmath>          
#include <string>
#include <vector>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// ----------------------------------------------------------------------------
//  Levenshtein with early‐exit (banded DP)
// ----------------------------------------------------------------------------
static inline int lv_threshold(const std::string &a,
                               const std::string &b,
                               int                thr)
{
  int n = a.size(), m = b.size();
  if (std::abs(n - m) > thr) return thr + 1;
  
  std::vector<int> row(m+1);
  std::iota(row.begin(), row.end(), 0);
  
  for (int i = 1; i <= n; ++i) {
    int prev = row[0];
    row[0]   = i;
    int row_min = i;
    
    for (int j = 1; j <= m; ++j) {
      int cur = row[j];
      int cost = (a[i-1] == b[j-1] ? 0 : 1);
      row[j] = std::min({ row[j-1] + 1,
                        row[j]   + 1,
                        prev     + cost });
      prev = cur;
      row_min = std::min(row_min, row[j]);
    }
    if (row_min > thr) return thr + 1;
  }
  return row[m];
}

// ----------------------------------------------------------------------------
//  Main export 
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
DataFrame fast_edge_list(CharacterVector           seqs,
                         double                    thresh   = 1.0,
                         Nullable<CharacterVector> v_gene   = R_NilValue,
                         Nullable<CharacterVector> j_gene   = R_NilValue,
                         bool                      match_v  = false,
                         bool                      match_j  = false,
                         Nullable<CharacterVector> ids      = R_NilValue)
{
  int n = seqs.size();
  if (n < 2) stop("At least two sequences are required.");
  if (thresh < 0) stop("`thresh` must be >= 0.");
  
  // 1) Read sequences + lengths
  std::vector<std::string> s(n);
  std::vector<int>         lens(n);
  for (int i = 0; i < n; ++i) {
    s[i]    = as<std::string>(seqs[i]);
    lens[i] = s[i].size();
  }
  
  // 2) Optional V/J
  std::vector<std::string> v(n), J(n);
  if (match_v) {
    if (v_gene.isNull()) stop("match_v=TRUE but no v_gene provided.");
    CharacterVector vv(v_gene);
    if (vv.size()==1) vv = CharacterVector(n, vv[0]);
    if ((int)vv.size() != n) stop("Length of v_gene != number of sequences.");
    v = as< std::vector<std::string> >(vv);
  }
  if (match_j) {
    if (j_gene.isNull()) stop("match_j=TRUE but no j_gene provided.");
    CharacterVector jj(j_gene);
    if (jj.size()==1) jj = CharacterVector(n, jj[0]);
    if ((int)jj.size() != n) stop("Length of j_gene != number of sequences.");
    J = as< std::vector<std::string> >(jj);
  }
  
  // 3) IDs / labels
  std::vector<std::string> lbl(n);
  if (ids.isNull()) {
    for (int i = 0; i < n; ++i) lbl[i] = std::to_string(i+1);
  } else {
    CharacterVector ix(ids);
    if ((int)ix.size() != n) stop("Length of ids != number of sequences.");
    lbl = as< std::vector<std::string> >(ix);
  }
  
  // 4) Reserve output buffers
  size_t approx = (size_t)n * (n-1) / 2;
  std::vector<std::string> out_from; out_from.reserve(approx);
  std::vector<std::string> out_to;   out_to  .reserve(approx);
  std::vector<double>      out_d;    out_d   .reserve(approx);
  
#ifdef _OPENMP
  int nthread = omp_get_max_threads();
#else
  int nthread = 1;
#endif
  
#ifdef _OPENMP
#pragma omp parallel
#endif
{
  std::vector<std::string> loc_from; loc_from.reserve(approx / nthread);
  std::vector<std::string> loc_to;   loc_to  .reserve(approx / nthread);
  std::vector<double>      loc_d;    loc_d   .reserve(approx / nthread);
  
#ifdef _OPENMP
#pragma omp for schedule(dynamic,32)
#endif
  for (int i = 0; i < n; ++i) {
    for (int k = i + 1; k < n; ++k) {
      if (match_v && v[i] != v[k]) continue;
      if (match_j && J[i] != J[k]) continue;
      
      int max_len = std::max(lens[i], lens[k]);
      if (max_len == 0) continue; // Skip empty sequences
      
      int maxd; // Maximum allowed absolute Levenshtein distance
      
      if (thresh >= 1.0) {
        maxd = std::lround(thresh);
      } else {
        double dist_thresh_rel = 1.0 - thresh;
        maxd = std::floor(dist_thresh_rel * max_len);
      }
      
      // if length difference > max allowed distance, skip calculation.
      if (std::abs(lens[i] - lens[k]) > maxd) continue;
      
      // Calculate Levenshtein distance with early-exit optimization.
      int d = lv_threshold(s[i], s[k], maxd);
      
      // Filter: Return only edges that meet the threshold criteria.
      if (d <= maxd) {
        loc_from.push_back(lbl[i]);
        loc_to.push_back(lbl[k]);
        
        // Store the appropriate distance value in the output data frame.
        if (thresh >= 1.0) {
          loc_d.push_back(static_cast<double>(d)); // Store absolute distance
        } else {
          double rel_dist = static_cast<double>(d) / max_len;
          loc_d.push_back(rel_dist); // Store relative distance
        }
      }
    }
  }
  
#ifdef _OPENMP
#pragma omp critical
#endif
{
  out_from.insert(out_from.end(), loc_from.begin(), loc_from.end());
  out_to  .insert(out_to  .end(),   loc_to  .begin(),   loc_to  .end());
  out_d   .insert(out_d   .end(),   loc_d   .begin(),   loc_d   .end());
}
}

return DataFrame::create(
  _["from"]               = out_from,
  _["to"]                 = out_to,
  _["dist"]               = out_d,
  _["stringsAsFactors"]   = false
);
}