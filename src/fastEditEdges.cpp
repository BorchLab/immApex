// ── src/fastEditEdges.cpp ───────────────────────────────────────────────────
#include <Rcpp.h>
#include <string>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// ============================================================================
//  Data Structures & Buffers (For Speed)
// ============================================================================

struct SubstMatrix {
  std::vector<std::vector<int>> scores;
  std::vector<int> char_map; // Faster than unordered_map for lookups
  int gap_open;
  int gap_extend;
  
  SubstMatrix() : char_map(256, -1) {} // Map ASCII chars to indices
  
  inline int get_score(unsigned char a, unsigned char b) const {
    int i = char_map[a];
    int j = char_map[b];
    if (i < 0 || j < 0) return -4; 
    return scores[i][j];
  }
};

// Thread-local workspace to prevent re-allocation in loops
struct DPWorkspace {
  std::vector<int> v1;
  std::vector<std::vector<int>> mat;
  std::vector<int> last_row_tracker; // For Damerau
  
  void resize_mat(int n, int m) {
    if ((int)mat.size() <= n) mat.resize(n + 5);
    for(auto &r : mat) {
      if ((int)r.size() <= m) r.resize(m + 5);
    }
  }
};

// ============================================================================
//  Distance Implementations (Buffer-aware)
// ============================================================================

// 1. Levenshtein (Banded)
static inline int levenshtein_opt(const std::string &a, const std::string &b, int thr, std::vector<int>& row) {
  int n = a.size(), m = b.size();
  if (std::abs(n - m) > thr) return thr + 1;
  if ((int)row.size() <= m) row.resize(m + 1);
  
  std::iota(row.begin(), row.begin() + m + 1, 0);
  
  for (int i = 1; i <= n; ++i) {
    int prev = row[0];
    row[0] = i;
    int row_min = i;
    
    int start_j = std::max(1, i - thr);
    int end_j   = std::min(m, i + thr);
    
    // Fill valid band
    if (start_j > 1) row[start_j - 1] = thr + 2; // Guard boundary
    
    for (int j = start_j; j <= end_j; ++j) {
      int cur = row[j];
      int cost = (a[i-1] == b[j-1] ? 0 : 1);
      row[j] = std::min({ row[j-1] + 1, cur + 1, prev + cost });
      prev = cur;
      row_min = std::min(row_min, row[j]);
    }
    if (row_min > thr) return thr + 1;
  }
  return row[m];
}

// 3. Damerau-Levenshtein (Buffered)
static inline int damerau_opt(const std::string &a, const std::string &b, int thr, DPWorkspace &ws) {
  int n = a.size(), m = b.size();
  if (std::abs(n - m) > thr) return thr + 1;
  
  ws.resize_mat(n + 2, m + 2);
  auto &H = ws.mat;
  
  // Reset borders
  const int INF = thr + 2; // Optimization: Don't need full n+m
  H[0][0] = INF;
  for (int i = 0; i <= n + 1; ++i) { H[i][0] = INF; H[i][1] = i > 0 ? i-1 : 0; }
  for (int j = 0; j <= m + 1; ++j) { H[0][j] = INF; H[1][j] = j > 0 ? j-1 : 0; }
  
  // Char tracker (ASCII 0-255)
  if (ws.last_row_tracker.empty()) ws.last_row_tracker.resize(256, 0);
  std::fill(ws.last_row_tracker.begin(), ws.last_row_tracker.end(), 0);
  
  for (int i = 1; i <= n; ++i) {
    int last_match_col = 0;
    int row_min = INF;
    
    for (int j = 1; j <= m; ++j) {
      int i1 = i + 1; // 1-based index in H for string char i
      int j1 = j + 1;
      
      int last_match_row = ws.last_row_tracker[(unsigned char)b[j-1]];
      int cost = (a[i-1] == b[j-1] ? 0 : 1);
      
      int d_sub = H[i1-1][j1-1] + cost;
      int d_ins = H[i1][j1-1] + 1;
      int d_del = H[i1-1][j1] + 1;
      int d_trans = INF;
      
      if (last_match_row > 0 && last_match_col > 0) {
        d_trans = H[last_match_row][last_match_col] + (i - last_match_row - 1) + 1 + (j - last_match_col - 1);
      }
      
      H[i1][j1] = std::min({d_sub, d_ins, d_del, d_trans});
      if (cost == 0) last_match_col = j;
      row_min = std::min(row_min, H[i1][j1]);
    }
    ws.last_row_tracker[(unsigned char)a[i-1]] = i;
    if (row_min > thr) return thr + 1; 
  }
  return H[n+1][m+1];
}

// 4. Needleman-Wunsch (Buffered)
static inline int nw_opt(const std::string &a, const std::string &b, const SubstMatrix &mat, DPWorkspace &ws) {
  int n = a.size(), m = b.size();
  ws.resize_mat(n + 1, m + 1);
  auto &dp = ws.mat;
  
  // Init with proper gap penalties
  dp[0][0] = 0;
  for (int i = 1; i <= n; ++i) dp[i][0] = mat.gap_open + (i-1) * mat.gap_extend;
  for (int j = 1; j <= m; ++j) dp[0][j] = mat.gap_open + (j-1) * mat.gap_extend;
  
  for (int i = 1; i <= n; ++i) {
    for (int j = 1; j <= m; ++j) {
      int match = dp[i-1][j-1] + mat.get_score(a[i-1], b[j-1]);
      int del = dp[i-1][j] + mat.gap_extend; // Linear approximation for speed (vs affine state machine)
      int ins = dp[i][j-1] + mat.gap_extend;
      dp[i][j] = std::max({match, del, ins});
    }
  }
  // Convert to distance-like (inverted score)
  int max_s = std::max(n, m) * 5; 
  return std::max(0, max_s - dp[n][m]);
}

// 5. Smith-Waterman (Buffered)
static inline int sw_opt(const std::string &a, const std::string &b, const SubstMatrix &mat, DPWorkspace &ws) {
  int n = a.size(), m = b.size();
  ws.resize_mat(n + 1, m + 1);
  auto &dp = ws.mat;
  
  int max_val = 0;
  // Init (zeros for local alignment)
  for (int i = 0; i <= n; ++i) dp[i][0] = 0;
  for (int j = 0; j <= m; ++j) dp[0][j] = 0;
  
  for (int i = 1; i <= n; ++i) {
    for (int j = 1; j <= m; ++j) {
      int match = dp[i-1][j-1] + mat.get_score(a[i-1], b[j-1]);
      int del = dp[i-1][j] + mat.gap_extend;
      int ins = dp[i][j-1] + mat.gap_extend;
      dp[i][j] = std::max({0, match, del, ins});
      if (dp[i][j] > max_val) max_val = dp[i][j];
    }
  }
  int max_possible = std::min(n, m) * 5; 
  return std::max(0, max_possible - max_val);
}

// ============================================================================
//  Matrix Parser Helper
// ============================================================================
SubstMatrix parse_matrix(NumericMatrix mat, CharacterVector rnames, int go, int ge) {
  SubstMatrix sm;
  sm.gap_open = go;
  sm.gap_extend = ge;
  int n = mat.nrow();
  sm.scores.resize(n, std::vector<int>(n));
  
  for(int i=0; i<n; ++i) {
    std::string s = as<std::string>(rnames[i]);
    if(!s.empty()) sm.char_map[(unsigned char)s[0]] = i;
    for(int j=0; j<n; ++j) sm.scores[i][j] = (int)mat(i,j);
  }
  return sm;
}

// ============================================================================
//  Main Export
// ============================================================================
// [[Rcpp::export]]
DataFrame fast_edge_list(CharacterVector seqs,
                         double thresh = 1.0,
                         Nullable<CharacterVector> v_gene = R_NilValue,
                         Nullable<CharacterVector> j_gene = R_NilValue,
                         bool match_v = false,
                         bool match_j = false,
                         Nullable<CharacterVector> ids = R_NilValue,
                         std::string metric = "levenshtein",
                         std::string normalize = "none",
                         Nullable<NumericMatrix> subst_matrix = R_NilValue,
                         int gap_open = -10,
                         int gap_extend = -1) 
{
  int n = seqs.size();
  if (n < 2) stop("Need at least 2 sequences.");
  
  // 1. Prepare Data
  std::vector<std::string> s(n);
  std::vector<int> lens(n);
  for(int i=0; i<n; ++i) { s[i] = as<std::string>(seqs[i]); lens[i] = s[i].size(); }
  
  std::vector<std::string> v_vec(n), j_vec(n);
  if(match_v && v_gene.isNotNull()) v_vec = as<std::vector<std::string>>(v_gene);
  if(match_j && j_gene.isNotNull()) j_vec = as<std::vector<std::string>>(j_gene);
  
  std::vector<std::string> lbl(n);
  if(ids.isNotNull()) lbl = as<std::vector<std::string>>(ids);
  else for(int i=0; i<n; ++i) lbl[i] = std::to_string(i+1);
  
  // 2. Prepare Grouping
  std::map<std::pair<std::string,std::string>, std::vector<int>> group_map;
  for(int i=0; i<n; ++i) {
    group_map[{match_v ? v_vec[i] : "", match_j ? j_vec[i] : ""}].push_back(i);
  }
  
  // 3. Prepare Substitution Matrix
  SubstMatrix smat;
  if (metric == "nw" || metric == "sw") {
    if (subst_matrix.isNotNull()) {
      NumericMatrix nm(subst_matrix);
      // FIXED: Ambiguous operator[] access
      List dimnames = nm.attr("dimnames");
      smat = parse_matrix(nm, dimnames[0], gap_open, gap_extend);
    } else {
      // Fallback Identity
      CharacterVector abc = CharacterVector::create("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y");
      NumericMatrix ident(20,20);
      
      // FIXED: Incorrect loop syntax/logic for matrix fill
      std::fill(ident.begin(), ident.end(), -1); // Fill all with mismatch (-1)
      for(int i=0; i<20; ++i) ident(i,i) = 1;    // Set diagonal to match (1)
      
      smat = parse_matrix(ident, abc, gap_open, gap_extend);
    }
  }
  
  // 4. Output Storage
  std::vector<std::string> out_from, out_to;
  std::vector<double> out_d;
  
  // 5. Parallel Processing
#ifdef _OPENMP
  int nthreads = omp_get_max_threads();
#else
  int nthreads = 1;
#endif
  
  std::vector<std::vector<int>> all_groups;
  for(auto &kv : group_map) if(kv.second.size() > 1) all_groups.push_back(kv.second);
  
  for (const auto& grp : all_groups) {
    int g_size = grp.size();
    
#pragma omp parallel num_threads(nthreads)
{
  // Thread-Local Storage
  DPWorkspace ws; 
  std::vector<std::string> loc_f, loc_t;
  std::vector<double> loc_d_vec;
  
#pragma omp for schedule(dynamic, 64)
  for (int i_idx = 0; i_idx < g_size; ++i_idx) {
    int i = grp[i_idx];
    for (int k_idx = i_idx + 1; k_idx < g_size; ++k_idx) {
      int k = grp[k_idx];
      
      int max_l = std::max(lens[i], lens[k]);
      if (max_l == 0) continue;
      
      int maxd;
      if (thresh >= 1.0) {
        // Absolute distance threshold
        maxd = (int)thresh;
      } else {
        // Normalized distance threshold - use same length as normalization
        double norm_len;
        if (normalize == "maxlen") {
          norm_len = (double)max_l;
        } else if (normalize == "length") {
          norm_len = (lens[i] + lens[k]) / 2.0;
        } else {
          // normalize == "none": interpret threshold as fraction of max_l
          norm_len = (double)max_l;
        }
        maxd = (int)(thresh * norm_len);
      }
      // =====================================================================
      
      if (metric != "nw" && metric != "sw" && std::abs(lens[i] - lens[k]) > maxd) continue;
      
      int d = 0;
      if (metric == "levenshtein") d = levenshtein_opt(s[i], s[k], maxd, ws.v1);
      else if (metric == "hamming") {
        if (lens[i] != lens[k]) d = maxd + 1;
        else {
          for(int p=0; p<lens[i]; ++p) if(s[i][p] != s[k][p]) { d++; if(d>maxd) break; }
        }
      }
      else if (metric == "damerau") d = damerau_opt(s[i], s[k], maxd, ws);
      else if (metric == "nw")      d = nw_opt(s[i], s[k], smat, ws);
      else if (metric == "sw")      d = sw_opt(s[i], s[k], smat, ws);
      
      if (d <= maxd) {
        double fd = (double)d;
        if (normalize == "maxlen") fd /= max_l;
        else if (normalize == "length") fd /= ((lens[i]+lens[k])/2.0);
        
        if (thresh < 1.0 && normalize != "none" && fd > thresh) {
          continue;  // Skip edges that exceed normalized threshold
        }
        // =====================================================================
        
        loc_f.push_back(lbl[i]);
        loc_t.push_back(lbl[k]);
        loc_d_vec.push_back(fd);
      }
    }
  }
#pragma omp critical
{
  out_from.insert(out_from.end(), loc_f.begin(), loc_f.end());
  out_to.insert(out_to.end(), loc_t.begin(), loc_t.end());
  out_d.insert(out_d.end(), loc_d_vec.begin(), loc_d_vec.end());
}
}
  }
  
  return DataFrame::create(_["from"] = out_from, _["to"] = out_to, _["dist"] = out_d, _["stringsAsFactors"] = false);
}
