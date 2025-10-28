// ── src/fastEditEdges.cpp ───────────────────────────────────────────────────
#include <Rcpp.h>
#include <numeric>
#include <cmath>          
#include <string>
#include <vector>
#include <algorithm>
#include <unordered_map>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// ============================================================================
//  Substitution Matrix Support
// ============================================================================

struct SubstMatrix {
  std::unordered_map<char, int> char_index;
  std::vector<std::vector<int>> scores;
  int gap_open;
  int gap_extend;
  
  int get_score(char a, char b) const {
    auto it_a = char_index.find(a);
    auto it_b = char_index.find(b);
    if (it_a == char_index.end() || it_b == char_index.end()) {
      return -4; // Default penalty for unknown characters
    }
    return scores[it_a->second][it_b->second];
  }
};

// Parse R matrix into SubstMatrix structure
SubstMatrix parse_r_matrix(NumericMatrix mat, 
                           CharacterVector rownames,
                           int gap_open,
                           int gap_extend) {
  SubstMatrix smat;
  smat.gap_open = gap_open;
  smat.gap_extend = gap_extend;
  
  int n = mat.nrow();
  
  // Build character index from row names
  for (int i = 0; i < n; ++i) {
    std::string name = as<std::string>(rownames[i]);
    if (name.length() > 0) {
      smat.char_index[name[0]] = i;
    }
  }
  
  // Copy matrix values
  smat.scores.resize(n, std::vector<int>(n));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      smat.scores[i][j] = static_cast<int>(mat(i, j));
    }
  }
  
  return smat;
}

// Fallback identity matrix if no matrix provided
SubstMatrix get_identity_matrix() {
  SubstMatrix mat;
  mat.gap_open = -2;
  mat.gap_extend = -1;
  
  std::string alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  for (size_t i = 0; i < alphabet.size(); ++i) {
    mat.char_index[alphabet[i]] = i;
  }
  
  mat.scores.resize(26, std::vector<int>(26));
  for (int i = 0; i < 26; ++i) {
    for (int j = 0; j < 26; ++j) {
      mat.scores[i][j] = (i == j) ? 1 : -1;
    }
  }
  
  return mat;
}

// ============================================================================
//  Distance Metrics
// ============================================================================

// ----------------------------------------------------------------------------
//  1. Levenshtein with early-exit (banded DP)
// ----------------------------------------------------------------------------
static inline int levenshtein_distance(const std::string &a,
                                       const std::string &b,
                                       int thr)
{
  int n = a.size(), m = b.size();
  if (std::abs(n - m) > thr) return thr + 1;
  
  std::vector<int> row(m + 1);
  std::iota(row.begin(), row.end(), 0);
  
  for (int i = 1; i <= n; ++i) {
    int prev = row[0];
    row[0] = i;
    int row_min = i;
    
    for (int j = 1; j <= m; ++j) {
      int cur = row[j];
      int cost = (a[i-1] == b[j-1] ? 0 : 1);
      row[j] = std::min({ row[j-1] + 1,
                        row[j] + 1,
                        prev + cost });
      prev = cur;
      row_min = std::min(row_min, row[j]);
    }
    if (row_min > thr) return thr + 1;
  }
  return row[m];
}

// ----------------------------------------------------------------------------
//  2. Hamming distance (requires equal length)
// ----------------------------------------------------------------------------
static inline int hamming_distance(const std::string &a,
                                   const std::string &b,
                                   int thr)
{
  if (a.size() != b.size()) return thr + 1;
  
  int d = 0;
  for (size_t i = 0; i < a.size(); ++i) {
    if (a[i] != b[i]) {
      d++;
      if (d > thr) return thr + 1; // Early exit
    }
  }
  return d;
}

// ----------------------------------------------------------------------------
//  3. Damerau-Levenshtein distance (with transpositions)
// ----------------------------------------------------------------------------
static inline int damerau_distance(const std::string &a,
                                   const std::string &b,
                                   int thr)
{
  int n = a.size(), m = b.size();
  if (std::abs(n - m) > thr) return thr + 1;
  
  std::vector<std::vector<int>> H(n + 2, std::vector<int>(m + 2));
  
  const int INF = n + m;
  H[0][0] = INF;
  for (int i = 0; i <= n; ++i) { H[i+1][0] = INF; H[i+1][1] = i; }
  for (int j = 0; j <= m; ++j) { H[0][j+1] = INF; H[1][j+1] = j; }
  
  std::unordered_map<char, int> last_row;
  
  for (int i = 1; i <= n; ++i) {
    int last_match_col = 0;
    for (int j = 1; j <= m; ++j) {
      int last_match_row = last_row[b[j-1]];
      int cost = (a[i-1] == b[j-1] ? 0 : 1);
      
      H[i+1][j+1] = std::min({
        H[i][j] + cost,                                    // substitution
        H[i+1][j] + 1,                                     // insertion
        H[i][j+1] + 1,                                     // deletion
        H[last_match_row][last_match_col] + (i - last_match_row - 1) + 1 + (j - last_match_col - 1) // transposition
      });
      
      if (cost == 0) last_match_col = j;
    }
    last_row[a[i-1]] = i;
    
    // Early exit check
    int row_min = H[i+1][1];
    for (int j = 2; j <= m + 1; ++j) {
      row_min = std::min(row_min, H[i+1][j]);
    }
    if (row_min > thr) return thr + 1;
  }
  
  return H[n+1][m+1];
}

// ----------------------------------------------------------------------------
//  4. Needleman-Wunsch (global alignment)
// ----------------------------------------------------------------------------
static inline int needleman_wunsch(const std::string &a,
                                   const std::string &b,
                                   const SubstMatrix &mat)
{
  int n = a.size(), m = b.size();
  std::vector<std::vector<int>> dp(n + 1, std::vector<int>(m + 1));
  
  // Initialize
  dp[0][0] = 0;
  for (int i = 1; i <= n; ++i) dp[i][0] = dp[i-1][0] + mat.gap_extend;
  for (int j = 1; j <= m; ++j) dp[0][j] = dp[0][j-1] + mat.gap_extend;
  
  // Fill
  for (int i = 1; i <= n; ++i) {
    for (int j = 1; j <= m; ++j) {
      int match = dp[i-1][j-1] + mat.get_score(a[i-1], b[j-1]);
      int del = dp[i-1][j] + mat.gap_extend;
      int ins = dp[i][j-1] + mat.gap_extend;
      dp[i][j] = std::max({match, del, ins});
    }
  }
  
  // Convert score to distance (higher score = lower distance)
  // Normalize to 0..max_len range
  int max_score = std::max(n, m) * 5; // Approximate max possible score
  int distance = max_score - dp[n][m];
  return std::max(0, distance);
}

// ----------------------------------------------------------------------------
//  5. Smith-Waterman (local alignment)
// ----------------------------------------------------------------------------
static inline int smith_waterman(const std::string &a,
                                 const std::string &b,
                                 const SubstMatrix &mat)
{
  int n = a.size(), m = b.size();
  std::vector<std::vector<int>> dp(n + 1, std::vector<int>(m + 1, 0));
  
  int max_score = 0;
  
  // Fill
  for (int i = 1; i <= n; ++i) {
    for (int j = 1; j <= m; ++j) {
      int match = dp[i-1][j-1] + mat.get_score(a[i-1], b[j-1]);
      int del = dp[i-1][j] + mat.gap_extend;
      int ins = dp[i][j-1] + mat.gap_extend;
      dp[i][j] = std::max({0, match, del, ins});
      max_score = std::max(max_score, dp[i][j]);
    }
  }
  
  // Convert score to distance
  int max_possible = std::min(n, m) * 5;
  int distance = max_possible - max_score;
  return std::max(0, distance);
}

// ============================================================================
//  Memory-Efficient Grouping Helper
// ============================================================================

struct GroupInfo {
  std::vector<int> indices;
  std::string v_gene;
  std::string j_gene;
};

std::vector<GroupInfo> create_groups(
    int n,
    const std::vector<std::string>& v,
    const std::vector<std::string>& J,
    bool match_v,
    bool match_j)
{
  std::vector<GroupInfo> groups;
  
  if (!match_v && !match_j) {
    // No grouping - all sequences in one group
    GroupInfo g;
    g.indices.resize(n);
    std::iota(g.indices.begin(), g.indices.end(), 0);
    groups.push_back(g);
    return groups;
  }
  
  // Create groups based on V/J genes
  std::map<std::pair<std::string, std::string>, GroupInfo> group_map;
  
  for (int i = 0; i < n; ++i) {
    std::string v_key = match_v ? v[i] : "";
    std::string j_key = match_j ? J[i] : "";
    
    auto key = std::make_pair(v_key, j_key);
    group_map[key].indices.push_back(i);
    group_map[key].v_gene = v_key;
    group_map[key].j_gene = j_key;
  }
  
  for (auto& kv : group_map) {
    groups.push_back(kv.second);
  }
  
  return groups;
}

// ============================================================================
//  Main export 
// ============================================================================
// [[Rcpp::export]]
DataFrame fast_edge_list(CharacterVector           seqs,
                         double                    thresh        = 1.0,
                         Nullable<CharacterVector> v_gene        = R_NilValue,
                         Nullable<CharacterVector> j_gene        = R_NilValue,
                         bool                      match_v       = false,
                         bool                      match_j       = false,
                         Nullable<CharacterVector> ids           = R_NilValue,
                         std::string               metric        = "levenshtein",
                         std::string               normalize     = "none",
                         Nullable<NumericMatrix>   subst_matrix  = R_NilValue,
                         int                       gap_open      = -10,
                         int                       gap_extend    = -1)
{
  int n = seqs.size();
  if (n < 2) stop("At least two sequences are required.");
  if (thresh < 0) stop("`thresh` must be >= 0.");
  
  // Validate metric
  std::vector<std::string> valid_metrics = {"levenshtein", "hamming", "damerau", "nw", "sw"};
  if (std::find(valid_metrics.begin(), valid_metrics.end(), metric) == valid_metrics.end()) {
    stop("Invalid metric. Must be one of: levenshtein, hamming, damerau, nw, sw");
  }
  
  // Validate normalize
  std::vector<std::string> valid_norms = {"none", "length", "maxlen"};
  if (std::find(valid_norms.begin(), valid_norms.end(), normalize) == valid_norms.end()) {
    stop("Invalid normalize. Must be one of: none, length, maxlen");
  }
  
  // 1) Read sequences + lengths
  std::vector<std::string> s(n);
  std::vector<int> lens(n);
  for (int i = 0; i < n; ++i) {
    s[i] = as<std::string>(seqs[i]);
    lens[i] = s[i].size();
  }
  
  // 2) Optional V/J
  std::vector<std::string> v(n), J(n);
  if (match_v) {
    if (v_gene.isNull()) stop("match_v=TRUE but no v_gene provided.");
    CharacterVector vv(v_gene);
    if (vv.size() == 1) vv = CharacterVector(n, vv[0]);
    if ((int)vv.size() != n) stop("Length of v_gene != number of sequences.");
    v = as<std::vector<std::string>>(vv);
  }
  if (match_j) {
    if (j_gene.isNull()) stop("match_j=TRUE but no j_gene provided.");
    CharacterVector jj(j_gene);
    if (jj.size() == 1) jj = CharacterVector(n, jj[0]);
    if ((int)jj.size() != n) stop("Length of j_gene != number of sequences.");
    J = as<std::vector<std::string>>(jj);
  }
  
  // 3) IDs / labels
  std::vector<std::string> lbl(n);
  if (ids.isNull()) {
    for (int i = 0; i < n; ++i) lbl[i] = std::to_string(i + 1);
  } else {
    CharacterVector ix(ids);
    if ((int)ix.size() != n) stop("Length of ids != number of sequences.");
    lbl = as<std::vector<std::string>>(ix);
  }
  
  // 4) Get substitution matrix (if needed for alignment-based metrics)
  SubstMatrix smat;
  if (metric == "nw" || metric == "sw") {
    if (subst_matrix.isNull()) {
      // Use identity matrix as fallback
      smat = get_identity_matrix();
      smat.gap_open = gap_open;
      smat.gap_extend = gap_extend;
    } else {
      NumericMatrix mat(subst_matrix);
      if (!mat.hasAttribute("dimnames")) {
        stop("Substitution matrix must have row/column names (amino acid labels)");
      }
      List dimnames = mat.attr("dimnames");
      CharacterVector rownames = dimnames[0];
      
      smat = parse_r_matrix(mat, rownames, gap_open, gap_extend);
    }
  }
  
  // 5) Create groups for memory efficiency
  auto groups = create_groups(n, v, J, match_v, match_j);
  
  // 6) Estimate output size based on groups
  size_t total_pairs = 0;
  for (const auto& group : groups) {
    size_t g_size = group.indices.size();
    total_pairs += g_size * (g_size - 1) / 2;
  }
  
  std::vector<std::string> out_from;
  std::vector<std::string> out_to;
  std::vector<double> out_d;
  out_from.reserve(total_pairs / 2); // Conservative estimate
  out_to.reserve(total_pairs / 2);
  out_d.reserve(total_pairs / 2);
  
#ifdef _OPENMP
  int nthread = omp_get_max_threads();
#else
  int nthread = 1;
#endif
  
  // Process each group
  for (const auto& group : groups) {
    int g_size = group.indices.size();
    if (g_size < 2) continue;
    
#ifdef _OPENMP
#pragma omp parallel
#endif
{
  std::vector<std::string> loc_from;
  std::vector<std::string> loc_to;
  std::vector<double> loc_d;
  loc_from.reserve(g_size * g_size / (2 * nthread));
  loc_to.reserve(g_size * g_size / (2 * nthread));
  loc_d.reserve(g_size * g_size / (2 * nthread));
  
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 32)
#endif
  for (int gi = 0; gi < g_size; ++gi) {
    int i = group.indices[gi];
    
    for (int gk = gi + 1; gk < g_size; ++gk) {
      int k = group.indices[gk];
      
      int max_len = std::max(lens[i], lens[k]);
      if (max_len == 0) continue;
      
      // Determine threshold
      int maxd;
      if (thresh >= 1.0) {
        maxd = std::lround(thresh);
      } else {
        double dist_thresh_rel = 1.0 - thresh;
        maxd = std::floor(dist_thresh_rel * max_len);
      }
      
      // Skip if length difference too large (for edit distances)
      if (metric != "nw" && metric != "sw") {
        if (std::abs(lens[i] - lens[k]) > maxd) continue;
      }
      
      // Compute distance
      int d = 0;
      if (metric == "levenshtein") {
        d = levenshtein_distance(s[i], s[k], maxd);
      } else if (metric == "hamming") {
        d = hamming_distance(s[i], s[k], maxd);
      } else if (metric == "damerau") {
        d = damerau_distance(s[i], s[k], maxd);
      } else if (metric == "nw") {
        d = needleman_wunsch(s[i], s[k], smat);
      } else if (metric == "sw") {
        d = smith_waterman(s[i], s[k], smat);
      }
      
      // Filter by threshold
      if (d > maxd) continue;
      
      // Normalize if requested
      double final_dist = static_cast<double>(d);
      if (normalize == "maxlen") {
        final_dist = final_dist / max_len;
      } else if (normalize == "length") {
        // Use mean length
        double mean_len = (lens[i] + lens[k]) / 2.0;
        final_dist = final_dist / mean_len;
      }
      // else normalize == "none", keep raw distance
      
      // Store edge
      loc_from.push_back(lbl[i]);
      loc_to.push_back(lbl[k]);
      loc_d.push_back(final_dist);
    }
  }
  
#ifdef _OPENMP
#pragma omp critical
#endif
{
  out_from.insert(out_from.end(), loc_from.begin(), loc_from.end());
  out_to.insert(out_to.end(), loc_to.begin(), loc_to.end());
  out_d.insert(out_d.end(), loc_d.begin(), loc_d.end());
}
}
  }
  
  return DataFrame::create(
    _["from"] = out_from,
    _["to"] = out_to,
    _["dist"] = out_d,
    _["stringsAsFactors"] = false
  );
}