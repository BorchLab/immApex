// -*- mode: C++; c-indent-level: 2; indent-tabs-mode: nil; -*-
// encodeSequences.cpp
//
// [[Rcpp::depends(Rcpp)]]
// We ask only for C++11 – newer standards are accepted automatically.
// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
#include <unordered_map>
#include <string>
#include <vector>
#include <algorithm>   // std::max, std::transform

// --- OpenMP is used **only** if the compiler defines _OPENMP ---------------
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// -----------------------------------------------------------------------------
// 3-D index helper:  dim = c(depth, length, samples)   (column-major - default)
// -----------------------------------------------------------------------------
inline std::size_t idx3d(const int s,
                         const int pos,
                         const int d,
                         const int S,
                         const int L,
                         const int D) {
  // offset = d + D*pos + D*L*s
  return static_cast<std::size_t>(d) +
    static_cast<std::size_t>(D) * (pos + L * s);
}

// -----------------------------------------------------------------------------
//  Main encoder
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List encodeSequences_cpp(const CharacterVector &sequences,
                               std::string           mode              = "onehot",
                               const CharacterVector alphabet          = CharacterVector::create("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"),
                               Nullable<Rcpp::NumericMatrix> prop_mat_ = R_NilValue,
                               const char             pad_token        = '.',
                               std::string            summary          = "",
                               int                    max_len          = -1,
                               int                    nthreads         = 1)
{
  // -------------------------------------------------------------------------
  // 0.  Validate input & build lookup table
  // -------------------------------------------------------------------------
  const int S = sequences.size();
  if (S == 0) Rcpp::stop("`sequences` is empty.");
  
  std::transform(mode.begin(), mode.end(), mode.begin(), ::tolower);
  if (mode != "onehot" && mode != "property")
    Rcpp::stop("`mode` must be either \"onehot\" or \"property\".");
  
  const int K_in = alphabet.size();
  if (K_in == 0) Rcpp::stop("`alphabet` cannot be empty.");
  
  std::unordered_map<char,int> aa2idx;
  aa2idx.reserve(static_cast<std::size_t>(K_in) + 1);
  
  for (int i = 0; i < K_in; ++i) {
    const std::string aa = Rcpp::as<std::string>(alphabet[i]);
    if (aa.size() != 1)
      Rcpp::stop("`alphabet` must contain single-character strings; "
                   "element %d is \"%s\".", i + 1, aa.c_str());
    aa2idx[ aa[0] ] = i;
  }
  
  // ensure padding symbol exists
  if (aa2idx.find(pad_token) == aa2idx.end())
    aa2idx[ pad_token ] = static_cast<int>(aa2idx.size());
  
  const int K = static_cast<int>(aa2idx.size());     // alphabet (+pad) size
  
  // -------------------------------------------------------------------------
  // 1.  Property matrix (only for property mode)
  // -------------------------------------------------------------------------
  Rcpp::NumericMatrix prop_mat;    // empty unless mode=="property"
  int P = 0;
  
  if (mode == "property") {
    if (prop_mat_.isNull())
      Rcpp::stop("`prop_mat` is required when mode == \"property\".");
    
    prop_mat = prop_mat_.get();
    
    if (prop_mat.nrow() != K_in && prop_mat.nrow() != K)
      Rcpp::stop("`prop_mat` rows (%d) must match `alphabet` length (%d).",
                 prop_mat.nrow(), K_in);
    
    P = prop_mat.ncol();
    if (P == 0) Rcpp::stop("`prop_mat` has zero columns.");
  }
  
  // -------------------------------------------------------------------------
  // 2.  Determine max_len if not supplied
  // -------------------------------------------------------------------------
  if (max_len <= 0) {
    max_len = 0;
    for (int i = 0; i < S; ++i) {
      const std::string seq = Rcpp::as<std::string>(sequences[i]);
      max_len = std::max<int>(max_len, static_cast<int>(seq.size()));
    }
  }
  
  // -------------------------------------------------------------------------
  // 3.  Thread setup – only effective if compiled with OpenMP
  // -------------------------------------------------------------------------
#ifdef _OPENMP
  nthreads = std::max(1, nthreads);
  omp_set_num_threads(nthreads);
#else
  nthreads = 1;            // guarantee single-thread fall-back
#endif
  
  // -------------------------------------------------------------------------
  // 4.  Allocate outputs
  // -------------------------------------------------------------------------
  const int D = (mode == "onehot") ? K : P;          // depth dimension
  
  Rcpp::NumericVector cube( static_cast<R_xlen_t>(S) * max_len * D );
  cube.attr("dim") = Rcpp::IntegerVector::create(D, max_len, S);
  
  Rcpp::NumericMatrix flat(S, max_len * D);          // row-major 2-D view
  
  const bool need_summary = (mode == "property" && summary == "mean");
  Rcpp::NumericMatrix summaryMat;
  if (need_summary) summaryMat = Rcpp::NumericMatrix(S, P);
  
  // -------------------------------------------------------------------------
  // 5.  Main (possibly parallel) loop
  // -------------------------------------------------------------------------
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int s = 0; s < S; ++s) {
    
    const std::string seq = Rcpp::as<std::string>(sequences[s]);
    const int L_obs       = static_cast<int>(seq.size());
    
    std::vector<double> acc( (mode == "property") ? P : 0, 0.0 );
    
    // Get the index for the padding token once
    const int pad_idx = aa2idx.at(pad_token);
    
    for (int pos = 0; pos < max_len; ++pos) {
      
      const char aa_char = (pos < L_obs) ? seq[pos] : pad_token;
      
      // Find the index, defaulting to the padding index if the character is not in the alphabet
      const int  aa_idx  = (aa2idx.count(aa_char)) ? aa2idx.at(aa_char) : pad_idx;
      
      if (mode == "onehot") {
        // This logic is correct as `D` (depth) is K, which includes the padding symbol
        const std::size_t off = idx3d(s, pos, aa_idx, S, max_len, D);
        cube[ off ] = 1.0;
        flat(s, pos * D + aa_idx) = 1.0;
        
      } else { // property mode
        
        // **FIX**: Check if the character is a pad/unknown symbol BEFORE accessing prop_mat
        if (aa_idx == pad_idx) {
          // For padded/unknown positions, encode with zeros
          for (int p = 0; p < P; ++p) {
            const std::size_t off = idx3d(s, pos, p, S, max_len, D);
            cube[ off ] = 0.0;
            flat(s, pos * D + p) = 0.0;
          }
        } else {
          // For valid amino acids, access the property matrix (this is now safe)
          for (int p = 0; p < P; ++p) {
            const double val = prop_mat(aa_idx, p);
            const std::size_t off = idx3d(s, pos, p, S, max_len, D);
            cube[ off ] = val;
            flat(s, pos * D + p) = val;
            if (need_summary) acc[p] += val;
          }
        }
      }
    } // end position loop
    
    if (need_summary && L_obs > 0) {
      for (int p = 0; p < P; ++p)
        summaryMat(s, p) = acc[p] / static_cast<double>(L_obs);
    }
  } // end sample loop
  
  // -------------------------------------------------------------------------
  // 6.  Assemble return list
  // -------------------------------------------------------------------------
  Rcpp::List out = Rcpp::List::create(
    Rcpp::_["cube"]      = cube,
    Rcpp::_["flattened"] = flat
  );
  if (need_summary) out["summary"] = summaryMat;
  
  out["sequence.dictionary"]  = alphabet;
  out["padding.symbol"] = std::string(1, pad_token);
#ifdef _OPENMP
  out["threads"]   = nthreads;
#else
  out["threads"]   = 1;
#endif
  
  return out;
}
