// -*- mode: C++; c-indent-level: 2; indent-tabs-mode: nil; -*-
// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(cpp17)]]
// Enable OpenMP where available; CRAN/Windows fall back to serial
// [[Rcpp::plugins(openmp)]]

#include <Rcpp.h>
#include <unordered_map>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// -----------------------------------------------------------------------------
// Internal helper: thread-safe increment of unordered_map
// -----------------------------------------------------------------------------
inline void add_motif(std::unordered_map<std::string, std::size_t> &tbl,
                      const std::string &motif) {
#ifdef _OPENMP
#pragma omp critical
#endif
{
  ++tbl[motif];
}
}

// -----------------------------------------------------------------------------
// Exported engine
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List calculateMotif_cpp(const CharacterVector &sequences,
                        const IntegerVector  &motif_lengths,
                        const bool  discontinuous      = false,
                        const char  gap_char           = '.',
                        const int   nthreads_requested = 1) {
  
  std::unordered_map<std::string, std::size_t> counts;
  
  // Decide on number of threads
  int nthreads = std::max(1, nthreads_requested);
#ifdef _OPENMP
  omp_set_num_threads(nthreads);
#else
  nthreads = 1;  // ensure consistent return
#endif
  
  // ---------------------------------------------------------------------------
  // Main parallel loop
  // ---------------------------------------------------------------------------
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < sequences.size(); ++i) {
    const std::string seq = Rcpp::as<std::string>(sequences[i]);
    const std::size_t L   = seq.size();
    
    for (int k : motif_lengths) {
      if (k <= 0 || static_cast<std::size_t>(k) > L) continue;
      
      for (std::size_t pos = 0; pos + k <= L; ++pos) {
        const std::string motif = seq.substr(pos, k);
        add_motif(counts, motif);
        
        if (discontinuous) {
          // generate all single-gap variants at this window
          for (int g = 0; g < k; ++g) {
            std::string gmotif = motif;
            gmotif[g] = gap_char;
            add_motif(counts, gmotif);
          }
        }
      }
    }
  }
  
  // ---------------------------------------------------------------------------
  // Transfer results back to R
  // ---------------------------------------------------------------------------
  const std::size_t N = counts.size();
  CharacterVector motifs(N);
  IntegerVector   freqs(N);
  
  std::size_t idx = 0;
  for (const auto &kv : counts) {
    motifs[idx] = kv.first;
    freqs[idx]  = static_cast<int>(kv.second);
    ++idx;
  }
  
  return List::create(
    _["motif"]     = motifs,
    _["frequency"] = freqs      
  );
}