#' Motif Enumeration and Counting 
#'
#' @description
#' Rapidly enumerates and quantifies **contiguous** (and, optionally, 
#' single-gap discontinuous) amino-acid motifs across a set of sequences. 
#'
#' @details
#' For every input sequence the algorithm slides windows of length *k*
#' (`motif.lengths`) and increments a motif counter (`unordered_map`).
#' If `discontinuous = TRUE`, each window is additionally copied *k*
#' times, substituting one position at a time with `discontinuous.symbol`
#' (default `"."`), yielding gapped motif patterns such as `"C.S"`.
#'
#' @param input.sequences Character vector of sequences (amino acid or 
#' nucleotide)
#' @param motif.lengths Integer vector of motif sizes (>= 1). **Default:** `2:5`.
#' @param min.depth Minimum count a motif must reach to be retained in the 
#' output (`>= 1`). **Default:** `3`.
#' @param discontinuous Logical; include single-gap motifs as well? 
#' **Default:** `FALSE`.
#' @param discontinuous.symbol Single character representing the gap when
#' `discontinuous = TRUE`. **Default:** `"."`.
#' @param nthreads Integer number of OpenMP threads to use. `1` forces serial 
#' execution. **Default:** `1`.
#'
#' @return A `data.frame` with two columns:
#' \describe{
#'   \item{motif}{Motif string (contiguous or gapped).}
#'   \item{frequency}{Integer occurrence count across all sequences.}
#' }
#'
#' @examples
#' seqs <- c("CASSLGQDTQYF", "CASSAGQDTQYF", "CASSLGEDTQYF")
#' calculateMotif(seqs, motif.lengths = 3, min.depth = 2)
#'
#'
#' @export
#' @importFrom Rcpp evalCpp
calculateMotif <- function(input.sequences,
                           motif.lengths      = 2:5,
                           min.depth          = 3,
                           discontinuous      = FALSE,
                           discontinuous.symbol = ".",
                           nthreads           = 1) {
  
  input.sequences <- as.character(input.sequences)
  if (length(input.sequences) == 0L)
    stop("`input.sequences` is empty.")
  
  if (any(motif.lengths <= 0))
    stop("`motif.lengths` must be positive integers.")
  
  if (!is.numeric(min.depth) || length(min.depth) != 1L ||
      min.depth < 1 || min.depth != as.integer(min.depth))
    stop("`min.depth` must be a single positive integer (>= 1).", call. = FALSE)
  
  if (!is.character(discontinuous.symbol) || length(discontinuous.symbol) != 1L ||
      nchar(discontinuous.symbol) != 1L)
    stop("`discontinuous.symbol` must be a single character.", call. = FALSE)
  
  res <- calculateMotif_cpp(input.sequences,
                            as.integer(motif.lengths),
                            discontinuous,
                            substr(discontinuous.symbol, 1L, 1L),
                            as.integer(nthreads))
  
  # Convert C++ list → data.frame and apply depth filter
  df <- as.data.frame(res, stringsAsFactors = FALSE)
  if (min.depth > 1)
    df <- df[df$frequency >= min.depth, , drop = FALSE]
  
  rownames(df) <- NULL
  df
}
