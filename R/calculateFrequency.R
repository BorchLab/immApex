#' Relative Residue Frequencies at Every Position
#'
#' Quickly computes the per-position relative frequency of each symbol
#' (amino-acid or nucleotide) in a set of biological sequences.  Variable-length
#' strings are padded to a common width so the calculation is entirely
#' vectorized (one logical comparison + one `colSums()` per residue).
#'
#' @param input.sequences  Character vector of sequences (amino acid or 
#' nucleotide)
#' @param max.length Integer.  Pad/trim to this length. Defaults to
#'  `max(nchar(sequences))`.
#' @param sequence.dictionary Vector of valid residue symbols that should be
#' tracked (defaults to the 20 canonical amino acids; supply
#' `c("A","C","G","T","N")` etc. for nucleotides).
#' @param padding.symbol Single character used for right-padding. **Must not**
#' be present in `sequence.dictionary`.
#' @param summary.fun Character string choosing the summary statistic:
#'   * `"proportion"` (default) – each cell sums to 1 over the table.  
#'   * `"count"`      – raw counts.  
#'   * `"percent"`    – proportion × 100.
#' @param tidy Logical; if `TRUE` a long-format `data.frame` is returned
#' instead of a matrix (useful for plotting with *ggplot2*).
#'
#' @return Either
#' \itemize{
#'   \item A numeric matrix of dimension `length(sequence.dictionary)` ×
#'         `max.length`, whose columns sum to 1, **or**
#'   \item A `data.frame` with columns *position*, *residue*, *frequency* when
#'         `tidy = TRUE`.
#' }
#' @importFrom stats na.omit
#' @examples
#' # Amino Acid example
#' seqs <- c("CASSLGQGAETQYF", "CASSPGQGDYEQYF", "CASSQETQYF")
#' rel.freq <- calculateFrequency(seqs)
#' head(rel.freq[, 1:5])                  
#'
#' # Nucleotide example
#' dna <- c("ATGCC", "ATGAC", "ATGGC")
#' calculateFrequency(dna,
#'                    sequence.dictionary = c("A","C","G","T"),
#'                    padding.symbol = "-",
#'                    tidy = TRUE)
#'
#' @export
calculateFrequency <- function(input.sequences,
                               max.length = NULL,
                               sequence.dictionary = amino.acids,
                               padding.symbol = ".",
                               summary.fun = c("proportion", "count", "percent"),
                               tidy = FALSE) {
  # Preflight checks-----------------------------------------------------------
  stopifnot(is.character(input.sequences),
            nchar(padding.symbol) == 1L,
            padding.symbol %!in% sequence.dictionary)
  summary.fun <- match.arg(summary.fun)
  
  # Handling multiple sequences per scRep environment
  if(any(grep(";", input.sequences))) {
    input.sequences <- unlist(strsplit(input.sequences, ";"))
  }
  
  # Remove NAs
  input.sequences <- na.omit(input.sequences)
  # Remove empty strings 
  input.sequences <- input.sequences[nchar(input.sequences) > 0]
  
  # 1. Pad to a rectangular character matrix 
  if (is.null(max.length))
    max.length <- max(nchar(input.sequences), 1L)
  
  padded <- .padded.strings(input.sequences,
                            max.length = max.length,
                            pad = padding.symbol,
                            collapse  = TRUE)
  
  seq_mat <- as.matrix(do.call(rbind, strsplit(unlist(padded), ""))[,seq_len(max.length)])  # nSeq × max.length
  nSeq <- length(padded)
  
  # 2.  Fast frequency calculation (one pass per residue) 
  res_mat <- matrix(0,
                    nrow = length(sequence.dictionary) + 1,
                    ncol = max.length,
                    dimnames = list(c(sequence.dictionary, padding.symbol),
                                    paste0("Pos.", seq_len(max.length))))
  
  denominator <- nrow(seq_mat)
  for (residue in c(sequence.dictionary, padding.symbol)) {
    # logical comparison is vectorized; colSums is C-level
    res_mat[residue, ] <- .scale.counts(colSums(seq_mat == residue), summary.fun, denominator)
  }
  
  # 3.  Optional tidy reshaping 
  if (tidy) {
    res_mat <- as.data.frame(as.table(res_mat),
                             stringsAsFactors = FALSE,
                             responseName = "frequency")
    names(res_mat) <- c("residue", "position", "frequency")
    res_mat$position <- as.integer(sub("Pos\\.", "", res_mat$position))
  }
  
  res_mat
}
