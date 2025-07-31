#' Positional Entropy / Diversity Biological Sequences
#'
#' @description
#' Computes residue-wise diversity for a set of aligned (right-padded)
#' CDR3 amino-acid sequences using *any* supported diversity estimator
#' in **immApex**.  The following metrics are recognized:
#'
#' * **Shannon entropy:**           \code{\link[=shannon_entropy]{shannon_entropy}}
#' * **Inverse Simpson:**           \code{\link[=inv_simpson]{inv_simpson}}
#' * **Gini–Simpson index:**        \code{\link[=gini_simpson]{gini_simpson}}
#' * **Normalized entropy:**        \code{\link[=norm_entropy]{norm_entropy}}
#' * **Pielou evenness:**           \code{\link[=pielou_evenness]{pielou_evenness}}
#' * * **Hill numbers** (orders 0, 1, 2): \code{\link[=hill_q]{hill_q(0)}}, 
#' \code{\link[=hill_q]{hill_q(1)}}, \code{\link[=hill_q]{hill_q(2)}}
#'
#' You may also supply a **custom function** to `method`; it must take a
#' numeric vector of clone counts and return a single numeric value.
#'
#' @param input.sequences `character()`. Vector of CDR3 AA strings.
#' @param max.length `integer(1)`. Target length to align / pad to.
#' *Default* = `max(nchar(sequences))`.
#' @param method Either the name of a built-in metric  (`"shannon"`, 
#' `"inv.simpson"`, `"gini.simpson"`, `"norm.entropy"`, `"pielou"`, `"hill0"`, 
#' `"hill1"`, `"hill2"`) **or** a custom function as described above.
#' @param padding.symbol Symbol to use for padding at the end of sequences.
#'
#' @return Named `numeric()` vector of diversity scores,
#'         one value per position (Pos1 … Pos*L*).
#'         
#' @examples
#' seqs <- c("CASSLGQDTQYF", "CASSIRSSYNEQFF", "CASSTGELFF")
#' calculateEntropy (seqs, method = "shannon")
#' @export
calculateEntropy <- function(input.sequences,
                             max.length = NULL,
                             method    = c("shannon", 
                                           "inv.simpson", 
                                           "gini.simpson", 
                                           "norm.entropy",
                                           "pielou", 
                                           "hill0", 
                                           "hill1", 
                                           "hill2"),
                             padding.symbol = ".") {
  # Preflight checks-----------------------------------------------------------
  stopifnot(is.character(input.sequences))
  if (is.null(max.length))
    max.length <- max(nchar(input.sequences))
  if (nchar(padding.symbol) != 1L)
    stop("'padding.symbol' must be a single character")
  
  method <- match.arg(method)
  
  # 1. Pad sequences to equal length and split into a char matrix
  pad_seq  <- paste0(
    input.sequences,
    vapply(max.length - nchar(input.sequences),
           function(x) if (x > 0) strrep(padding.symbol, x) else "",
           character(1))
  )
  mat <- matrix(
    unlist(strsplit(pad_seq, "")),
    ncol = max.length, byrow = TRUE,
    dimnames = list(NULL, paste0("Pos", seq_len(max.length)))
  )
  
  ## pick the scoring function 
  div_fun <- if (is.function(method)) method else {
    if (!(method %in% names(.div.registry)))
      stop("Unknown method; choose one of ",
           paste(names(.div.registry), collapse = ", "),
           " or supply a function.")
    .div.registry[[method]]
  }
  
  ## 3. Vectorised per-column calculation 
  res <- vapply(seq_len(max.length), function(i) {
    cnt <- table(mat[, i])
    cnt <- cnt[names(cnt) != padding.symbol]             # drop padding
    if (length(cnt) <= 1L) return(0)          # no variability
    div_fun(cnt)
  }, numeric(1L))
  
  names(res) <- colnames(mat)
  res
}
