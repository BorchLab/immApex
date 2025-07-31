#' Generate Sinusoidal Positional Encodings
#'
#' Creates a matrix of sinusoidal positional encodings as described in the
#' "Attention Is All You Need" paper. This provides a way to inject information
#' about the relative or absolute position of tokens in a sequence.
#'
#' @section Details:
#' The implementation uses the standard formulas:
#' `PE(pos, 2i) = sin(pos / base^(2i / d.model))`
#' `PE(pos, 2i+1) = cos(pos / base^(2i / d.model))`
#' where `pos` is the position, `i` is the dimension pair, `d.model` is the
#' embedding dimension, and `base` is a user-definable base, typically 10000.
#'
#' @examples
#' pos_encoding <- positionalEncoder(max.length = 50, 
#'                                   d.model = 64)
#'
#' my_sequences <- c("SEQVENCE", "ANOTHERSEQ")
#' pos_enc_auto <- positionalEncoder(input.sequences = my_sequences, 
#'                                   d.model = 32)
#'
#' @param max.length The maximum sequence length (number of positions) to encode.
#' This is the primary way to specify the output size.
#' @param d.model The dimensionality of the embedding. Must be an even number.
#' @param input.sequences Optional. A character vector of sequences. If provided,
#' `max.length` is automatically determined from the longest sequence,
#' unless `max.length` is also explicitly set to a larger value.
#' @param base The base for the geometric progression of frequencies. The default
#' is 10000, as used in the original paper.
#' @param position.offset An integer offset for position numbering. Defaults to 1
#' (1-based indexing common in R). Set to 0 for 0-based indexing.
#'
#' @export
#' @return A matrix of shape `max.length` x `d.model` containing the positional
#'   encodings.
positionalEncoder <- function(max.length = NULL,
                              d.model = NULL,
                              input.sequences = NULL,
                              base = 10000,
                              position.offset = 1L) {
  
  # Preflight checks-----------------------------------------------------------
  if (is.null(d.model)) {
    stop("`d.model` (embedding dimension) must be specified.")
  }
  if (d.model %% 2 != 0) {
    stop("`d.model` must be an even number.")
  }
  
  # Determine max.length automatically if not provided
  if (is.null(max.length)) {
    if (is.null(input.sequences)) {
      stop("Must provide either `max.length` or a set of `sequences`.")
    }
    max.length <- max(nchar(input.sequences))
  }
  
  # 1. Efficient Calculation of Angle Rads 
  positions <- seq.int(from = position.offset, 
                       to = max.length + position.offset - 1)
  i <- seq.int(from = 0, to = d.model - 2, by = 2)
  denominators <- base^(i / d.model)
  angle_rads <- outer(positions, denominators, FUN = "/")
  
  # 2. Assemble Final Encoding Matrix
  pos_encoding <- matrix(NA, nrow = max.length, ncol = d.model)
  pos_encoding[, seq.int(1, d.model, by = 2)] <- sin(angle_rads)
  pos_encoding[, seq.int(2, d.model, by = 2)] <- cos(angle_rads)
  
  return(pos_encoding)
}