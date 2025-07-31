#' Generate Tokenized Sequences from Amino Acid String
#' 
#' Use this to transform amino acid sequences into tokens in preparing for 
#' deep learning models. 
#' 
#' @examples
#' new.sequences <- generateSequences(prefix.motif = "CAS",
#'                                    suffix.motif = "YF",
#'                                    number.of.sequences = 100,
#'                                    min.length = 8,
#'                                    max.length = 16)
#'                           
#' sequence.matrix <- tokenizeSequences(new.sequences, 
#'                                     add.startstop = TRUE,
#'                                     start.token = "!",
#'                                     stop.token = "^", 
#'                                     convert.to.matrix = TRUE)
#'                         
#' @param input.sequences The amino acid or nucleotide sequences to use
#' @param add.startstop Add start and stop tokens to the sequence
#' @param start.token The character to use for the start token
#' @param stop.token The character to use for the stop token
#' @param max.length Additional length to pad, NULL will pad sequences 
#' to the max length of input.sequences
#' @param convert.to.matrix Return a matrix (TRUE) or a vector (FALSE)
#' @param padding.symbol Single character used for right-padding. 
#' @param verbose Print messages corresponding to the processing step
#' 
#' @export
#' @importFrom stats setNames
#' @return Integer matrix (rows = sequences, cols = positions) or list of vectors.
tokenizeSequences <- function(input.sequences, 
                              add.startstop = TRUE,
                              start.token = "!",
                              stop.token = "^", 
                              max.length = NULL,
                              convert.to.matrix = TRUE,
                              padding.symbol = NULL,
                              verbose = TRUE) {
  # Preflight checks-----------------------------------------------------------
  if (!length(input.sequences))
    return(if (convert.to.matrix) 
      matrix(integer(0), nrow = 0, ncol = 0) else list())
  
  # Build character set -------------------------------------------------------
  if (add.startstop) {
    sequences <- paste0(start.token, input.sequences, stop.token)
    char_set  <- c(start.token, amino.acids, stop.token)
    if (verbose) message("Added start/stop tokens...")
  } else {
    sequences <- input.sequences
    char_set  <- amino.acids
  }
  if (any(nchar(char_set) != 1))
    stop("All tokens in `char_set` must be single characters.")
  
  char_set <- unique(char_set)                 # safeguard 
  char_to_int <- setNames(seq_along(char_set), char_set)
  
  # Length bookkeeping --------------------------------------------------------
  lens <- nchar(sequences)
  if (is.null(max.length)) max.length <- max(lens)
  if (max(lens) > max.length)
    stop("`max.length` is smaller than the longest sequence.")
  
  if (verbose) message("Padding Sequences...")
  pad.id <- padding.symbol %||% (length(char_set) + 1L)
  
  # Tokenise and pad in one sweep ---------------------------------------------
  if (convert.to.matrix) {
    if (verbose) message("Converting to Matrix...")
    N <- length(sequences)
    mat <- matrix(pad.id, nrow = N, ncol = max.length) # integer matrix
    
    for (i in seq_len(N)) {
      ints <- char_to_int[strsplit(sequences[i], "", fixed = TRUE)[[1L]]]
      if (anyNA(ints))
        stop("Unknown character in sequence ", i, call. = FALSE)
      mat[i, seq_along(ints)] <- ints
    }
    return(mat)
  }
  
  ## list-mode ----------------------------------------------------------------
  lapply(sequences, function(seq) {
    ints <- char_to_int[strsplit(seq, "", fixed = TRUE)[[1L]]]
    if (length(ints) < max.length)
      ints <- c(ints, rep(pad.id, max.length - length(ints)))
    ints
  })
}




