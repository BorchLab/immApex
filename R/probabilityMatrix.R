#' Position Probability Matrix for Amino Acid or Nucleotide Sequences
#'
#' Generates a position-probability (PPM) or position-weight (PWM) matrix
#' from a set of biological sequences.
#'
#' @param input.sequences Character vector of sequences.
#' @param max.length Integer; sequences will be right-padded to this length. 
#'   If NULL (default), pads to the length of the longest sequence in the input.
#' @param convert.PWM Logical; if TRUE, converts the matrix into a PWM.
#' @param background.frequencies Named vector of background frequencies for PWM 
#'   calculation. If NULL, a uniform distribution is assumed. Names must 
#'   correspond to characters in `sequence.dictionary`.
#' @param sequence.dictionary Character vector of residues to include in the matrix.
#' @param pseudocount A small number added to raw counts for PWM calculation to 
#'   avoid zero probabilities. Defaults to 1.
#' @param padding.symbol Single character for right-padding. Must not be in `sequence.dictionary`.
#' 
#' @examples
#' new.sequences <- generateSequences(prefix.motif = "CAS",
#'                                    suffix.motif = "YF",
#'                                    number.of.sequences = 100,
#'                                    min.length = 8,
#'                                    max.length = 16)
#'                           
#' PPM.matrix <- probabilityMatrix(new.sequences)
#'
#' @export
#' @return A matrix with position-specific probabilities (PPM) or weights (PWM).
probabilityMatrix <- function(input.sequences,
                              max.length = NULL,
                              convert.PWM = FALSE,
                              background.frequencies = NULL,
                              sequence.dictionary = amino.acids,
                              pseudocount = 1,
                              padding.symbol = ".") {
  
  #  1. Preflight Checks & Setup --------------------------------------------
  if (padding.symbol %in% sequence.dictionary) {
    stop("`padding.symbol` cannot be present in `sequence.dictionary`.")
  }
  if (length(input.sequences) == 0) {
    return(matrix(0, nrow = length(sequence.dictionary), ncol = 0, 
                  dimnames = list(sequence.dictionary, NULL)))
  }
  
  # Determine padding length
  if (is.null(max.length)) {
    max.length <- max(nchar(input.sequences))
  }
  
  # 2. Pad Sequences 
  padded_sequences <- .padded.strings(input.sequences, 
                                      max.length = max.length, 
                                      pad = ".", collapse = TRUE)
  sequence_matrix <- do.call(rbind, strsplit(padded_sequences, ""))
  
  # 3. Count Occurrences 
  count_matrix <- matrix(0, nrow = length(sequence.dictionary), ncol = max.length)
  rownames(count_matrix) <- sequence.dictionary
  colnames(count_matrix) <- paste0("Pos.", seq_len(max.length))
  
  # Vectorized counting for each character
  for (letter in sequence.dictionary) {
    count_matrix[letter, ] <- colSums(sequence_matrix == letter)
  }
  
  #  4. PPM Calculation (Default) 
  # Calculate column totals, excluding pads, for correct normalization
  col_totals <- colSums(count_matrix)
  
  # Avoid division by zero for columns that are all padding
  col_totals[col_totals == 0] <- 1 
  
  # Normalize to get probabilities
  prob_matrix <- sweep(count_matrix, 2, col_totals, FUN = "/")
  
  # 5. PWM Conversion 
  if (convert.PWM) {
    # Validate and prepare background frequencies
    if (is.null(background.frequencies)) {
      # Assume uniform background if not provided
      background.frequencies <- rep(1 / length(sequence.dictionary), length(sequence.dictionary))
      names(background.frequencies) <- sequence.dictionary
    } else {
      if (!setequal(names(background.frequencies), sequence.dictionary)) {
        stop("Names of `background.frequencies` must match `sequence.dictionary`.")
      }
    }
    
    # Add pseudocounts to the *raw counts*
    smoothed_counts <- count_matrix + pseudocount
    
    # Recalculate column totals with pseudocounts
    smoothed_col_totals <- colSums(smoothed_counts)
    
    # Re-normalize to get smoothed probabilities
    smoothed_probs <- sweep(smoothed_counts, 2, smoothed_col_totals, FUN = "/")
    
    # Calculate log-likelihood scores, ensuring background vector alignment
    prob_matrix <- log2(smoothed_probs / background.frequencies[rownames(smoothed_probs)])
  }
  
  return(prob_matrix)
}
