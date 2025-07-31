#' Randomly Generate Amino Acid Sequences
#'
#' Use this to make synthetic amino acid sequences for purposes of testing code,
#' training models, or providing noise.
#'
#' @param prefix.motif A defined amino acid/nucleotide sequence to add to the
#'   start of the generated sequences.
#' @param suffix.motif A defined amino acid/nucleotide sequence to add to the
#'   end of the generated sequences.
#' @param number.of.sequences The number of sequences to generate.
#' @param min.length The minimum length of the final sequence. If this value is
#'   too short to fit the motifs, it will be automatically increased.
#' @param max.length The maximum length of the final sequence. If it is less
#'   than the final `min.length`, it will also be adjusted.
#' @param sequence.dictionary A character vector of the letters to use in
#'   random sequence generation.
#' @param verbose Logical. If TRUE, prints messages when arguments like
#'   `min.length` or `max.length` are automatically adjusted.
#'   
#' @examples
#' generateSequences(prefix.motif = "CAS",
#'                   suffix.motif = "YF",
#'                   number.of.sequences = 100,
#'                   min.length = 8,
#'                   max.length = 16)
#'
#' @export
#' @return A character vector of generated sequences.
generateSequences <- function(prefix.motif = NULL,
                              suffix.motif = NULL,
                              number.of.sequences = 100,
                              min.length = 1,
                              max.length = 10,
                              verbose = TRUE, 
                              sequence.dictionary = amino.acids) {
  
  # Preflight checks-----------------------------------------------------------
  prefix.motif <- if (is.null(prefix.motif)) "" else as.character(prefix.motif)
  suffix.motif <- if (is.null(suffix.motif)) "" else as.character(suffix.motif)
  
  if (!is.numeric(number.of.sequences) || length(number.of.sequences) != 1L ||
      number.of.sequences < 0 || number.of.sequences %% 1 != 0) {
    stop("`number.of.sequences` must be a single non-negative integer.")
  }
  if (number.of.sequences == 0) return(character(0))
  
  # 1. Automatic Length Adjustment with User Notification
  motif_len <- nchar(prefix.motif) + nchar(suffix.motif)
  
  if (min.length < motif_len) {
    if (verbose) {
      message(paste("`min.length` of", min.length, "is too short for motifs (total length",
                    motif_len, "). Adjusting min.length to", motif_len))
    }
    min.length <- motif_len
  }
  
  if (max.length < min.length) {
    if (verbose) {
      message(paste("`max.length` of", max.length, "is less than the final min.length of",
                    min.length, "). Adjusting max.length to", min.length))
    }
    max.length <- min.length
  }
  
  # Final validation check after adjustments
  all_motif_chars <- c(strsplit(prefix.motif, "")[[1]], strsplit(suffix.motif, "")[[1]])
  dict <- unique(as.character(sequence.dictionary))
  if (any(!all_motif_chars %in% dict)) {
    bad_chars <- unique(all_motif_chars[!all_motif_chars %in% dict])
    stop("Prefix/suffix contain letters not in `sequence.dictionary`: ",
         paste(bad_chars, collapse = ", "))
  }
  
  # 3. Generate Sequence Lengths 
  range_len <- min.length:max.length
  if(length(range_len) > 1) {
    target_lengths <- sample(range_len, number.of.sequences, replace = TRUE)
  } else {
    target_lengths <- rep(range_len, number.of.sequences)
  }
  random_part_lengths <- target_lengths - motif_len
  
  # 4. Generate Random Sequences 
  random_pieces <- vapply(random_part_lengths, function(n) {
    paste0(sample(dict, n, replace = TRUE), collapse = "")
  }, FUN.VALUE = character(1))
  
  # 5. Combine and Return 
  sequences <- paste0(prefix.motif, random_pieces, suffix.motif)
  
  return(sequences)
}