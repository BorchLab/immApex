#' Randomly Mutate Sequences of Amino Acids
#' 
#' Use this to mutate or mask sequences for purposes of testing code, 
#' training models, or noise.
#' 
#' @examples
#' sequences <- generateSequences(prefix.motif = "CAS",
#'                                suffix.motif = "YF",
#'                                number.of.sequences = 100,
#'                                min.length = 8,
#'                                max.length = 16)
#'                                 
#' mutated_sequences <- mutateSequences(sequences, 
#'                                      number.of.sequences = 1,
#'                                      position.start = 3,
#'                                      position.end = 8)
#' 
#' @param input.sequences The amino acid or nucleotide sequences to use
#' @param number.of.sequences The number of mutated sequences to return
#' @param mutation.rate The rate of mutations to introduce into sequences
#' @param position.start The starting position to mutate along the sequence
#' \strong{Default} = NULL will start the random mutations at position 1
#' @param position.end The ending position to mutate along the sequence
#' \strong{Default} = NULL will end the random mutations at the last position
#' @param sequence.dictionary The letters to use in sequence mutation
#' (default are all amino acids)
#' 
#' @export 
#' @return A vector of mutated sequences
mutateSequences <- function(input.sequences, 
                            number.of.sequences = 1, 
                            mutation.rate = 0.01,
                            position.start = NULL,
                            position.end = NULL,
                            sequence.dictionary = amino.acids) {
  
  # Preflight checks-----------------------------------------------------------
  stopifnot(is.character(input.sequences),
            length(number.of.sequences) == 1L,
            number.of.sequences >= 0L,
            is.numeric(mutation.rate),
            mutation.rate >= 0, mutation.rate <= 1)
  
  dict <- unique(as.character(sequence.dictionary))
  bad  <- sprintf("[^%s]", paste(dict, collapse = ""))
  if (any(grepl(bad, input.sequences)))
    stop("`input.sequences` contain letters not in `sequence.dictionary`.")
  
  # Internal helper -----------------------------------------------------------
  mutate_one <- function(seq) {
    
    L <- nchar(seq)
    if (L == 0 || mutation.rate == 0) return(seq)
    
    # mutation window ---------------------------------------------------------
    p0 <- if (is.null(position.start)) 1L else max(1L, position.start)
    p1 <- if (is.null(position.end))   L  else min(L, position.end)
    if (p0 > p1) stop("position.start > position.end")
    
    # how many mutations ------------------------------------------------------
    m <- ceiling((p1 - p0 + 1L) * mutation.rate)
    if (m == 0) return(seq)  
    
    pos <- sample.int(p1 - p0 + 1L, m, replace = FALSE) + (p0 - 1L)
    
    # split once, mutate in place ---------------------------------------------
    aa <- strsplit(seq, "", fixed = TRUE)[[1L]]
    for (idx in pos) {
      new <- sample(dict, 1L)
      while (new == aa[idx]) new <- sample(dict, 1L)
      aa[idx] <- new
    }
    paste0(aa, collapse = "")
  }
  
  ## Generate all mutants -----------------------------------------------------
  originals <- rep(input.sequences, each = number.of.sequences)
  mutated.sequences   <- vapply(originals, mutate_one, FUN.VALUE = character(1))
  return(mutated.sequences)
}
  
