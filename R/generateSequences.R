#' Randomly Generate Amino Acid Sequences
#' 
#' Use this to make synthetic amino acid sequences
#' for purposes of testing code, training models,
#' or noise.
#' 
#' @examples
#' generateSequences(prefix.motif = "CAS",
#'                   suffix.motif = "YF",
#'                   number.of.sequences = 100,
#'                   min.length = 8,
#'                   max.length = 16)
#'                         
#' @param prefix.motif Add defined amino acid/nucleotide sequence to the start of the generated sequences.
#' @param suffix.motif Add defined amino acid/nucleotide sequence to the end of the generated sequences
#' @param number.of.sequences Number of sequences to generate
#' @param min.length Minimum length of the final sequence. The min.length may be adjusted if 
#' incongruent with prefix.motif/suffix.motif lengths
#' @param max.length Maximum length of the final sequence
#' @param sequence.dictionary The letters to use in sequence generation (default are all amino acids)
#' 
#' @export
#' @importFrom stringi stri_rand_strings
#' 
#' @return A vector of generated sequences
generateSequences <- function(prefix.motif = NULL,
                              suffix.motif = NULL,
                              number.of.sequences = 100,
                              min.length = 1,
                              max.length = 10,
                              sequence.dictionary = amino.acids) {
  #TODO Add check for prefix/suffix motif with sequence dictionary
  if(!is.null(prefix.motif) | !is.null(suffix.motif)) {
    motif.length <- sum(nchar(c(prefix.motif,suffix.motif)))
    if(min.length < motif.length) {
      message("New min.length = ", motif.length, " based on the selection of motifs.")
      min.length <- 0
    }
    max.length <- max.length-motif.length
  }
  length.range <- seq(min.length, max.length)
  seq.length <- sample(length.range, number.of.sequences, replace = TRUE)
  
  # Generate the random sequences
  lapply(seq_len(number.of.sequences), function(x) {
      sequence <- stri_rand_strings(n = 1, length = seq.length[x], pattern = paste0("[", paste(sequence.dictionary, collapse = ""), "]"))
      #add prefix and suffix
      if(!is.null(prefix.motif) | !is.null(suffix.motif)) {
        sequence <- paste(prefix.motif, sequence, suffix.motif, sep = "")
      }
      sequence
  }) -> sequences
  sequences <- unlist(sequences)
  return(sequences)
}
