#' Adjacency Matrix From Amino Acid or Nucleotide Sequences
#'
#' Calculate frequency of adjacency between residues
#' along a set of biological sequences.
#'
#' @examples
#' # new.sequences <- generateSequences(prefix.motif = "CAS",
#' #                                      suffix.motif = "YF",
#' #                                      number.of.sequences = 100,
#' #                                      min.length = 8,
#' #                                      max.length = 16)
#' #
#' # adj.matrix <- adjacencyMatrix(new.sequences,
#' #                               normalize = TRUE)
#'
#' @param input.sequences Character vector of sequences (amino acid or
#' nucleotide)
#' @param normalize Return the values as a normalized frequency (TRUE) or raw
#' counts (FALSE).
#' @param sequence.dictionary The letters to use in the matrix
#' (defaults to a standard 20 amino acids).
#' @param directed Logical; if FALSE (default) the matrix is symmetrised.
#'
#' @export
#' @importFrom stats setNames
#' @return An adjacency matrix.

adjacencyMatrix <- function(input.sequences,
                            normalize = TRUE,
                            sequence.dictionary = amino.acids,
                            directed  = FALSE) {
  
  # Preflight checks-----------------------------------------------------------
  stopifnot(is.character(input.sequences),
            is.character(sequence.dictionary),
            length(sequence.dictionary) >= 2L)
  
  dict <- unique(sequence.dictionary)
  dict_lookup <- setNames(seq_along(dict), dict)
  M <- length(dict)
  
  # Handle Edge Case: Empty Input 
  if (length(input.sequences) == 0L) {
    return(matrix(0, nrow = M, ncol = M,
                  dimnames = list(dict, dict)))
  }
  
  # Process Sequences Individually 
  # This avoids creating false adjacencies between sequences.
  all_adjacencies <- lapply(input.sequences, function(seq) {
    if (nchar(seq) < 2) return(NULL) # Skip sequences too short to have pairs
    
    chars <- strsplit(seq, "", fixed = TRUE)[[1]]
    
    # Check for any characters not in the dictionary
    bad_chars <- setdiff(chars, dict)
    if (length(bad_chars) > 0) {
      stop("Unknown letters found: ", paste(unique(bad_chars), collapse = ", "))
    }
    
    idx <- dict_lookup[chars]
    
    # Create from->to pairs for this sequence
    list(from = idx[-length(idx)], to = idx[-1L])
  })
  
  # Remove NULLs from sequences with < 2 characters
  all_adjacencies <- all_adjacencies[!sapply(all_adjacencies, is.null)]
  
  if (length(all_adjacencies) == 0L) {
    stop("All sequences have length < 2 - no adjacencies to compute.")
  }
  
  # Combine all 'from' and 'to' vectors
  from_vec <- unlist(sapply(all_adjacencies, `[[`, "from"))
  to_vec   <- unlist(sapply(all_adjacencies, `[[`, "to"))
  
  # Tabulate Adjacencies
  bins <- (to_vec - 1L) * M + from_vec
  counts <- tabulate(bins, nbins = M * M)
  
  adj_matrix <- matrix(counts, nrow = M, ncol = M, byrow = FALSE, # byrow must be FALSE with this indexing
                       dimnames = list(dict, dict))
  
  
  # Make Symmetrical (if requested) 
  if (!directed) {
    adj_matrix <- adj_matrix + t(adj_matrix)
  }
  
  # Normalize (if requested) 
  if (normalize) {
    total <- sum(adj_matrix)
    if (total == 0) {
      warning("Adjacency matrix is all zeros.")
      return(adj_matrix) 
    }
    adj_matrix <- adj_matrix / total
  }
  
  return(adj_matrix)
}