#' Decode Amino Acid or Nucleotide Sequences
#'
#' Transforms one-hot or property-encoded sequences back into their original
#' character representation. This function serves as the inverse to
#' `sequenceEncoder`.
#'
#' @param encoded.object A `list` object produced by `sequenceEncoder`, or a
#'   numeric `matrix` (flattened 2D) or `array` (3D cube) from it.
#' @param mode The encoding mode used for decoding: `"onehot"` or `"property"`.
#'   This is typically inferred if `encoded.object` is a list from `sequenceEncoder`.
#' @param property.set For `mode = "property"`, a character vector of property
#'   names (e.g., `"atchleyFactors"`) that were used for the original encoding.
#'   See `?sequenceEncoder`. This is ignored if `property.matrix` is supplied.
#' @param property.matrix For `mode = "property"`, the exact numeric matrix
#'   (with dimensions `20 x P`) that was used for encoding. This overrides
#'   `property.set`.
#' @param call.threshold A numeric confidence threshold for making a call.
#'   - In `"onehot"` mode, this is the minimum required value in the vector (e.g., `0.9`).
#'   - In `"property"` mode, this is the maximum allowable Euclidean distance.
#'   Positions with scores not meeting the threshold are assigned the `padding.symbol`.
#' @param sequence.dictionary A character vector of the alphabet (e.g., amino acids).
#'   Must match the one used during encoding.
#' @param padding.symbol The single character used to represent padding or
#'   low-confidence positions.
#' @param remove.padding Logical. If `TRUE`, trailing padding symbols are
#'   removed from the end of the decoded sequences.
#' @return A character vector of the decoded sequences.
#'
#' @examples
#' # Example sequences
#' aa.sequences <- c("CAR", "YMD", "ACAC")
#'
#' # Encode the sequences
#' encoded.onehot <- sequenceEncoder(aa.sequences, 
#'                                   mode = "onehot")
#' encoded.prop <- sequenceEncoder(aa.sequences, 
#'                                 mode = "property", 
#'                                 property.set = "atchleyFactors")
#'
#' # Decode the sequences
#' # 1. Decode from the full list object
#' decoded.1 <- sequenceDecoder(encoded.onehot, 
#'                              mode = "onehot")
#'
#' # 2. Decode from just the 3D cube array
#' decoded.2 <- sequenceDecoder(encoded.prop$cube,
#'                              mode = "property",
#'                              property.set = "atchleyFactors")
#'
#'
#' @export
sequenceDecoder <- function(encoded.object,
                            mode = c("onehot", "property"),
                            property.set = NULL,
                            property.matrix = NULL,
                            call.threshold = 0.5,
                            sequence.dictionary = amino.acids,
                            padding.symbol = ".",
                            remove.padding = TRUE) {
  
  # Preflight checks-----------------------------------------------------------
  if (is.list(encoded.object) && all(c("cube", "sequence.dictionary") %in% names(encoded.object))) {
    if (missing(sequence.dictionary) && !is.null(encoded.object$sequence.dictionary)) {
      sequence.dictionary <- encoded.object$sequence.dictionary
    }
    if (missing(padding.symbol) && !is.null(encoded.object$pad_token)) {
      padding.symbol <- encoded.object$pad_token
    }
    cube <- encoded.object$cube
  } else if (is.array(encoded.object) && length(dim(encoded.object)) == 3) {
    cube <- encoded.object
    mode <- match.arg(mode)
  } else if (is.matrix(encoded.object)) {
    mode <- match.arg(mode)
    n_seq <- nrow(encoded.object)
    
    if (mode == "onehot") {
      depth <- length(c(sequence.dictionary, padding.symbol))
    } else {
      if (!is.null(property.matrix)) {
        depth <- nrow(property.matrix)
      } else if (!is.null(property.set)) {
        depth <- nrow(.aa.property.matrix(property.set))
      } else {
        stop("For flattened matrix input in 'property' mode, supply 'property.set' or 'property.matrix'.")
      }
    }
    
    if ((ncol(encoded.object) %% depth) != 0) {
      stop("Cannot reshape flattened matrix: total columns not a multiple of encoding depth.")
    }
    max_len <- ncol(encoded.object) / depth
    cube <- array(t(encoded.object), dim = c(depth, max_len, n_seq))
  } else {
    stop("'encoded.object' must be a list from sequenceEncoder, a 3D array, or a 2D matrix.")
  }
  
  if (call.threshold <= 0) {
    stop("'call.threshold' must be a positive number.")
  }
  
  # 1) Dispatch to Efficient Decoder 
  if (mode == "onehot") {
    decoded_sequences <- .onehotDecoder(cube, sequence.dictionary, padding.symbol, call.threshold)
  } else if (mode == "property") {
    if (is.null(property.matrix)) {
      if (is.null(property.set)) {
        stop("In 'property' mode, you must supply either 'property.set' or 'property.matrix'.")
      }
      property.matrix <- .aa.property.matrix(property.set)
    }
    if (ncol(property.matrix) != length(sequence.dictionary)) {
      stop("Rows in 'property.matrix' must match the length of 'sequence.dictionary'.")
    }
    decoded_sequences <- .propertyDecoder(cube, property.matrix, sequence.dictionary, padding.symbol, call.threshold)
  }
  
  # 2) Final Processing 
  if (remove.padding) {
    decoded_sequences <- sub(paste0("\\", padding.symbol, "*$"), "", decoded_sequences)
  }
  
  return(decoded_sequences)
}


# Efficient One-Hot Decoder 
.onehotDecoder <- function(cube, 
                           sequence.dictionary, 
                           padding.symbol, 
                           call.threshold) {
  permuted_cube <- aperm(cube, c(3, 2, 1))
  
  # Get the index of the max value for each position (returns n_seq x max_len matrix)
  max_indices <- apply(permuted_cube, c(1, 2), which.max)
  
  # Get the max value and confident calls
  max_values <- apply(permuted_cube, c(1, 2), max)
  is_unique_max <- apply(permuted_cube, c(1, 2), function(vec) sum(vec == max(vec)) == 1)
  is_confident <- (max_values >= call.threshold) & is_unique_max
  
  # Create a matrix of characters
  char_matrix <- matrix(padding.symbol, nrow = nrow(is_confident), ncol = ncol(is_confident))
  char_matrix[is_confident] <- c(sequence.dictionary, padding.symbol)[max_indices[is_confident]]
  
  # Collapse each row into a final sequence string
  apply(char_matrix, 1, function(x) {
    paste0(x[!is.na(x)], collapse = "")
  })
}


# Efficient Property Decoder 
#' @importFrom stats dist
.propertyDecoder <- function(cube, 
                             property.matrix, 
                             sequence.dictionary, 
                             padding.symbol, 
                             call.threshold) {
  ref_mat <- property.matrix
  
  # Permute cube for easier iteration: n_seq x max_len x P (depth)
  permuted_cube <- aperm(cube, c(3, 2, 1))
  
  # Pre-calculate squared column sums of the reference matrix for speed
  ref_mat_col_sq_sums <- colSums(ref_mat^2)
  
  # Apply decoding logic to each sequence matrix (max_len x P)
  apply(permuted_cube, 1, function(seq_matrix) {
    chars <- apply(seq_matrix, 1, function(pos_vector) {
      if (all(pos_vector == 0)) {
        return(padding.symbol)
      }
      combined_vectors <- rbind(pos_vector, t(ref_mat))
      distance_matrix <- as.matrix(dist(combined_vectors))
      distances <- distance_matrix[1, -1]
      min_dist <- min(distances)
      if (min_dist <= call.threshold) {
        return(sequence.dictionary[which.min(distances)])
      } else {
        return(padding.symbol)
      }
    })
    paste0(chars, collapse = "")
  })
}