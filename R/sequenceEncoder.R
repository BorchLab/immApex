#' Universal Amino-acid Sequence Encoder 
#'
#' `sequenceEncoder()` is a high-level function that converts a character vector
#' of amino-acid sequences into one of three representations:
#' 1.  **one-hot**: A binary representation for each amino acid position.
#' 2.  **property-based**: A numerical representation based on amino acid properties
#'     (e.g., atchleyFactors, kideraFactors, etc).
#' 3.  **geometric**: A fixed-length 20-dimensional vector for each sequence,
#'     derived from a substitution matrix and geometric rotation.
#'
#' The function acts as a wrapper for either the C++ backend (for one-hot and
#' property modes) or the R-based geometric transformation.
#'
#' @section Property Mode:
#' If you supply `property.matrix` directly, it **must** be a numeric matrix
#' whose **rows correspond to the 20 canonical amino acids in the order of
#' `sequence.dictionary`** and whose columns are the property scales.
#'
#' @section Geometric Mode:
#' This mode projects sequences into a 20D space. It calculates the average
#' vector for each sequence using a substitution matrix (e.g., "BLOSUM62")
#' and then applies a planar rotation to the resulting vector.
#'
#' @param input.sequences `character` vector. Sequences (uppercase
#'   single-letter code).
#' @param mode Either `"onehot"`, `"property"`, or `"geometric"`.
#' @param property.set Character string (one of the supported names) 
#'  Defaults to `"atchleyFactors"`, but includes: `"crucianiProperties"`, 
#' `"FASGAI"`, `"kideraFactors"`, `"MSWHIM"`, `"ProtFP"`, `"stScales"`, 
#' `"tScales"`, `"VHSE"`, `"zScales"` Ignored if `property.matrix` is supplied.
#' @param property.matrix *Optional numeric matrix (`20 × P`)*. Overrides
#'   `property.set` in `"property"` mode.
#' @param method *(For geometric mode)* Character key for a built-in substitution
#'   matrix (e.g., "BLOSUM62"), or a 20x20 numeric matrix itself.
#' @param theta *(For geometric mode)* Rotation angle in radians (default `pi/3`).
#' @param sequence.dictionary Character vector of the alphabet (default = 20
#'   standard amino acids).
#' @param padding.symbol Single character for right-padding (non-geometric modes).
#' @param summary.fun For property mode only: `"mean"` or `""` (none).
#' @param max.length Integer for truncation/padding. If `NULL` (default), the
#'   longest sequence sets the maximum. Not used in geometric mode.
#' @param nthreads Number of threads for C++ backend. Not used in geometric mode.
#' @param verbose Logical. If `TRUE` (default), prints a progress message.
#' @param ... Additional arguments passed to `sequenceEncoder()` when using 
#' wrapper functions (`onehotEncoder`, `propertyEncoder`, `geometricEncoder`).
#'
#' @return A named `list` containing the encoded data and metadata.
#' \describe{
#'   \item{`cube`}{3D Numeric array. `NULL` in geometric mode.}
#'   \item{`flattened`}{2D Numeric matrix. `NULL` in geometric mode.}
#'   \item{`summary`}{2D Numeric matrix containing sequence-level representations.
#'     This is the primary output for geometric mode.}
#'   \item{...}{Other metadata related to the encoding process.}
#' }
#' @importFrom utils data
#'
#' @examples
#' aa <- c("CARDRST", "YYYGMD", "ACACACAC")
#'
#' # One-hot encoding
#' enc_onehot <- sequenceEncoder(aa, 
#'                               mode = "onehot")
#'
#' # Property-based encoding
#' enc_prop <- sequenceEncoder(aa, 
#'                             mode = "property", 
#'                             property.set = "atchleyFactors")
#'
#' # Geometric encoding
#' enc_geo <- sequenceEncoder(aa, 
#'                            mode = "geometric", 
#'                            method = "BLOSUM62")
#'
#' @export
sequenceEncoder <- function(input.sequences,
                            mode             = c("onehot", "property", "geometric"),
                            property.set     = NULL,
                            property.matrix  = NULL,
                            method           = "BLOSUM62",
                            theta            = pi / 3,
                            sequence.dictionary = amino.acids,
                            padding.symbol   = ".",
                            summary.fun      = "",
                            max.length       = NULL,
                            nthreads         = parallel::detectCores(),
                            verbose          = TRUE,
                            ...) {
  
  mode <- match.arg(mode)
  
  if (verbose)
    message(sprintf("[sequenceEncoder] Encoding %d sequence%s (%s mode)...",
                    length(input.sequences), ifelse(length(input.sequences) == 1L, "", "s"), mode))
  
  # Mode: Geometric
  if (mode == "geometric") {
    # Hoist helper function from the original geometricEncoder.R
    fetch_matrix <- function(m) {
      if (is.matrix(m)) {
        if (!all(dim(m) == 20) || is.null(rownames(m))) {
          stop("If `method` is a matrix, it must be 20x20 with amino acid rownames.")
        }
        return(m)
      }
      data("immapex_blosum.pam.matrices", package = "immApex", envir = environment())
      mat <- immapex_blosum.pam.matrices[[m]]
      if (is.null(mat)) stop("Cannot find matrix for method '", m, "'.")
      return(mat)
    }
    
    S <- fetch_matrix(method)
    aa_lookup <- setNames(seq_len(nrow(S)), rownames(S))
    seq_lengths <- nchar(input.sequences)
    group_id <- rep.int(seq_along(seq_lengths), seq_lengths)
    all_chars <- unlist(strsplit(input.sequences, "", fixed = TRUE), use.names = FALSE)
    all_indices <- aa_lookup[all_chars]
    
    if (anyNA(all_indices)) {
      bad_chars <- unique(all_chars[is.na(all_indices)])
      stop("Non-canonical amino acid(s) found: ", paste(bad_chars, collapse = ", "))
    }
    
    summed_vectors <- rowsum(S[all_indices, , drop = FALSE], group_id, reorder = FALSE)
    avg_vectors <- summed_vectors / seq_lengths
    
    R2 <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2, 2)
    R20 <- diag(20)
    for (i in seq(1, 19, 2)) {
      R20[i:(i + 1), i:(i + 1)] <- R2
    }
    rotated_vectors <- avg_vectors %*% t(R20)
    
    return(list(
      cube = NULL,
      flattened = NULL,
      summary = rotated_vectors,
      mode = "geometric",
      method = ifelse(is.character(method), method, "custom matrix")
    ))
  }
  
  # Modes: onehot, property (handled by C++ backend)
  prop_mat <- NULL
  if (mode == "property") {
    if (!is.null(property.set)) {
      prop_mat <- t(.aa.property.matrix(property.set))
    } else if (!is.null(property.matrix)) {
      if (!is.matrix(property.matrix) || !is.numeric(property.matrix) ||
          nrow(property.matrix) != length(sequence.dictionary)) {
        stop(sprintf(
          "`property.matrix` rows (%d) must match `sequence.dictionary` length (%d).",
          nrow(property.matrix), length(sequence.dictionary)
        ))
      }
      prop_mat <- property.matrix
    } else {
      stop("In `property` mode, supply either `property.set` or `property.matrix`.")
    }
  }
  
  if (is.null(max.length))
    max.length <- max(nchar(input.sequences), 1L)
  
  out <- encodeSequences_cpp(
    sequences      = input.sequences,
    mode           = mode,
    alphabet       = sequence.dictionary,
    prop_mat_      = prop_mat,
    pad_token      = padding.symbol,
    summary        = summary.fun,
    max_len        = max.length,
    nthreads       = nthreads
  )
  
  # Add mode to output for clarity
  out$mode <- mode
  
  # Adding dimension names
  feature_names <- NULL
  if (mode == "onehot") {
    # For one-hot, features are the AAs + the padding symbol
    feature_names <- out$sequence.dictionary
    if (!(out$padding.symbol %in% feature_names)) {
      feature_names <- c(feature_names, out$padding.symbol)
    }
  } else if (mode == "property") {
    # For properties, features are the property names (from colnames)
    if (!is.null(prop_mat)) {
      if (!is.null(colnames(prop_mat))) {
        feature_names <- colnames(prop_mat)
      } else {
        # Fallback to generic names if no column names are found
        num_props <- ncol(prop_mat)
        feature_names <- paste0("Prop", seq_len(num_props))
      }
    }
  }
  
  # Naming flattened element
  if (!is.null(out$flattened) && !is.null(feature_names) && max.length > 0) {
    col_names <- as.vector(sapply(seq_len(max.length), function(pos) {
      paste(feature_names, pos, sep = "_")
    }))
    colnames(out$flattened) <- col_names
  }
  
  # Naming summary element
  if (!is.null(out$summary) && !is.null(feature_names)) {
    colnames(out$summary) <- feature_names
  }
  
  # Naming array element
  if (!is.null(out$cube) && !is.null(feature_names) && max.length > 0) {
    # Use names from the input vector if they exist, otherwise create generic ones
    sample_names <- names(input.sequences)
    if (is.null(sample_names)) {
      sample_names <- paste0("S", seq_along(input.sequences))
    }
    dimnames(out$cube) <- list(
      Feature = feature_names,
      Position = seq_len(max.length),
      Sample = sample_names
    )
  }
  return(out)
}

#' @rdname sequenceEncoder
#' @aliases onehotEncoder
#' @export
onehotEncoder <- function(..., mode = "onehot") {
  sequenceEncoder(..., mode = "onehot")
}

#' @rdname sequenceEncoder
#' @aliases propertyEncoder
#' @export
propertyEncoder <- function(..., mode = "property") {
  sequenceEncoder(..., mode = "property")
}

#' @rdname sequenceEncoder
#' @aliases geometricEncoder
#' @export
geometricEncoder <- function(..., mode = "geometric") {
  sequenceEncoder(..., mode = "geometric")
}