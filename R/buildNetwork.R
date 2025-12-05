#' Build Edit Distance Network 
#'
#' Build a sequence similarity network using various distance metrics and
#' normalization options. Supports Levenshtein, Hamming, Damerau-Levenshtein,
#' Needleman-Wunsch, and Smith-Waterman distances.
#'
#' @param input.data `data.frame`/`tibble` with sequence & metadata  
#' (optional - omit if you supply `sequences` directly).
#' @param input.sequences Character vector of sequences **or** column name
#' inside `input.data`. Ignored when `NULL` and `seq_col` is non-`NULL`.
#' @param seq_col,v_col,j_col Column names to use when `input.data` is given. 
#' By default the function looks for common AIRR names (`junction_aa`, 
#' `cdr3`, `v_call`, `j_call`).
#' @param threshold >= 1 for absolute distance **or** 0 < x <= 1 for relative.
#' When using normalized distances (`normalize != "none"`), this typically 
#' should be a value between 0 and 1 (e.g., 0.9 for 10 percent dissimilarity).
#' @param filter.v Logical; require identical V when `TRUE`.
#' @param filter.j Logical; require identical J when `TRUE`.
#' @param ids Optional character labels; recycled from row-names if missing.
#' @param output `"edges"` (default) or `"sparse"` - return an edge-list
#' `data.frame` **or** a symmetric `Matrix::dgCMatrix` adjacency matrix.
#' @param weight `"dist"` (store the edit distance) **or** `"binary"`
#' (all edges get weight 1). Ignored when `output = "edges"`.
#' @param dist_type Character string specifying the distance metric to use:
#'  \itemize{
#'    \item{`"levenshtein"`}  - Standard edit distance (default, backward compatible)
#'    \item{`"hamming"`}      - Hamming distance (requires equal-length sequences)
#'    \item{`"damerau"`}      - Damerau-Levenshtein (allows transpositions)
#'    \item{`"nw"`}           - Needleman-Wunsch global alignment score
#'    \item{`"sw"`}           - Smith-Waterman local alignment score
#'  }
#' @param dist_mat Character string specifying which substitution matrix 
#' to use for alignment-based metrics (`"nw"`, `"sw"`). Options include:
#'  \itemize{
#'    \item{`"BLOSUM45"`}     - BLOSUM45 matrix (distantly related)
#'    \item{`"BLOSUM50"`}     - BLOSUM50 matrix
#'    \item{`"BLOSUM62"`}     - BLOSUM62 matrix (default, good for proteins)
#'    \item{`"BLOSUM80"`}     - BLOSUM80 matrix (closely related)
#'    \item{`"BLOSUM100"`}    - BLOSUM100 matrix (very closely related)
#'    \item{`"PAM30"`}        - PAM30 matrix (closely related sequences)
#'    \item{`"PAM40"`}        - PAM40 matrix
#'    \item{`"PAM70"`}        - PAM70 matrix
#'    \item{`"PAM120"`}       - PAM120 matrix
#'    \item{`"PAM250"`}       - PAM250 matrix (distantly related)
#'  }
#' @param normalize Character string specifying how to normalize distances:
#'  \itemize{
#'    \item{`"none"`}         - Raw distance values (default, backward compatible)
#'    \item{`"maxlen"`}       - Normalize by max(length(seq1), length(seq2))
#'    \item{`"length"`}       - Normalize by mean sequence length
#'  }
#' @param gap_open Gap opening penalty for alignment-based metrics (default: -10).
#' Only used when `metric` is "nw" or "sw".
#' @param gap_extend Gap extension penalty for alignment-based metrics (default: -1).
#' Only used when `metric` is "nw" or "sw".
#' 
#' @examples
#' data(immapex_example.data)
#' 
#' # Levenshtein distance
#' edges <- buildNetwork(input.data = immapex_example.data[["AIRR"]],
#'                       seq_col    = "junction_aa",
#'                       threshold  = 0.9,     
#'                       filter.v   = TRUE)
#'
#' # Using Hamming distance with normalization
#' edges <- buildNetwork(input.data = immapex_example.data[["AIRR"]],
#'                       seq_col    = "junction_aa",
#'                       threshold  = 0.1,
#'                       dist_type  = "hamming",
#'                       normalize  = "maxlen",
#'                       filter.v   = TRUE)
#'
#' # Using Needleman-Wunsch with BLOSUM62
#' edges <- buildNetwork(input.data = immapex_example.data[["AIRR"]],
#'                       seq_col    = "junction_aa",
#'                       threshold  = 0.2,
#'                       dist_type  = "nw",
#'                       normalize  = "maxlen",
#'                       dist_mat   = "BLOSUM62",
#'                       filter.v   = TRUE)
#'
#' # Using PAM30 for closely related sequences
#' edges <- buildNetwork(input.data = immapex_example.data[["AIRR"]],
#'                       seq_col    = "junction_aa",
#'                       threshold  = 0.15,
#'                       dist_type  = "nw",
#'                       normalize  = "maxlen",
#'                       dist_mat   = "PAM30",
#'                       filter.v   = TRUE)
#'
#' # Damerau-Levenshtein (allows transpositions)
#' edges <- buildNetwork(input.data = immapex_example.data[["AIRR"]],
#'                       seq_col    = "junction_aa",
#'                       threshold  = 2,
#'                       dist_type  = "damerau",
#'                       filter.v   = TRUE)
#'
#' @return edge-list `data.frame` **or** sparse adjacency `dgCMatrix`
#' @importFrom Matrix sparseMatrix
#' @export
buildNetwork <- function(input.data        = NULL,
                         input.sequences   = NULL,
                         seq_col           = NULL,
                         v_col             = NULL,
                         j_col             = NULL,
                         threshold         = 2,
                         dist_type         = "levenshtein",
                         dist_mat          = NULL,
                         normalize         = c("none", "length", "maxlen"), 
                         gap_open          = -10,
                         gap_extend        = -1,
                         filter.v          = FALSE,
                         filter.j          = FALSE,
                         ids               = NULL,
                         output            = c("edges", "sparse"),
                         weight            = c("dist", "binary")) {
  
  output    <- match.arg(output)
  weight    <- match.arg(weight)
  normalize <- match.arg(normalize) # Validate input
  
  ## 1. Decide where sequences come from 
  if (is.null(input.data)) {
    if (is.null(input.sequences))
      stop("Provide either `input.data` *or* a `sequences` vector.")
    seq_vec <- as.character(input.sequences)
    n       <- length(seq_vec)
    v_vec <- j_vec <- NULL         
  } else {
    if (!is.data.frame(input.data))
      stop("`input.data` must be a data.frame / tibble.")
    
    dat <- input.data
    guess_column <- function(x, choices) choices[choices %in% names(x)][1]
    
    if (is.null(seq_col))
      seq_col <- guess_column(dat, c("junction_aa", "cdr3", "sequence", "seq"))
    if (is.null(seq_col) || !seq_col %in% names(dat))
      stop("Could not find a sequence column. Please supply `seq_col`.")
    
    seq_vec <- as.character(dat[[seq_col]])
    n       <- length(seq_vec)
    
    ## 2. V / J columns 
    if (filter.v || !is.null(v_col)) {
      if (is.null(v_col)) v_col <- guess_column(dat, c("v_call", "v_gene", "v"))
      if (is.null(v_col) || !v_col %in% names(dat))
        stop("`v_col` not found (needed for V filtering).")
      v_vec <- as.character(dat[[v_col]])
    } else v_vec <- NULL
    
    if (filter.j || !is.null(j_col)) {
      if (is.null(j_col)) j_col <- guess_column(dat, c("j_call", "j_gene", "j"))
      if (is.null(j_col) || !j_col %in% names(dat))
        stop("`j_col` not found (needed for J filtering).")
      j_vec <- as.character(dat[[j_col]])
    } else j_vec <- NULL
    
    ## 3. ids 
    if (is.null(ids))
      ids <- rownames(dat) %||% paste0("cell", seq_len(n))
  }
  
  ## 4. Input sanity checks 
  if (length(threshold) != 1 || !is.numeric(threshold) || threshold < 0)
    stop("`threshold` must be >= 0.")
  
  if (filter.v && is.null(v_vec)) stop("`filter.v = TRUE` requires V gene information.")
  if (filter.j && is.null(j_vec)) stop("`filter.j = TRUE` requires J gene information.")
  if (!is.null(ids) && length(ids) != n) stop("`ids` must have the same length as the sequence vector.")
  
  # Check for Hamming length requirement 
  if (dist_type == "hamming") {
    # Efficiently check if all lengths are the same as the first one
    if (n > 0) {
      first_len <- nchar(seq_vec[1])
      if (any(nchar(seq_vec) != first_len)) {
        warning("Hamming distance requires equal-length sequences. ", 
                "Pairs of unequal length will be assigned max distance.")
      }
    }
  }
  
  # Prepare Matrix for NW/SW
  numeric_mat <- NULL
  if (dist_type %in% c("nw", "sw")) {
    # If NULL, default to BLOSUM80
    if (is.null(dist_mat)) dist_mat <- "BLOSUM80"
    
    # Convert the string name (or custom matrix) into the actual matrix object
    numeric_mat <- .fetch.matrix(dist_mat)
  }
  
  ## 5. Call the C++ engine 
  edge_df  <- fast_edge_list(
    seqs         = seq_vec,
    thresh       = threshold,
    v_gene       = v_vec,
    j_gene       = j_vec,
    match_v      = filter.v,
    match_j      = filter.j,
    ids          = ids,
    metric       = dist_type,
    subst_matrix = numeric_mat,
    gap_open     = gap_open,
    gap_extend   = gap_extend,
    normalize    = normalize  
  )
  
  if (output == "edges")
    return(edge_df)
  
  ## 6. Convert edge list to sparse adjacency 
  if (!requireNamespace("Matrix", quietly = TRUE))
    stop("Matrix package required for sparse output.")
  
  all_ids <- sort(unique(c(edge_df$from, edge_df$to)))
  idx_from <- match(edge_df$from, all_ids)
  idx_to   <- match(edge_df$to,   all_ids)
  
  x <- if (weight == "binary") rep(1L, nrow(edge_df)) else edge_df$dist
  
  A <- Matrix::sparseMatrix(
    i = c(idx_from, idx_to),          
    j = c(idx_to,   idx_from),
    x = c(x,         x),
    dims = c(length(all_ids), length(all_ids)),
    dimnames = list(all_ids, all_ids)
  )
  return(A)
}



#' Get substitution matrix from package data or custom input
#' 
#' @param matrix_name Character string or numeric matrix
#' @return Numeric matrix with amino acid row/column names
#' @keywords internal
get_substitution_matrix <- function(matrix_name) {
  # If already a matrix, validate and return
  if (is.matrix(matrix_name)) {
    if (is.null(rownames(matrix_name)) || is.null(colnames(matrix_name))) {
      stop("Custom substitution matrix must have row and column names (amino acid codes)")
    }
    return(matrix_name)
  }
  
  # Otherwise, load from package data
  if (!is.character(matrix_name) || length(matrix_name) != 1) {
    stop("`subst_matrix` must be a character string or numeric matrix")
  }
  
  # Special case: identity matrix
  if (matrix_name == "identity") {
    return(make_identity_matrix())
  }
  
  # Valid matrix names from the package
  valid_matrices <- c("BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM100",
                      "PAM30", "PAM40", "PAM70", "PAM120", "PAM250")
  
  if (!matrix_name %in% valid_matrices) {
    stop("`subst_matrix` must be one of: ", 
         paste(c(valid_matrices, "identity", "or a custom matrix"), collapse = ", "))
  }
  
  # Load the data (lazy loading should handle this)
  # The data object is named immapex_blosum.pam.matrices
  if (!exists("immapex_blosum.pam.matrices", envir = .GlobalEnv)) {
    # Try to load from package
    tryCatch({
      data("immapex_blosum.pam.matrices", envir = environment())
    }, error = function(e) {
      stop("Could not load substitution matrices from package data. ",
           "Please ensure 'immapex_blosum.pam.matrices' is available.")
    })
  }
  
  # Get the specific matrix
  if (exists("immapex_blosum.pam.matrices")) {
    matrices <- get("immapex_blosum.pam.matrices")
  } else {
    stop("Substitution matrix data not found. Please load the package properly.")
  }
  
  if (!matrix_name %in% names(matrices)) {
    stop("Matrix '", matrix_name, "' not found in immapex_blosum.pam.matrices")
  }
  
  mat <- matrices[[matrix_name]]
  
  # Validate matrix structure
  if (!is.matrix(mat)) {
    stop("Invalid matrix format for '", matrix_name, "'")
  }
  
  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    stop("Matrix '", matrix_name, "' lacks row/column names")
  }
  
  return(mat)
}

#' Create a simple identity substitution matrix
#' @keywords internal
make_identity_matrix <- function() {
  aa <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
          "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
  n <- length(aa)
  mat <- matrix(0, nrow = n, ncol = n, dimnames = list(aa, aa))
  diag(mat) <- 1
  mat[mat == 0] <- -1
  return(mat)
}
