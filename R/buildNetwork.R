#' Build Edit Distance Network
#'
#' @param input.data `data.frame`/`tibble` with sequence & metadata  
#' (optional - omit if you supply `sequences` directly).
#' @param input.sequences Character vector of sequences **or** column name
#' inside `input.data`. Ignored when `NULL` and `seq_col` is non-`NULL`.
#' @param seq_col,v_col,j_col Column names to use when `input.data` is given. 
#' By default the function looks for common AIRR names (`junction_aa`, 
#' `cdr3`, `v_call`, `j_call`).
#' @param threshold >= 1 for absolute distance **or** 0 < x <= 1 for relative.
#' @param filter.v,filter.j Logical; require identical V/J when `TRUE`.
#' @param ids Optional character labels; recycled from row-names if missing.
#' @param output `"edges"` (default) or `"sparse"` - return an edge-list
#' `data.frame` **or** a symmetric `Matrix::dgCMatrix` adjacency matrix.
#' @param weight `"dist"` (store the edit distance) **or** `"binary"`
#' (all edges get weight 1). Ignored when `output = "edges"`.
#' 
#' @examples
#' data(immapex_example.data)
#' 
#' # Build Edge List
#' edges <- buildNetwork(input.data = immapex_example.data[["AIRR"]],
#'                       seq_col    = "junction_aa",
#'                       threshold  = 0.9,     
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
                         filter.v          = FALSE,
                         filter.j          = FALSE,
                         ids               = NULL,
                         output            = c("edges", "sparse"),
                         weight            = c("dist", "binary")) {
  
  output <- match.arg(output)
  weight <- match.arg(weight)
  
  ## 1. Decide where sequences come from 
  if (is.null(input.data)) {
    # user gave a bare vector
    if (is.null(input.sequences))
      stop("Provide either `input.data` *or* a `sequences` vector.")
    seq_vec <- as.character(input.sequences)
    n       <- length(seq_vec)
    v_vec <- j_vec <- NULL         
  } else {
    # user gave a data.frame / tibble
    if (!is.data.frame(input.data))
      stop("`input.data` must be a data.frame / tibble.")
    
    dat <- input.data
    # column auto-detection helpers
    guess_column <- function(x, choices)
      choices[choices %in% names(x)][1]
    
    if (is.null(seq_col))
      seq_col <- guess_column(dat, c("junction_aa", "cdr3", "sequence", "seq"))
    if (is.null(seq_col) || !seq_col %in% names(dat))
      stop("Could not find a sequence column.  Please supply `seq_col`.")
    
    seq_vec <- as.character(dat[[seq_col]])
    n       <- length(seq_vec)
    
    ## 2. V / J columns 
    if (filter.v || !is.null(v_col)) {
      if (is.null(v_col))
        v_col <- guess_column(dat, c("v_call", "v_gene", "v"))
      if (is.null(v_col) || !v_col %in% names(dat))
        stop("`v_col` not found (needed for V filtering).")
      v_vec <- as.character(dat[[v_col]])
    } else v_vec <- NULL
    
    if (filter.j || !is.null(j_col)) {
      if (is.null(j_col))
        j_col <- guess_column(dat, c("j_call", "j_gene", "j"))
      if (is.null(j_col) || !j_col %in% names(dat))
        stop("`j_col` not found (needed for J filtering).")
      j_vec <- as.character(dat[[j_col]])
    } else j_vec <- NULL
    
    ## 3. ids 
    if (is.null(ids))
      ids <- rownames(dat) %||% paste0("cell", seq_len(n))
  }
  
  ##  4. Input sanity checks 
  if (length(threshold) != 1 || !is.numeric(threshold) || threshold <= 0)
    stop("`threshold` must be > 0 (integer or 0-1).")
  if (threshold <= 1 && threshold <= 0)
    stop("Relative threshold must be 0 < x <= 1.")
  if (filter.v && is.null(v_vec))
    stop("`filter.v = TRUE` requires V gene information.")
  if (filter.j && is.null(j_vec))
    stop("`filter.j = TRUE` requires J gene information.")
  
  if (!is.null(ids) && length(ids) != n)
    stop("`ids` must have the same length as the sequence vector.")
  
  ##  5. Call the C++ engine 
  edge_df  <- fast_edge_list(
    seqs    = seq_vec,
    thresh  = threshold,
    v_gene  = v_vec,
    j_gene  = j_vec,
    match_v = filter.v,
    match_j = filter.j,
    ids     = ids)
  
  if (output == "edges")
    return(edge_df)
  
  ## convert edge list to sparse adjacency 
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

`%||%` <- function(x, y) if (is.null(x)) y else x   # tiny helper