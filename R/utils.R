`%||%` <- function(x, y) if (is.null(x)) y else x   # tiny helper

#' Standard 20 amino acids
#'
#' Vector of one-letter codes for the 20 standard amino acids.
#' @export
amino.acids <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

"%!in%" <- Negate("%in%")

.summary.function <- function(x) {
  x <- switch(x,
              "mean"   = base::mean,
              "median" = stats::median,
              "sum"    = base::sum,
              "min"    = base::min,
              "max"    = base::max,
              stop("Unknown summary.fun keyword: ", x))
}

.scale.counts <- function(tab, mode, denominator = NULL) {
  # If a denominator is not provided, calculate it based on the input 'tab'
  if (is.null(denominator)) {
    denominator <- sum(tab)
  }
  switch(mode,
         count      = tab,
         proportion = {
           if (sum(denominator) == 0) {
             tab[] <- 0
             tab
           } else {
             tab / denominator
           }
         },
         percent    = {
           if (sum(denominator) == 0) {
             tab[] <- 0
             tab
           } else {
             100 * tab / denominator
           }
         })
}


#' @importFrom utils head
.check.sequences <- function(sequences, dict) {
  pat <- sprintf("[^%s]", paste(dict, collapse = ""))
  any(grepl(pat, head(sequences, 20L)))  
}

#Add additional sequence padding to max length
.padded.strings <- function(strings, max.length, pad = ".", collapse = TRUE) {
  # 1. Truncate all strings to be AT MOST max.length.
  truncated_strings <- substr(strings, 1, max.length)
  
  # 2. Calculate the needed padding for the  strings.
  needed_padding <- max.length - nchar(truncated_strings)
  
  # 3. Add the padding to the right of each string.
  final_strings <- paste0(truncated_strings, strrep(pad, needed_padding))
  
  if (collapse) {
    return(final_strings)
  } else {
    return(strsplit(final_strings, ""))
  }
}

.substring.extractor <- function(strings, k) {
  vapply(strings, function(s) {
    n <- nchar(s)
    if (n < k) return(NA_character_)
    starts <- seq_len(n - k + 1L)
    substring(s, starts, starts + k - 1L)
  }, FUN.VALUE = character(1), USE.NAMES = FALSE)
}


.min.max.normalize <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0, length(x)))   # avoids div-by-zero
  (x - rng[1L]) / diff(rng)
}


#' @importFrom matrixStats colMedians colMeans2 colSums2 colVars colMads
.get.stat.function <- function(method) {
  if (!requireNamespace("matrixStats", quietly = TRUE))
    stop("matrixStats is required for summary statistics.", call. = FALSE)
  switch(method,
         median = matrixStats::colMedians,
         mean   = matrixStats::colMeans2,
         sum    = matrixStats::colSums2,
         vars   = matrixStats::colVars,
         mads   = matrixStats::colMads,
         stop("Invalid `method`"))
}

.is_seurat_object <- function(obj) inherits(obj, "Seurat")
.is_se_object <- function(obj) inherits(obj, "SummarizedExperiment")
.is_seurat_or_se_object <- function(obj) {
  .is_seurat_object(obj) || .is_se_object(obj)
}

.get.genes.updated <- function(data, tech, region) {
  if (.is_seurat_or_se_object(data)) return(region)
  
  tech <- tech %||% "other"
  key  <- list(TenX = paste0(region, "_gene"),
               Adaptive = paste0(region, "_gene"),
               AIRR = paste0(region, "_call"))[[tech]]
  
  if (is.null(key) || !(key %in% colnames(data)))
    key <- paste0(region, "GeneName")
  
  key
}

.array.dimnamer <- function(arr) {
  comb <- expand.grid(dimnames(arr)[2:3], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  comb <- comb[order(as.integer(comb[,1])), ]    # numeric-aware sort
  paste(comb[,1], comb[,2], sep = "_")
}

.array.reshape <- function(x, dim, order = c("C", "F")) {
  
  order <- match.arg(order)
  dim   <- as.integer(dim)
  
  # Checks 
  if (!is.atomic(x))
    stop("`x` must be an atomic vector, matrix, or array.", call. = FALSE)
  if (any(dim <= 0))
    stop("All elements of `dim` must be positive.", call. = FALSE)
  if (length(x) != prod(dim))
    stop(sprintf("Input of length %d cannot be reshaped to %s (product = %d).",
                 length(x), paste(dim, collapse = "x"), prod(dim)),
         call. = FALSE)
  
  # row-major flatten 
  row_flatten <- function(a) {
    if (is.null(dim(a))) {
      as.vector(a)                                 # already 1D
    } else {
      as.vector(aperm(a, rev(seq_along(dim(a)))))  # reverse dims then flatten
    }
  }
  
  #  build result 
  if (order == "F") {
    array(as.vector(x), dim = dim)
    
  } else {  
    vals_C <- row_flatten(x)
    tmp    <- array(vals_C, dim = rev(dim))
    aperm(tmp, rev(seq_along(dim)))
  }
}

.get_v_column <- function(input.data, technology = c("TenX", "Adaptive", "AIRR")) {
  stopifnot(is.data.frame(input.data))
  technology <- match.arg(technology)
  
  if ("v_IMGT" %in% names(input.data)) {
    v.col <- "v_IMGT"
  } else {
    v.col <- switch(technology,
                    TenX     = if ("v_gene" %in% names(input.data)) "v_gene" else "vGeneName",
                    Adaptive = if ("v_gene" %in% names(input.data)) "v_gene" else "vGeneName",
                    AIRR     = "v_call")
  }
  
  if (!v.col %in% names(input.data)) {
    stop("Cannot find a V-gene column in `input.data` for technology = '", technology, "'.")
  }
  
  v.col
}

.match.gene <- function(x,
                       table,
                       suffix_rx = "\\*.*$") {
  
  ## 1. first try an exact 1-to-1 match -------------------------------
  out <- match(x, table, nomatch = NA_integer_)
  
  ## 2. second pass on the still-unmatched entries --------------------
  no_hit <- is.na(out) & !is.na(x)
  if (any(no_hit)) {
    
    # strip allele suffixes from both sets once
    tbl_core <- sub(suffix_rx, "", table)
    x_core   <- sub(suffix_rx, "", x[no_hit])
    
    # attempt a core-gene match
    out[no_hit] <- match(x_core, tbl_core, nomatch = NA_integer_)
  }
  
  out
}


.transform.apply <- function(x,
                             method = c("none", "sqrt", "log1p",
                                        "zscore", "minmax")) {
  
  stopifnot(is.matrix(x) || is.numeric(x))
  method <- match.arg(method)
  
  out <- switch(method,
                none   = x,
                
                sqrt   = sqrt(pmax(x, 0)),
                
                log1p  = log1p(pmax(x, 0)),
                
                zscore = {
                  center <- rowMeans(x, na.rm = TRUE)
                  scale  <- sqrt(pmax(rowMeans(x^2, na.rm = TRUE) - center^2,
                                      .Machine$double.eps))
                  sweep(sweep(x, 1L, center), 1L, scale, "/")
                },
                
                minmax = {
                  a <- apply(x, 1L, min,  na.rm = TRUE)
                  b <- apply(x, 1L, max,  na.rm = TRUE)
                  denom <- pmax(b - a, .Machine$double.eps)
                  sweep(sweep(x, 1L, a), 1L, denom, "/")
                }
  )
  
  # keep row/col names intact
  dimnames(out) <- dimnames(x)
  out
}
