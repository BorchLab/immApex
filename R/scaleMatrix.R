#' Fast Matrix Scaling or Transformation
#'
#' Applies a chosen transformation to every row *or* column of a numeric
#' matrix without altering its dimensions.  Designed for lightweight
#' pre-processing pipelines ahead of machine-learning models.
#'
#' @param x Numeric matrix (coerced with \code{as.matrix()}).
#' @param method Character scalar. One of:
#'   \itemize{
#'     \item \code{"minmax"}   – rescale linearly to [\code{range}].
#'     \item \code{"z"}        – mean 0 / sd 1 (per margin).
#'     \item \code{"robust_z"} – median 0 / MAD 1 (outlier-resistant).
#'     \item \code{"unit_var"} – divide by sd (keep mean shifts).
#'     \item \code{"l2"}, \code{"l1"} – divide by Euclidean / L1 norm.
#'     \item \code{"sqrt"}     – element-wise square-root.
#'     \item \code{"log1p"}    – element-wise \code{log1p(x + offset)}.
#'     \item \code{"log2"}, \code{"log10"} – logs with small offset.
#'     \item \code{"arcsinh"}  – \code{asinh(x / cofactor)} (Flow/CyTOF).
#'     \item \code{"none"}     – return unchanged.
#'   }
#' @param margin 1 = operate row-wise, 2 = column-wise (default 2).
#' @param range Numeric length-2 vector for \code{method = "minmax"}.
#' @param offset Non-negative scalar added before logs / sqrt
#' (\emph{ignored} otherwise). Default \code{1e-8}.
#' @param cofactor Numeric > 0 for \code{method = "arcsinh"} (default 5).
#' @param na.rm Logical; drop NAs when computing summaries.
#'
#' @return Matrix of identical dimension (dimnames preserved).
#' @importFrom stats mad median sd 
#' @export
#'
#' @examples
#' m <- matrix(rnorm(20), 4, 5,
#'             dimnames = list(paste0("g", 1:4), paste0("s", 1:5)))
#' scaleMatrix(m, "minmax")
#' scaleMatrix(m, "robust_z", margin = 1)
#' scaleMatrix(m, "l2")                          
#' scaleMatrix(abs(m), "arcsinh", cofactor = 150)
scaleMatrix <- function(x,
                        method   = c("minmax", "z", "robust_z", 
                                     "unit_var", "l2", "l1",
                                     "sqrt", "log1p", "log2",
                                     "log10", "arcsinh", "none"),
                        margin   = 2,
                        range    = c(0, 1),
                        offset   = 1e-8,
                        cofactor = 5,
                        na.rm    = TRUE) {
  # Preflight checks-----------------------------------------------------------
  stopifnot(margin %in% c(1, 2),
            length(range) == 2L,
            offset >= 0, cofactor > 0)
  method <- match.arg(method)
  if (method == "none") return(x)
  
  x  <- as.matrix(x)
  dn <- dimnames(x)
  
  # Helper factories ----------------------------------------------------------
  have_ms <- requireNamespace("matrixStats", quietly = TRUE)
  
  make_row_fun <- function(ms_sym, base_fun) {
    if (have_ms) getExportedValue("matrixStats", ms_sym)
    else function(m) apply(m, 1L, base_fun, na.rm = na.rm)
  }
  make_col_fun <- function(ms_sym, base_fun) {
    if (have_ms) getExportedValue("matrixStats", ms_sym)
    else function(m) apply(m, 2L, base_fun, na.rm = na.rm)
  }
  
  row_min   <- make_row_fun("rowMins",   min)
  row_max   <- make_row_fun("rowMaxs",   max)
  row_mean  <- make_row_fun("rowMeans2", mean)
  row_sd    <- make_row_fun("rowSds",    sd)
  row_med   <- make_row_fun("rowMedians", median)
  row_mad   <- make_row_fun("rowMads",     mad)
  
  col_min   <- make_col_fun("colMins",   min)
  col_max   <- make_col_fun("colMaxs",   max)
  col_mean  <- make_col_fun("colMeans2", mean)
  col_sd    <- make_col_fun("colSds",    sd)
  col_med   <- make_col_fun("colMedians", median)
  col_mad   <- make_col_fun("colMads",     mad)
  
  # Dispatch -----------------------------------------------------------------
  eps <- .Machine$double.eps               
  
  if (method == "minmax") {
    rng_low  <- range[1L]; rng_span <- diff(range)
    
    if (margin == 1) {                 
      mn   <- row_min(x);  mx <- row_max(x)
      rng  <- pmax(mx - mn, eps)
      scl  <- sweep(x, 1L, mn, "-")     
      scl  <- sweep(scl, 1L, rng, "/")   
      x    <- rng_low + sweep(scl, 1L, rng_span, "*")
      
    } else {                      
      mn   <- col_min(x);  mx <- col_max(x)
      rng  <- pmax(mx - mn, eps)
      scl  <- sweep(x, 2L, mn, "-") 
      scl <- sweep(scl, 2L, rng, "/")
      x    <- rng_low + sweep(scl, 2L, rng_span, "*")
    }
  } else if (method == "z") {
    stats <- if (margin == 1) list(mu = row_mean(x), sd = row_sd(x))
    else               list(mu = col_mean(x), sd = col_sd(x))
    x <- sweep(sweep(x, margin, stats$mu, "-"), margin,
               pmax(stats$sd, eps), "/")
    
  } else if (method == "robust_z") {
    stats <- if (margin == 1) list(mu = row_med(x), sd = row_mad(x))
    else               list(mu = col_med(x), sd = col_mad(x))
    x <- sweep(sweep(x, margin, stats$mu, "-"), margin,
               pmax(stats$sd, eps), "/")
    
  } else if (method == "unit_var") {
    sdev <- if (margin == 1) row_sd(x) else col_sd(x)
    x <- sweep(x, margin, pmax(sdev, eps), "/")
    
  } else if (method %in% c("l2", "l1")) {
    p <- if (method == "l2") 2L else 1L
    norm_fun <- function(v) (sum(abs(v)^p, na.rm = na.rm))^(1/p)
    norms <- if (margin == 1) apply(x, 1L, norm_fun) else apply(x, 2L, norm_fun)
    x <- sweep(x, margin, pmax(norms, eps), "/")
    
  } else if (method == "sqrt") {
    x <- sqrt(pmax(x, 0) + offset)
    
  } else if (method == "log1p") {
    x <- log1p(x + offset)
    
  } else if (method %in% c("log2", "log10")) {
    log_fun <- if (method == "log2") base::log2 else base::log10
    x <- log_fun(pmax(x, 0) + offset)
    
  } else if (method == "arcsinh") {
    x <- asinh(x / cofactor)
  }
  
  dimnames(x) <- dn
  return(x)
}
