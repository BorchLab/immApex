#' Fast Matrix Summaries
#'
#' Computes a comprehensive panel of univariate statistics for every
#' **row** *or* **column** of a numeric matrix.  It is designed for
#' lightweight feature-engineering pipelines where many summaries are
#' required up-front (e.g. before modeling).
#'
#' @param x Numeric matrix (will be coerced with \code{as.matrix()}).
#' @param margin Integer. 1 = operate row-wise; 2 = column-wise (default 2).
#' @param stats Character vector naming the statistics to return.  Any
#' combination of the following (case-insensitive):
#'  \itemize{
#'     \item \code{"min"} 
#'     \item \code{"max"}
#'     \item \code{"mean"} 
#'     \item \code{"median"}
#'     \item \code{"sd"} 
#'     \item \code{"var"}
#'     \item \code{"mad"} 
#'     \item \code{"sum"},
#'     \item \code{"iqr"}
#'     \item \code{"n"}
#'     \item \code{"na"}
#'     \item \code{"mode"}
#'     \item \code{"all"} 
#'    }
#' @param na.rm  Logical; ignore \code{NA}s when calculating statistics
#' default \code{TRUE}).
#'
#' @return A numeric matrix with one **row per object that was summarised**
#' (rows of the input when \code{margin = 1}, otherwise columns) and one 
#' **column per requested statistic**. Row-names (if present) are preserved; 
#' column names are the statistic labels.
#' @importFrom stats IQR mad median sd var
#'
#' @export
#' @examples
#' m <- matrix(rnorm(20), 4, 5,
#'             dimnames = list(paste0("g", 1:4), paste0("s", 1:5)))
#'
#' ## Column-wise summaries (default)
#' head(summaryMatrix(m))
#'
#' ## Row-wise summaries
#' head(summaryMatrix(m, margin = 1))
summaryMatrix <- function(x,
                          margin = 2,
                          stats  = "all",
                          na.rm  = TRUE) {
  
  # Preflight checks-----------------------------------------------------------
  stopifnot(margin %in% c(1, 2))
  
  x <- as.matrix(x)
  if (!is.numeric(x))
    stop("`x` must be a numeric matrix.")
  
  all_stats <- c("min","max","mean","median","sd","var",
                 "mad","sum","iqr","n","na","mode")
  if (identical(tolower(stats), "all"))
    stats <- all_stats
  
  stats <- tolower(stats)
  bad   <- setdiff(stats, all_stats)
  if (length(bad))
    stop("Unknown statistics requested: ",
         paste(bad, collapse = ", "))
  
  have_ms <- requireNamespace("matrixStats", quietly = TRUE)
  

  # Helper factories (row/col wrappers with {matrixStats} fall-back)  #
  make_row_fun <- function(ms_sym, base_fun) {
    if (have_ms) {
      f <- getExportedValue("matrixStats", ms_sym)
      force(f)                              
      function(m) f(m, na.rm = na.rm)             
    } else {
      function(m) apply(m, 1L, base_fun, na.rm = na.rm)
    }
  }
  make_col_fun <- function(ms_sym, base_fun) {
    if (have_ms) {
      f <- getExportedValue("matrixStats", ms_sym)
      force(f)
      function(m) f(m, na.rm = na.rm)             
    } else {
      function(m) apply(m, 2L, base_fun, na.rm = na.rm)
    }
  }
  
  row_min  <- make_row_fun("rowMins",    min)
  row_max  <- make_row_fun("rowMaxs",    max)
  row_mean <- make_row_fun("rowMeans2",  mean)
  row_med  <- make_row_fun("rowMedians", median)
  row_sd   <- make_row_fun("rowSds",     sd)
  row_var  <- make_row_fun("rowVars",    var)
  row_mad  <- make_row_fun("rowMads",    mad)
  row_sum  <- make_row_fun("rowSums2",   sum)
  
  # IQR – a tiny wrapper around rowQuantiles if available
  if (have_ms) {
    row_iqr <- function(m) {
      qs <- matrixStats::rowQuantiles(m, probs = c(0.25, 0.75), na.rm = na.rm)
      qs[, 2] - qs[, 1]
    }
  } else {
    row_iqr <- function(m) apply(m, 1L, IQR, na.rm = na.rm)
  }
  
  # n (non-NA) & na (NA count)
  row_n  <- function(m) rowSums(!is.na(m))
  row_na <- function(m) rowSums(is.na(m))
  
  # Mode – simple but O(n * rows); acceptable for modest matrices
  row_mode <- function(m) {
    apply(m, 1L, function(v) {
      v <- v[!is.na(v)]
      if (!length(v)) return(NA_real_)
      tb <- tabulate(match(v, sort(unique(v))))
      sort(unique(v))[which.max(tb)]
    })
  }
  
  # Column equivalents built from row versions --------------------------------
  to_col <- function(f) function(m) f(t(m))
  if (margin == 2) {
    col_min  <- to_col(row_min);  col_max  <- to_col(row_max)
    col_mean <- to_col(row_mean); col_med  <- to_col(row_med)
    col_sd   <- to_col(row_sd);   col_var  <- to_col(row_var)
    col_mad  <- to_col(row_mad);  col_sum  <- to_col(row_sum)
    col_iqr  <- to_col(row_iqr);  col_n    <- to_col(row_n)
    col_na   <- to_col(row_na);   col_mode <- to_col(row_mode)
  }
  
  # Dispatch helpers based on `margin` ----------------------------------------                              
  f <- if (margin == 1)                     # row-wise summaries
    list(min   = row_min,  max = row_max,   mean = row_mean,
         median= row_med,  sd  = row_sd,    var  = row_var,
         mad   = row_mad,  sum = row_sum,   iqr  = row_iqr,
         n     = row_n,    na  = row_na,    mode = row_mode)
  else                                      # column-wise
    list(min   = col_min,  max = col_max,   mean = col_mean,
         median= col_med,  sd  = col_sd,    var  = col_var,
         mad   = col_mad,  sum = col_sum,   iqr  = col_iqr,
         n     = col_n,    na  = col_na,    mode = col_mode)
  
  f <- f[stats]                          
  
  # Build result matrix --------------------------------------------------------
  objs <- if (margin == 1) rownames(x) else colnames(x)
  if (is.null(objs))
    objs <- seq_len(if (margin == 1) nrow(x) else ncol(x))
  
  res <- vapply(f, function(fun) fun(x), numeric(length(objs)))
  dimnames(res) <- list(obj = objs, stat = names(f))
  
  res                     
}
