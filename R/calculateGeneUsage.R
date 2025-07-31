#' Quantifcation of Gene-Locus Usage 
#'
#' Computes either the **counts**, **proportions** (default), or **percentages**
#' of one locus *or* a locus pair that are already present as columns in
#' `input.data`.  No external dependencies.
#'
#' @param input.data A data.frame whose rows are sequences / clones and whose
#' columns named in `loci` contain gene identifiers.
#' @param loci Character vector of length 1 or 2 giving the column names.
#' @param levels Optional list of length 1 or 2 with the full set of factor
#' levels to include.  Missing levels are filled with zeros. If `NULL`
#' (default) only observed levels appear.
#' @param summary.fun Character string choosing the summary statistic:
#'   * `"proportion"` (default) – each cell sums to 1 over the table.  
#'   * `"count"`      – raw counts.  
#'   * `"percent"`    – proportion × 100.
#'
#' @return Named numeric **vector** (single locus) or numeric **matrix**
#'   (paired loci). For `"proportion"` and `"percent"` results sum to 1 or 100.
#'
#' @examples
#' df <- data.frame(V = c("TRBV7-2","TRBV7-2","TRBV5-1"),
#'                  J = c("TRBJ2-3","TRBJ2-5","TRBJ2-3"))
#' calculateGeneUsage(df, "V",          summary = "count")
#' calculateGeneUsage(df, c("V","J"),   summary = "percent")
#'
#' @export
calculateGeneUsage <- function(input.data,
                               loci,
                               levels  = NULL,
                               summary.fun = c("proportion", "count", "percent")) {
  # Preflight checks-----------------------------------------------------------
  summary.fun <- match.arg(summary.fun)
  stopifnot(is.data.frame(input.data),
            is.character(loci), length(loci) %in% 1:2,
            all(loci %in% names(input.data)))

  
  ## single locus -------------------------------------------------------------
  if (length(loci) == 1L) {
    x <- input.data[[loci[1]]]
    tab <- if (is.null(levels)) table(x)
    else table(factor(x, levels = levels[[1]]))
    out <- .scale.counts(tab, summary.fun)
    return(out)          # return vector
  }
  
  ## paired locus -------------------------------------------------------------
  x <- input.data[[loci[1]]]
  y <- input.data[[loci[2]]]
  tab <- if (is.null(levels))
    table(x, y)
  else
    table(factor(x, levels = levels[[1]]),
          factor(y, levels = levels[[2]]))
  .scale.counts(tab, summary.fun)          # matrix
}
