#' Shannon Diversity Index (Entropy)
#'
#' Calculates Shannon’s information entropy (often denoted *H*) for a set
#' of clone or sequence counts.
#'
#' \deqn{H' \;=\; -\sum_{i = 1}^{S} p_i \ln p_i}
#'
#' where *p*<sub>*i*</sub> = *n*<sub>*i*</sub> / *N* are the relative
#' frequencies (proportions) of each of the *S* distinct categories.
#'
#' @param cnt Numeric vector of non-negative counts (one entry per clone/
#'   residue/OTU).  Zero counts are ignored.
#'
#' @return A single numeric value (>= 0).  When `cnt` contains exactly one
#'   positive entry the function returns `0`.
#'
#' @seealso [norm_entropy()], [inv_simpson()]
#' @examples
#' counts <- c(A = 12, B = 4, C = 4)
#' shannon_entropy(counts)
#' @export
shannon_entropy <- function(cnt) {
  cnt <- cnt[cnt > 0]
  p   <- cnt / sum(cnt)
  -sum(p * log(p))
}

#' Inverse Simpson Diversity
#'
#' Computes the inverse of Simpson’s concentration index, sometimes written
#' as *1/D*.  This metric emphasizes dominant categories.
#'
#' \deqn{1/D \;=\; \frac{1}{\sum_{i} p_i^{\,2}}}
#'
#' @inheritParams shannon_entropy
#' @return Numeric value >= 1.  Equals 1 when all observations belong to a
#'   single category.
#' @examples
#' inv_simpson(c(10, 5, 1))
#' @export
inv_simpson <- function(cnt) {
  cnt <- cnt[cnt > 0]
  p   <- cnt / sum(cnt)
  1 / sum(p * p)
}

#' Gini–Simpson Diversity
#'
#' Computes the complement of Simpson’s index (also called the
#' Gini–Simpson index or probability of interspecific encounter):
#'
#' \deqn{1 - \lambda = 1 - \sum_{i} p_i^{\,2}}
#'
#' @inheritParams shannon_entropy
#' @return Value in the interval [0, 1].  Higher numbers indicate greater
#'   heterogeneity.
#' @examples
#' gini_simpson(c(10, 5, 5))
#' @export
gini_simpson <- function(cnt) 1 - 1 / inv_simpson(cnt)

#' Normalised Shannon Entropy
#'
#' Shannon entropy scaled to the interval [0, 1] by its maximum possible
#' value given *S* observed categories:
#'
#' \deqn{H^* = \frac{H'}{\ln S}}
#'
#' (also known as “Shannon evenness”).
#'
#' @inheritParams shannon_entropy
#' @return Numeric value in [0, 1]; `0` when all observations are in a
#'   single category.
#' @examples
#' norm_entropy(c(40, 10, 10, 10))
#' @export
norm_entropy <- function(cnt) {
  cnt <- cnt[cnt > 0]
  S   <- length(cnt)
  if (S == 1L) 0 else shannon_entropy(cnt) / log(S)
}

#' Pielou’s Evenness
#'
#' Convenience wrapper for normalized Shannon entropy (*E* = *H* / ln *S*).
#'
#' @inheritParams shannon_entropy
#' @return Numeric evenness measure in [0, 1].
#' @examples
#' pielou_evenness(c(3, 3, 3))
#' @export
pielou_evenness <- function(cnt) norm_entropy(cnt)

#' Hill-Number Generator
#'
#' Returns a *function* that computes the Hill diversity of order *q*
#' (also called the “effective number of species”):
#'
#' \deqn{^{q}D \;=\;
#'    \left( \sum_{i} p_i^{\,q} \right)^{1/(1-q)}, \quad q \neq 1}
#'
#' For *q = 1* the formula is undefined; the limit is
#' \deqn{^{1}D = e^{H'}}.
#'
#' @param q Numeric order of diversity.  Common values:
#'   *0* (richness), *1* (exp(*H*)), *2* (inverse Simpson).
#'
#' @return A **closure**: `hill_q(q)` returns a function that takes a
#'   vector of counts and yields the corresponding ^qD.  The returned
#'   function is vectorised over its input.
#'
#' @section References:
#' Hill, M. O. (1973) *Diversity and Evenness: A Unifying Notation and its
#' Consequences.* Ecology **54** (2), 427–432.
#'
#' @examples
#' hill1 <- hill_q(1)   # q = 1
#' hill1(c(5, 1, 1, 1))
#'
#' hill2 <- hill_q(2)   # q = 2, inverse-Simpson
#' hill2(c(5, 1, 1, 1))
#' @export
hill_q <- function(q) {
  force(q)
  function(cnt) {
    cnt <- cnt[cnt > 0]
    p   <- cnt / sum(cnt)
    if (q == 1) {
      exp(-sum(p * log(p)))
    } else {
      (sum(p ^ q))^(1 / (1 - q))
    }
  }
}

#' ACE Richness Estimator
#'
#' Calculates the Abundance-based Coverage Estimator (ACE) of species richness.
#' This metric is particularly useful for datasets with a large number of rare
#' species.
#'
#' \deqn{S_{ace} = S_{abund} + \frac{S_{rare}}{C_{ace}} + \frac{F_1}{C_{ace}} \gamma^2_{ace}}
#'
#' where the classification of rare and abundant species is based on a
#' threshold of 10 individuals, *F*1 is the count of singletons,
#' *S*rare is the number of rare species, and *C*ace is the
#' sample coverage for rare species.
#'
#' @inheritParams shannon_entropy
#' @return A single numeric value representing the estimated total number of
#'   species. The estimate is constrained to be at least the number of
#'   observed species.
#'
#' @references
#' Chao, A., & Lee, S.-M. (1992). *Estimating the number of classes via sample coverage*.
#' Journal of the American Statistical Association, 87(417), 210-217.
#'
#' @examples
#' counts <- rpois(50, lambda=1.5)
#' ace_richness(counts)
#' @export
ace_richness <- function(cnt) {
  cnt <- cnt[cnt > 0]
  S_obs <- length(cnt)
  if (S_obs == 0) return(0)
  
  rare_thresh <- 10
  is_rare <- cnt <= rare_thresh
  S_rare <- sum(is_rare)
  S_abund <- sum(!is_rare)
  
  if (S_rare == 0) return(S_obs)
  
  cnt_rare <- cnt[is_rare]
  N_rare <- sum(cnt_rare)
  F1 <- sum(cnt_rare == 1)
  
  if (N_rare == 0) return(S_abund)
  
  C_ace <- 1 - F1 / N_rare
  if (C_ace == 0) return(S_obs) # Avoid division by zero
  
  sum_gamma <- sum(cnt_rare * (cnt_rare - 1))
  
  gamma_sq_num <- S_rare * sum_gamma
  gamma_sq_den <- C_ace * N_rare * (N_rare - 1)
  
  gamma_sq <- if (gamma_sq_den == 0) {
    0
  } else {
    max(0, (gamma_sq_num / gamma_sq_den) - 1)
  }
  
  est <- S_abund + (S_rare / C_ace) + (F1 / C_ace) * gamma_sq
  
  max(S_obs, est)
}

#' Gini Coefficient of Abundance Inequality
#'
#' Calculates the Gini coefficient, a measure of inequality, for a vector
#' of clone/sequence counts. It ranges from 0 (perfect equality) to nearly 1
#' (maximal inequality).
#'
#' \deqn{G = \frac{\sum_{i=1}^{S} (2i - S - 1) n_i}{S \sum_{i=1}^{S} n_i}}
#'
#' where *n*i are the counts of each of the *S* categories,
#' sorted in non-decreasing order.
#'
#' @inheritParams shannon_entropy
#' @return A numeric value in [0, 1]. Returns `0` if there is only one
#'   category.
#'
#' @seealso [gini_simpson()]
#' @examples
#' # High inequality
#' gini_coef(c(100, 1, 1, 1))
#' # Perfect equality
#' gini_coef(c(10, 10, 10, 10))
#' @export
gini_coef <- function(cnt) {
  cnt <- cnt[cnt > 0]
  S <- length(cnt)
  if (S <= 1) return(0)
  
  cnt_sorted <- sort(cnt)
  i <- 1:S
  
  num <- sum((2 * i - S - 1) * cnt_sorted)
  den <- S * sum(cnt_sorted)
  
  num / den
}

#' Dxx Dominance Index
#'
#' Calculates the minimum number of top clones/sequences (ranked by abundance)
#' that constitute a specified percentage of the total dataset. This function
#' allows the user to designate the percentage.
#'
#' @param cnt Numeric vector of non-negative counts.
#' @param pct A numeric value (0-100) for the target percentage.
#'
#' @return The smallest number of categories whose cumulative abundance is at
#'   least `pct` percent of the total abundance.
#'
#' @seealso [d50_dom()]
#' @examples
#' counts <- c(100, 50, 20, 10, 5, rep(1, 5)) 
#' dxx_dom(counts, 80)
#' @export
dxx_dom <- function(cnt, pct) {
  if (pct < 0 || pct > 100) stop("`pct` must be between 0 and 100.")
  cnt <- cnt[cnt > 0]
  if (length(cnt) == 0) return(0)
  
  cnt_sorted <- sort(cnt, decreasing = TRUE)
  cum_sum <- cumsum(cnt_sorted)
  target_sum <- sum(cnt) * pct / 100
  
  # which.max returns the index of the first TRUE
  which.max(cum_sum >= target_sum)
}

#' D50 Dominance Index
#'
#' A convenience wrapper for `dxx_dom(cnt, 50)`. Calculates the minimum
#' number of top clones required to constitute 50% of the total abundance.
#'
#' @inheritParams shannon_entropy
#' @return The smallest number of categories whose cumulative abundance is at
#'   least 50% of the total.
#' @examples
#' d50_dom(c(100, 50, 20, 10, 5, rep(1, 5)))
#' @export
d50_dom <- function(cnt) {
  dxx_dom(cnt, 50)
}

#' Chao1 Richness Estimator
#'
#' Calculates the Chao1 non-parametric estimator of species richness. 
#'
#' The bias-corrected formula is used:
#' \deqn{S_{chao1} = S_{obs} + \frac{F_1(F_1 - 1)}{2(F_2 + 1)}}
#' where *S*obs is the number of observed species, *F*1
#' is the count of singletons, and *F*2 is the count of doubletons.
#'
#' If the conditions for the formula are not met (*F*1 <= 1 or
#' *F*2 = 0), the function returns the observed richness (*S*obs).
#'
#' @inheritParams shannon_entropy
#' @return A single numeric value representing the estimated total number of
#'   species.
#'
#' @references
#' Chao, A. (1984). *Nonparametric estimation of the number of classes in a
#' population*. Scandinavian Journal of Statistics, 11(4), 265-270.
#'
#' @examples
#' # Sample with singletons and doubletons
#' counts <- c(rep(1, 10), rep(2, 5), 5, 8, 12)
#' chao1_richness(counts)
#'
#' # Sample without doubletons returns observed richness
#' chao1_richness(c(rep(1, 5), 3, 4, 5))
#' @export
chao1_richness <- function(cnt) {
  cnt <- cnt[cnt > 0]
  S_obs <- length(cnt)
  if (S_obs == 0) return(0)
  
  F1 <- sum(cnt == 1)
  F2 <- sum(cnt == 2)
  
  # Use bias-corrected Chao1 formula
  if (F1 > 1 && F2 > 0) {
    S_obs + (F1 * (F1 - 1)) / (2 * (F2 + 1))
  } else {
    # Fallback to observed richness if conditions not met
    S_obs
  }
}


## Diversity function list for switch
.div.registry <- list(
  shannon      = shannon_entropy,
  inv.simpson  = inv_simpson,
  gini.simpson = gini_simpson,
  norm.entropy = norm_entropy,
  pielou       = pielou_evenness,
  hill0        = hill_q(0),   # richness
  hill1        = hill_q(1),   # exp(H)
  hill2        = hill_q(2)    # 1/Simpson
)
