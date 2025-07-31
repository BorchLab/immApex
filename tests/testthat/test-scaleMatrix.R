# test script for scaleMatrix.R - testcases are NOT comprehensive!

# ------------------------------------------------------------------------
#  Utilities
# ------------------------------------------------------------------------

# helper to compute column SD 
col_sd   <- function(m) apply(m, 2L, sd)
row_sd   <- function(m) apply(m, 1L, sd)
row_mad  <- function(m) apply(m, 1L, mad)
col_mad  <- function(m) apply(m, 2L, mad)
l2_norm  <- function(v) sqrt(sum(v^2))
l1_norm  <- function(v) sum(abs(v))

set.seed(42)
base_mat <- matrix(rnorm(60, mean = 5, sd = 3), 10, 6,
                   dimnames = list(paste0("g", 1:10),
                                   paste0("s", 1:6)))

## ------------------------------------------------------------------------
##  Core structural guarantees
## ------------------------------------------------------------------------
test_that("scaleMatrix preserves dimensions and dimnames", {
  out <- scaleMatrix(base_mat, "z")        # arbitrary transform
  expect_identical(dim(out), dim(base_mat))
  expect_identical(dimnames(out), dimnames(base_mat))
})

test_that("method = 'none' returns identical matrix", {
  expect_identical(scaleMatrix(base_mat, "none"), base_mat)
})

## ------------------------------------------------------------------------
##  min–max scaling
## ------------------------------------------------------------------------
test_that("minmax rescales columns to [0,1] by default", {
  mm <- scaleMatrix(base_mat, "minmax")
  expect_true(all(abs(apply(mm, 2, min) - 0) < 1e-8))
  expect_true(all(abs(apply(mm, 2, max) - 1) < 1e-8))
})

test_that("minmax honours custom range and margin = 1", {
  rng <- c(-2, 2)
  mmr <- scaleMatrix(base_mat, "minmax", range = rng, margin = 1)
  expect_true(all(abs(apply(mmr, 1, min) - rng[1]) < 1e-8))
  expect_true(all(abs(apply(mmr, 1, max) - rng[2]) < 1e-8))
})

# ------------------------------------------------------------------------
#  Standard (z) and robust (median/MAD) scalers
# ------------------------------------------------------------------------
test_that("z-score: column means ~0 and sds ~1", {
  z <- scaleMatrix(base_mat, "z")
  expect_true(all(abs(colMeans(z)) < 1e-8))
  expect_true(all(abs(col_sd(z) - 1) < 1e-8))
})

test_that("robust_z: row medians ~0 and MADs ~1", {
  rz <- scaleMatrix(base_mat, "robust_z", margin = 1)
  expect_true(all(abs(apply(rz, 1, median)) < 1e-8))
  expect_true(all(abs(row_mad(rz) - 1) < 1e-8))
})

## ------------------------------------------------------------------------
##  unit_var  (variance only) 
## ------------------------------------------------------------------------
test_that("unit_var scales SD to 1 but keeps means", {
  uv <- scaleMatrix(base_mat, "unit_var")
  expect_true(all(abs(col_sd(uv) - 1) < 1e-8))
})

# ------------------------------------------------------------------------
#  Vector-norm scalers
# ------------------------------------------------------------------------
test_that("l2 normalisation gives column L2 norms of 1", {
  l2 <- scaleMatrix(base_mat, "l2")   # default margin = 2
  expect_true(all(abs(apply(l2, 2, l2_norm) - 1) < 1e-12))
})

test_that("l1 normalisation gives row L1 norms of 1", {
  l1 <- scaleMatrix(base_mat, "l1", margin = 1)
  expect_true(all(abs(apply(l1, 1, l1_norm) - 1) < 1e-12))
})

# ------------------------------------------------------------------------
#  Element-wise transforms
# ------------------------------------------------------------------------
test_that("sqrt transformation matches sqrt(pmax(x,0)+offset)", {
  off <- 0.01
  sq  <- scaleMatrix(base_mat, "sqrt", offset = off)
  expect_equal(sq,
               sqrt(pmax(base_mat, 0) + off),
               tolerance = 1e-12)
})

test_that("arcsinh transformation uses supplied cofactor", {
  cf <- 7
  asinh_mat <- scaleMatrix(base_mat, "arcsinh", cofactor = cf)
  expect_equal(asinh_mat,
               asinh(base_mat / cf),
               tolerance = 1e-12)
})

# ------------------------------------------------------------------------
#  Error handling
# ------------------------------------------------------------------------
test_that("invalid method or margin errors cleanly", {
  expect_error(scaleMatrix(base_mat, "nonexistent_method"))
  expect_error(scaleMatrix(base_mat, "z", margin = 3))
})
