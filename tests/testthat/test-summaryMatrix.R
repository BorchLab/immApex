# test script for summaryMatrix.R - testcases are NOT comprehensive!

test_that("basic shape, row/col names and default stats", {
  set.seed(1)
  m <- matrix(rnorm(12), 4, 3,
              dimnames = list(paste0("g", 1:4),
                              paste0("s", 1:3)))
  
  # column-wise (default)
  res_col <- summaryMatrix(m)
  expect_equal(dim(res_col), c(ncol(m), 12))
  expect_identical(rownames(res_col), colnames(m))
  expect_setequal(colnames(res_col),
                  c("min","max","mean","median","sd","var",
                    "mad","sum","iqr","n","na","mode"))
  
  # row-wise
  res_row <- summaryMatrix(m, margin = 1)
  expect_equal(dim(res_row), c(nrow(m), 12))
  expect_identical(rownames(res_row), rownames(m))
})

test_that("statistical values are correct for a simple matrix", {
  simple <- matrix(1:9, 3, 3,
                   dimnames = list(paste0("r", 1:3),
                                   paste0("c", 1:3)))
  
  # Row 1: 1,4,7
  r1 <- summaryMatrix(simple, margin = 1)[1,]
  expect_equal(as.numeric(r1["min"]),    min(c(1,4,7)))
  expect_equal(as.numeric(r1["max"]),    max(c(1,4,7)))
  expect_equal(as.numeric(r1["mean"]),   mean(c(1,4,7)))
  expect_equal(as.numeric(r1["median"]), median(c(1,4,7)))
  expect_equal(as.numeric(r1["sd"]),     sd(c(1,4,7)))
  expect_equal(as.numeric(r1["var"]),    var(c(1,4,7)))
  expect_equal(as.numeric(r1["sum"]),    sum(c(1,4,7)))
  expect_equal(as.numeric(r1["iqr"]),    IQR(c(1,4,7)))
  expect_equal(as.numeric(r1["n"]),      3)
  expect_equal(as.numeric(r1["na"]),     0)
  
  # Column 1: 1 2 3
  c1 <- summaryMatrix(simple, margin = 2)[1, ]
  expect_equal(as.numeric(c1["min"]),    1)
  expect_equal(as.numeric(c1["max"]),    3)
  expect_equal(as.numeric(c1["mean"]),   mean(c(1, 2, 3)))
})

test_that("sub-setting `stats` returns only requested columns", {
  m <- matrix(1:6, 2)
  res <- summaryMatrix(m, stats = c("mean", "sd"))
  expect_equal(ncol(res), 2)
  expect_identical(colnames(res), c("mean", "sd"))
})

test_that("`mode` calculation with ties chooses the smallest value", {
  mm <- matrix(c(1, 1, 2,
                 2, 2, 3,
                 4, 5, 5), 3, byrow = TRUE)
  res <- summaryMatrix(mm, margin = 1)["1", "mode"]  # first row
  expect_equal(res, 1)                               # tie 1 vs 2 → 1
})

test_that("handling of NA values respects `na.rm` and counts", {
  m_na <- matrix(c(1, NA, 3,
                   4, 5, 6), 2, byrow = TRUE)
  # default na.rm = TRUE
  res1 <- summaryMatrix(m_na, margin = 1)
  expect_equal(res1["1", "mean"], mean(c(1, NA, 3), na.rm = TRUE))
  expect_equal(res1["1", "na"],   1)
  expect_equal(res1["1", "n"],    2)
  
  # na.rm = FALSE should propagate NA to stats that rely on it
  res2 <- summaryMatrix(m_na, margin = 1, stats = "mean", na.rm = FALSE)
  expect_true(is.na(res2[1, "mean"]))
})

test_that("errors are thrown for bad input", {
  non_num <- matrix(letters[1:4], 2)
  expect_error(summaryMatrix(non_num),
               "`x` must be a numeric matrix")
  
  expect_error(summaryMatrix(matrix(1:4), stats = "unknown"),
               "Unknown statistics")
  
  expect_error(summaryMatrix(matrix(1:4), margin = 3),
               "margin %in% c\\(1, 2\\)")
})

test_that("row/column names are generated if absent", {
  m <- matrix(1:6, 2, 3)              # no dimnames
  res <- summaryMatrix(m)             # margin = 2
  expect_identical(rownames(res), as.character(seq_len(ncol(m))))
})
