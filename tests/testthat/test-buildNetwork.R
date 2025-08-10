# test script for calculateMotif.R - testcases are NOT comprehensive!


## helper: tiny sequence set
seq3 <- c("AAA", "AAB", "AAC")

test_that("plain vector input returns valid edge list", {
  res <- buildNetwork(input.sequences = seq3, threshold = 1)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("from", "to", "dist"))
  expect_true(all(res$dist <= 1))
  expect_setequal(unique(c(res$from, res$to)), c("1", "2", "3"))
})


toy <- data.frame(
  cdr3   = c("CASSLG", "CASSLA", "CATSLA", "CARALA"),
  v_call = c("TRBV12", "TRBV12", "TRBV20", "TRBV12"),
  j_call = c("TRBJ2-7", "TRBJ2-7", "TRBJ2-3", "TRBJ2-7")
)

test_that("V filtering removes cross-V edges", {
  no_v  <- buildNetwork(input.data = toy, seq_col = "cdr3", threshold = 2)
  yes_v <- buildNetwork(input.data = toy, seq_col = "cdr3",
                        threshold = 2, filter.v = TRUE)
  expect_true(nrow(yes_v) < nrow(no_v))
  expect_true(all(toy$v_call[as.integer(yes_v$from)] ==
                    toy$v_call[as.integer(yes_v$to)]))
})

test_that("sparse output is symmetric dgCMatrix", {
  A <- buildNetwork(input.data = toy, seq_col = "cdr3",
                    threshold = 2, output = "sparse")
  expect_s4_class(A, "dgCMatrix")
  expect_equal(dim(A)[1], dim(A)[2])        # symmetry
  expect_equal(dimnames(A)[[1]], as.character(seq_len(nrow(toy))))
})

test_that("binary vs distance weights differ", {
  Ab <- buildNetwork(toy, seq_col = "cdr3", threshold = 2,
                     output = "sparse", weight = "binary")
  Ad <- buildNetwork(toy, seq_col = "cdr3", threshold = 2,
                     output = "sparse", weight = "dist")
  expect_true(all(Ad@x >= Ab@x))      # binary = 1, dist >= 1
})

test_that("invalid threshold throws", {
  expect_error(buildNetwork(input.sequences = seq3, threshold = 0), "threshold")
  expect_error(buildNetwork(input.sequences = seq3, threshold = -1), "threshold")
})

test_that("ids length mismatch errors", {
  expect_error(buildNetwork(input.sequences = seq3, ids = c("a", "b")), "ids")
})
