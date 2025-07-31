# test script for calculateEntropy.R - testcases are NOT comprehensive!

toy <- c("AA", "AB")          
ref <- list(
  shannon      = log(2),                      # ln 2 ≃ 0.693147
  inv.simpson  = 2,                           # 1 / (0.5² + 0.5²)
  gini.simpson = 0.5,                         # 1 − 1/2
  norm.entropy = 1,                           # H / ln S  (S = 2)
  pielou       = 1,                           # identical to norm_entropy
  hill0        = 2,                           # richness
  hill1        = exp(log(2)),                # e^{H} = 2
  hill2        = 2                            # inverse-Simpson
)

metrics <- names(ref)

test_that("output vector has correct length and names", {
  out <- calculateEntropy(toy, method = "shannon")
  expect_length(out, 2)
  expect_named(out, c("Pos1", "Pos2"))
})

test_that("each built-in metric matches analytic reference", {
  for (m in metrics) {
    val <- calculateEntropy(toy, method = m)
    expect_equal(val[["Pos1"]], 0)                 # invariant column
    expect_equal(val[["Pos2"]], ref[[m]], tolerance = 1e-10,
                 info = paste("metric =", m))
  }
})

test_that("hill metrics are self-consistent", {
  h0 <- calculateEntropy(toy, method = "hill0")[["Pos2"]]
  h1 <- calculateEntropy(toy, method = "hill1")[["Pos2"]]
  h2 <- calculateEntropy(toy, method = "hill2")[["Pos2"]]
  expect_equal(h0, 2)                  # richness
  expect_equal(h1, 2)                  # exp(H)
  expect_equal(h2, 2)                  # inv-Simpson
})


test_that("explicit max.length pads sequences correctly", {
  out <- calculateEntropy(toy, 
                          max.length = 4, 
                          method = "shannon")
  expect_length(out, 4)
  # padded columns (no real AA) should return 0
  expect_true(all(out[3:4] == 0))
})

test_that("alternative padding.symbol is honoured", {
  out.dot  <- calculateEntropy(toy, 
                               padding.symbol = ".", 
                               method = "shannon")
  out.star <- calculateEntropy(toy, 
                               padding.symbol = "*", 
                               method = "shannon")
  expect_equal(out.dot, out.star)    
})

test_that("invalid inputs raise informative errors", {
  expect_error(calculateEntropy(1:5),           "character")   # non-char input
  expect_error(calculateEntropy(toy, 
                                method = "banana"),
               regexp =  "should be one of")
  expect_error(calculateEntropy(toy, 
                                padding.symbol = "XX"),
               "single character")
})
