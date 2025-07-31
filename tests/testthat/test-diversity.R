# test script for diversity.R - testcases are NOT comprehensive!

## ------------------------------------------------------------------ 
## 1. Analytic reference checks                                       
## ------------------------------------------------------------------ 
counts <- c(A = 2, B = 2)   

test_that("Shannon matches ln2 for two equal categories", {
  expect_equal(shannon_entropy(counts), log(2), tolerance = 1e-12)
})

test_that("Inverse Simpson gives 2 for equal split", {
  expect_equal(inv_simpson(counts), 2)
})

test_that("Gini–Simpson complements inverse Simpson", {
  gs <- gini_simpson(counts)
  expect_equal(gs, 1 - 1 / inv_simpson(counts))
  expect_equal(gs, 0.5)
})

test_that("Normalised entropy and Pielou evenness equal 1 under even split", {
  expect_equal(norm_entropy(counts), 1)
  expect_equal(pielou_evenness(counts), 1)
})

## ------------------------------------------------------------------ 
## 2. Edge-case handling                                              
## ------------------------------------------------------------------ 
single <- c(Solo = 7)

test_that("All metrics return 0/1 for single category", {
  expect_equal(shannon_entropy(single), 0)
  expect_equal(inv_simpson(single),     1)
  expect_equal(gini_simpson(single),    0)
  expect_equal(norm_entropy(single),    0)
  expect_equal(pielou_evenness(single), 0)
  expect_equal(hill_q(0)(single),       1)  # richness
})

zeros <- c(A = 3, B = 0, C = 2)

test_that("Zero counts are ignored", {
  # effective counts = c(3,2)
  expect_equal(shannon_entropy(zeros),
               shannon_entropy(c(3,2)))
})

neg <- c(A = 5, B = -1, C = 1)

test_that("Negative counts are silently discarded (cnt > 0 filter)", {
  expect_equal(inv_simpson(neg), inv_simpson(c(5,1)))
})

## ------------------------------------------------------------------ 
## 3. Hill-number relationships                                       
## ------------------------------------------------------------------ 
hill_vec <- c(10, 5, 1)   # uneven distribution

H  <- shannon_entropy(hill_vec)
p2 <- inv_simpson(hill_vec)

test_that("Hill numbers reproduce known equivalents", {
  expect_equal(hill_q(0)(hill_vec), length(hill_vec[hill_vec > 0]))   
  expect_equal(hill_q(1)(hill_vec), exp(H),   tolerance = 1e-12)
  expect_equal(hill_q(2)(hill_vec), p2,       tolerance = 1e-12)
})

## ------------------------------------------------------------------ 
## 4. Vectorization & length-one inputs                               
## ------------------------------------------------------------------ 
test_that("Functions accept length-one numeric and give scalar output", {
  expect_equal(shannon_entropy(5), 0)
  expect_equal(norm_entropy(5),    0)
})

## ------------------------------------------------------------------ 
## 5. Large random vectors                        
## ------------------------------------------------------------------ 
set.seed(42)
large_cnt <- sample(1:100, 100, TRUE)

test_that("All metrics return finite positive numbers for large counts", {
  expect_true(is.finite(shannon_entropy(large_cnt)))
  expect_true(inv_simpson(large_cnt)  > 0)
  expect_true(gini_simpson(large_cnt) >= 0 && gini_simpson(large_cnt) <= 1)
  expect_true(norm_entropy(large_cnt) >= 0 && norm_entropy(large_cnt) <= 1)
  expect_true(pielou_evenness(large_cnt) >= 0 && pielou_evenness(large_cnt) <= 1)
})
