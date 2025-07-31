# test script for inferCDR.R - testcases are NOT comprehensive!

# A mock AA sequence where CDR1 is all 'C' and CDR2 is all 'D'
mock_seq_aa_1 <- paste0(paste0(rep("A", 26), collapse = ""),  # FRs
                        paste0(rep("C", 12), collapse = ""), # CDR1 (27-38)
                        paste0(rep("A", 17), collapse = ""),  # FRs
                        paste0(rep("D", 10), collapse = ""), # CDR2 (56-65)
                        paste0(rep("A", 39), collapse = "")) # FRs

# A second mock AA sequence where CDR1 is 'H' and CDR2 is 'I'
mock_seq_aa_2 <- paste0(paste0(rep("G", 26), collapse = ""),
                        paste0(rep("H", 12), collapse = ""), # CDR1
                        paste0(rep("G", 17), collapse = ""),
                        paste0(rep("I", 10), collapse = ""), # CDR2
                        paste0(rep("G", 39), collapse = ""))

# Mock Amino Acid (aa) reference list
mock_reference_aa <- list(
  sequences = c("TRBV6-4*01" = mock_seq_aa_1,
                "TRBV5-1*01" = mock_seq_aa_2),
  misc = list(species = "human", chain = "TRB", region = "v", sequence.type = "aa")
)

# Mock Nucleotide (nt) reference list (AAA -> A, CCC -> C, GGG -> D)
mock_seq_nt_1 <- paste0(paste0(rep("AAA", 26), collapse = ""),
                        paste0(rep("CCC", 12), collapse = ""), # CDR1
                        paste0(rep("AAA", 17), collapse = ""),
                        paste0(rep("GGG", 10), collapse = ""), # CDR2
                        paste0(rep("AAA", 39), collapse = ""))

mock_reference_nt <- list(
  sequences = c("TRBV6-4*01" = mock_seq_nt_1),
  misc = list(species = "human", chain = "TRB", region = "v", sequence.type = "nt")
)

mock_input_data <- data.frame(
  barcode = c("bc1", "bc2", "bc3", "bc4"),
  v_gene = c("TRBV6-4*01",            # Standard case
             "TRBV5-1*01,TRBV5-1*02", # Multiple alleles, should match first
             "TRBV-nonexistent*01",   # Gene not in reference
             NA),                     # NA gene call
  stringsAsFactors = FALSE
)


test_that("Input validation and error handling", {
  # Test for NULL or invalid reference
  expect_error(
    inferCDR(mock_input_data, reference = NULL),
    "`reference` must be a V-gene list from getIMGT()."
  )
  
  # Test for sequence type mismatch
  expect_error(
    inferCDR(mock_input_data, reference = mock_reference_aa, sequence.type = "nt"),
    "`sequence.type` mismatch with supplied reference."
  )
  
  # UPDATED: Test for missing V-gene column with the new error message
  bad_input <- data.frame(x = 1)
  expect_error(
    inferCDR(bad_input, reference = mock_reference_aa, technology = "TenX"),
    "Cannot find a V-gene column in `input.data` for technology = 'TenX'."
  )
  expect_error(
    inferCDR(bad_input, reference = mock_reference_aa, technology = "AIRR"),
    "Cannot find a V-gene column in `input.data` for technology = 'AIRR'."
  )
})

test_that("Correct CDR extraction for AA sequences (TenX/Adaptive)", {
  result <- inferCDR(mock_input_data, reference = mock_reference_aa, 
                     technology = "TenX", verbose = FALSE)
  
  expected_cdr1 <- c("CCCCCCCCCCCC", "HHHHHHHHHHHH", NA, NA)
  expected_cdr2 <- c("DDDDDDDDDD", "IIIIIIIIII", NA, NA)
  
  expect_true("CDR1_IMGT" %in% names(result))
  expect_true("CDR2_IMGT" %in% names(result))
  expect_equal(result$CDR1_IMGT, expected_cdr1)
  expect_equal(result$CDR2_IMGT, expected_cdr2)
  expect_equal(nrow(result), nrow(mock_input_data))
})

test_that("Correct CDR extraction for NT sequences", {
  input_nt <- data.frame(v_gene = "TRBV6-4*01")
  result <- inferCDR(input_nt, reference = mock_reference_nt, 
                     sequence.type = "nt", technology = "TenX", verbose = FALSE)
  
  expected_cdr1_nt <- paste0(rep("CCC", 12), collapse = "")
  expected_cdr2_nt <- paste0(rep("GGG", 10), collapse = "")
  
  expect_equal(result$CDR1_IMGT, expected_cdr1_nt)
  expect_equal(result$CDR2_IMGT, expected_cdr2_nt)
})



test_that("V-gene column selection logic is correct", {
  # 1. `v_IMGT` has highest priority
  input_v_imgt <- data.frame(v_IMGT = "TRBV6-4*01", v_call = "TRBV5-1*01")
  # `v_IMGT` should be used, ignoring `v_call` and `technology`
  result_imgt <- inferCDR(input_v_imgt, reference = mock_reference_aa, technology = "AIRR", verbose = FALSE)
  expect_equal(result_imgt$CDR1_IMGT, "CCCCCCCCCCCC") # From TRBV6-4*01
  
  # 2. TenX/Adaptive: `v_gene` is preferred over `vGeneName`
  input_tenx_pref <- data.frame(v_gene = "TRBV6-4*01", vGeneName = "TRBV5-1*01")
  result_tenx <- inferCDR(input_tenx_pref, reference = mock_reference_aa, technology = "TenX", verbose = FALSE)
  expect_equal(result_tenx$CDR1_IMGT, "CCCCCCCCCCCC") # From TRBV6-4*01
  
  # 3. TenX/Adaptive: `vGeneName` is used as a fallback
  input_tenx_fallback <- data.frame(vGeneName = "TRBV5-1*01")
  result_tenx_fallback <- inferCDR(input_tenx_fallback, reference = mock_reference_aa, technology = "Adaptive", verbose = FALSE)
  expect_equal(result_tenx_fallback$CDR1_IMGT, "HHHHHHHHHHHH") # From TRBV5-1*01
  
  # 4. AIRR: `v_call` is used
  input_airr <- data.frame(v_call = "TRBV5-1*01", v_gene = "TRBV6-4*01")
  result_airr <- inferCDR(input_airr, reference = mock_reference_aa, technology = "AIRR", verbose = FALSE)
  expect_equal(result_airr$CDR1_IMGT, "HHHHHHHHHHHH") # From TRBV5-1*01
})


test_that("Handles missing V-genes and NA gracefully", {
  result <- inferCDR(mock_input_data, reference = mock_reference_aa, technology="TenX", verbose = FALSE)
  
  # Row 3 has a V-gene not in the reference
  expect_true(is.na(result$CDR1_IMGT[3]))
  expect_true(is.na(result$CDR2_IMGT[3]))
  
  # Row 4 has an NA V-gene call
  expect_true(is.na(result$CDR1_IMGT[4]))
  expect_true(is.na(result$CDR2_IMGT[4]))
})


