#' Position-wise Amino-Acid Property Profiles
#'
#' Computes a range of summary statistics for property values of one or more AA 
#' property scales at every residue position of a set of protein (or peptide) 
#' sequences. The function is entirely vectorized: it first calls 
#' [`calculateFrequency()`] to obtain a residue-by-position **frequency** 
#' matrix *F* (each column sums to 1) and then performs a single matrix product.
#'
#' @param input.sequences Character vector of amino-acid strings.
#' @param property.set Character string (one of the supported names) 
#' Defaults to `"atchleyFactors"`, but includes: `"crucianiProperties"`, 
#' `"FASGAI"`, `"kideraFactors"`, `"MSWHIM"`, `"ProtFP"`, `"stScales"`, 
#' `"tScales"`, `"VHSE"`, `"zScales"`
#' @param summary.fun Character string (`"mean"`, `"median"`, `"sum"`,
#' `"min"`, `"max"`), **or** a function accepting a numeric vector and
#' returning length-1 numeric.  Defaults to `"mean"`.
#' @param transform Character string controlling a *post-summary*
#' transformation. One of  `"none"` (default), `"sqrt"`, `"log1p"`, 
#' `"zscore"` (row-wise), or `"minmax"` (row-wise).
#' @param max.length Integer. Pad/trim to this length
#'   (`max(nchar(sequences))` by default).
#' @param padding.symbol Single character used for right-padding. Must not be
#'   one of the 20 canonical residues.
#' @param tidy Logical; if `TRUE`, return a long-format `data.frame`
#'
#' @return A numeric matrix (*k* × *L*) **or** a tidy data.frame with columns
#' scale, position, value.
#' 
#' @examples
#' set.seed(1)
#' seqs <- c("CASSLGQGAETQYF", "CASSPGQGDYEQYF", "CASSQETQYF")
#' aa.Atchley <- calculateProperty(seqs, property.set = "atchleyFactors")
#' 
#' @export
calculateProperty <- function(input.sequences,
                              property.set  = "atchleyFactors",
                              summary.fun   = "mean",
                              transform     = "none",
                              max.length    = NULL,
                              padding.symbol = ".",
                              tidy           = FALSE) {
  
  # Preflight checks-----------------------------------------------------------
  stopifnot(is.character(input.sequences),
            padding.symbol %!in% amino.acids,
            nchar(padding.symbol) == 1L)
  
  if (is.character(summary.fun)) {
    summary.fun <- .summary.function(summary.fun)
  }
  
  if (!is.function(summary.fun))
    stop("'summary.fun' must be a one-argument numeric function ",
         "or one of the built-in keywords.")
  
  transform <- match.arg(transform,
                         c("none","sqrt","log1p","zscore","minmax"))
  
  if (is.null(max.length))
    max.length <- max(nchar(input.sequences), 1L)
  
  # 1.  Property matrix (k × 20)  
  S <- if (is.character(property.set)) {
    .aa.property.matrix(property.set)
  } else if (is.matrix(property.set) || is.data.frame(property.set)) {
    as.matrix(property.set)
  } else {
    stop("'property.set' must be a recognised name or a numeric matrix.")
  }
  
  S <- S[ , amino.acids, drop = FALSE]                 # enforce AA order
  k <- nrow(S)
  
  # 2. residue frequencies (20 × L)  
  nSeq  <- length(input.sequences)
  Freq  <- calculateFrequency(
    input.sequences     = input.sequences,
    max.length          = max.length,
    sequence.dictionary = amino.acids,
    padding.symbol      = padding.symbol,
    summary.fun         = "proportion",
    tidy                = FALSE)[amino.acids, , drop = FALSE]
  
  # 3. Summary per position 
  if (identical(summary.fun, base::mean)) {
    
    Summ <- S %*% Freq                                    # k × L  (mean)
    
  } else if (identical(summary.fun, base::sum)) {
    
    Summ <- S %*% (Freq * nSeq)                           # counts (= sum)
    
  } else {                                                # median, min, max, custom
    Counts <- Freq * nSeq                               
    Summ   <- matrix(NA_real_, k, max.length)
    
    for (j in seq_len(max.length)) {
      cnts <- Counts[, j]
      if (!any(cnts != 0)) next                               
      for (i in seq_len(k)) {
        vals <- rep(S[i, ], times = cnts)
        Summ[i, j] <- summary.fun(vals)
      }
    }
  }
  dimnames(Summ) <- list(scale    = rownames(S),
                         position = paste0("Pos.", seq_len(max.length)))

  # 4.Optional transform 
  Summ <- .transform.apply(Summ, method = transform) 
  
  # 4.Optional tidy
  if (tidy) {
    Summ <- as.data.frame.table(Summ,
                                stringsAsFactors = FALSE,
                                responseName     = "value")
    Summ$position <- as.integer(sub("^Pos\\.", "", Summ$position))
  }
  Summ
}

.aa.property.matrix <- function(key) {
  
  if (exists(key, envir = .builtin_scales, inherits = FALSE))
    return(.builtin_scales[[key]])
  
  if (requireNamespace("Peptides", quietly = TRUE)) {
    acc <- utils::getFromNamespace("AAdata", "Peptides")  
    if (key %in% names(acc)) {
      v <- do.call(rbind, acc[[key]])
      return(v)
    }
  }
  
  stop("Unknown property set: '", key, "'. ",
       "Use one of the built-ins, ",
       "or supply a custom numeric matrix.")
}

.builtin_scales <- new.env(parent = emptyenv())
amino.acids <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

.builtin_scales$atchleyFactors <- t(matrix(c(
  # A (Alanine)
  -0.591, -1.302, -0.733,  1.570, -0.146,
  # R (Arginine)
  1.538, -0.055,  1.502,  0.440,  2.897,
  # N (Asparagine)
  0.945,  0.828,  1.299, -0.169,  0.933,
  # D (Aspartic Acid)
  1.050,  0.302, -3.656, -0.259, -3.242,
  # C (Cysteine)
  -1.343,  0.465, -0.862, -1.020, -0.255,
  # Q (Glutamine)
  0.931, -0.179, -3.005, -0.503, -1.853,
  # E (Glutamic Acid)
  1.357, -1.453,  1.477,  0.113, -0.837,
  # G (Glycine)
  -0.384,  1.652,  1.330,  1.045,  2.064,
  # H (Histidine)
  0.336, -0.417, -1.673, -1.474, -0.078,
  # I (Isoleucine)
  -1.239, -0.547,  2.131,  0.393,  0.816,
  # L (Leucine)
  -1.019, -0.987, -1.505,  1.266, -0.912,
  # K (Lysine)
  1.831, -0.561,  0.533, -0.277,  1.648,
  # M (Methionine)
  -0.663, -1.524,  2.219, -1.005,  1.212,
  # F (Phenylalanine)
  -1.006, -0.590,  1.891, -0.397,  0.412,
  # P (Proline)
  0.189,  2.081, -1.628,  0.421, -1.392,
  # S (Serine)
  -0.228,  1.399, -4.760,  0.670, -2.647,
  # T (Threonine)
  -0.032,  0.326,  2.213,  0.908,  1.313,
  # W (Tryptophan)
  -0.595,  0.009,  0.672, -2.128, -0.184,
  # Y (Tyrosine)
  0.260,  0.830,  3.097, -0.838,  1.512,
  # V (Valine)
  -1.337, -0.279, -0.544,  1.242, -1.262
), nrow = 20, ncol = 5, byrow = TRUE,
dimnames = list(amino.acids, paste0("AF", 1:5))))