# Package-level constants.
#
# This file is named `aaa-constants.R` so it collates first. The built-in
# property scales in `calculateProperty.R` reference `amino.acids` at load
# time (in their `dimnames`), and the package has no `Collate:` field, so R
# sources files in alphabetical order -- `calculateProperty.R` before
# `utils.R`. Keeping the definition here lets there be exactly one copy.

#' Standard 20 amino acids
#'
#' Vector of one-letter codes for the 20 standard amino acids.
#' @export
amino.acids <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
