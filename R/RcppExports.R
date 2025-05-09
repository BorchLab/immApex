# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

symmetric_deletion_lookup_cpp <- function(sequences, threshold) {
    .Call(`_immApex_symmetric_deletion_lookup_cpp`, sequences, threshold)
}

edit_distance_threshold <- function(a, b, threshold) {
    .Call(`_immApex_edit_distance_threshold`, a, b, threshold)
}

post_filter_candidates_seq <- function(candidatePairs, sequences, vGenes, jGenes, threshold, filterV, filterJ) {
    .Call(`_immApex_post_filter_candidates_seq`, candidatePairs, sequences, vGenes, jGenes, threshold, filterV, filterJ)
}

