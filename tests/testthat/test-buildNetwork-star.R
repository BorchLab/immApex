# `expand = "star"` is NOT bit-identical to clique (it changes edge
# multiplicity by design). Its contract is: identical connected-components
# membership and identical vertex set. We verify that against clique mode on the
# duplicate-heavy fixtures, which is exactly what scRepertoire's default
# `cluster.method = "components"` relies on.

test_that("star expansion preserves connected components and vertex set", {
  skip_if_not_installed("igraph")
  sc <- .bn_scenarios()
  datasets <- sc$datasets

  # Each connected component as a sorted vertex-set signature; the set of these
  # signatures is the partition, compared order-independently.
  # NOTE: the integer component IDs igraph assigns are arbitrary labels that
  # depend on edge row order, and buildNetwork()'s row order is not stable
  # between expand modes (OpenMP merge order). vapply() carries those IDs
  # through as names, so they must be dropped -- otherwise expect_identical()
  # compares the labels, not the partition, and fails whenever clique and star
  # happen to number the same components differently. Vertex names are sorted
  # with method = "radix" for the same reason .bn_canon() does: byte order, so
  # the signature does not depend on LC_COLLATE.
  component_sets <- function(edge_df) {
    if (nrow(edge_df) == 0) return(character(0))
    g <- igraph::graph_from_data_frame(edge_df[, c("from", "to")], directed = FALSE)
    memb <- igraph::membership(igraph::components(g, mode = "weak"))
    groups <- split(names(memb), as.integer(memb))
    sigs <- vapply(groups,
                   function(v) paste(sort(v, method = "radix"), collapse = "\031"),
                   character(1), USE.NAMES = FALSE)
    sort(sigs, method = "radix")
  }

  cfgs <- list(
    list(data = "dup_heavy",  args = list(dist_type = "levenshtein", threshold = 2, filter.v = TRUE)),
    list(data = "dup_heavy",  args = list(dist_type = "levenshtein", threshold = 0.15, normalize = "length", filter.v = TRUE)),
    list(data = "equal_len",  args = list(dist_type = "hamming", threshold = 2, filter.v = TRUE)),
    list(data = "dup_heavy",  args = list(dist_type = "damerau", threshold = 2, filter.v = TRUE, filter.j = TRUE)),
    list(data = "unique_mixed", args = list(dist_type = "nw", dist_mat = "BLOSUM62", threshold = 5, filter.v = TRUE))
  )

  for (cfg in cfgs) {
    d <- datasets[[cfg$data]]
    base_args <- c(cfg$args, list(input.data = d, seq_col = "junction_aa",
                                  v_col = "v_call", j_col = "j_call",
                                  ids = d$sequence_id))
    e_clique <- suppressMessages(do.call(buildNetwork, c(base_args, list(expand = "clique"))))
    e_star   <- suppressMessages(do.call(buildNetwork, c(base_args, list(expand = "star"))))

    # vertex sets identical
    v_clique <- sort(unique(c(e_clique$from, e_clique$to)), method = "radix")
    v_star   <- sort(unique(c(e_star$from, e_star$to)), method = "radix")
    expect_identical(v_star, v_clique,
                     info = paste(cfg$data, "| vertex set"))

    # connected-component partition identical
    expect_identical(component_sets(e_star), component_sets(e_clique),
                     info = paste(cfg$data, "| component partition"))

    # star is never larger than clique
    expect_lte(nrow(e_star), nrow(e_clique))
  }
})
