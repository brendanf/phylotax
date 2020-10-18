test_that("consensus is not assumed for lower taxa when higher taxa are inconsistent", {
  consensus_test <- tibble::tribble(
    ~label, ~method, ~rank, ~taxon,
    "A",    "XTAX",  "genus", "G1",
    "A",    "XTAX",  "species", "G1 s1",
    "A",    "YTAX",  "genus", "G2")
  expect_equal(length(phylotax(taxa = consensus_test)$assigned$method), 0)
  expect_equal(length(lca_consensus(taxa = consensus_test)$assigned$method), 0)
})

test_that("consensus does not overassign with implicit incertae sedis", {
  consensus_test2 <- tibble::tribble(
    ~label, ~method, ~rank, ~taxon,
    "A",    "XTAX",  "family", "F1",
    "A",    "YTAX",  "genus",  "G2"
  )
  expect_equal(length(phylotax(taxa = consensus_test2)$assigned$method), 0)
  expect_equal(length(lca_consensus(taxa = consensus_test2)$assigned$method), 0)
})
