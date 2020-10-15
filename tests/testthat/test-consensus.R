consensus_test <- tibble::tribble(
  ~label, ~method, ~rank, ~taxon,
  "A",    "XTAX",  "genus", "G1",
  "A",    "XTAX",  "species", "G1 s1",
  "A",    "YTAX",  "genus", "G2")

test_that("consensus is not assumed for lower taxa when higher taxa are inconsistent", {
  expect_equal(length(phylotax(taxa = consensus_test)$tip_taxa$method), 0)
})
