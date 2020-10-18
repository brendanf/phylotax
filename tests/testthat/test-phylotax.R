phylotax_out <- phylotax(tree = example_tree(), taxa = example_taxa())

test_that("phylotax does not revert", {
  expect_known_value(phylotax_out$rejected, "rejected.rds", update = FALSE)
  expect_known_value(phylotax_out$retained, "retained.rds", update = FALSE)
  expect_known_value(phylotax_out$tip_taxa, "tip_taxa.rds", update = FALSE)
  expect_known_value(phylotax_out$node_taxa, "node_taxa.rds", update = FALSE)
})

incertae_taxon <- tibble::tribble(
  ~label, ~rank, ~taxon,
  "A",    "order", "ord",
  "B",    "order", "ord",
  "C",    "order",  "ord",
  "C",    "family", "fam",
  "A",    "genus", "gen",
  "B",    "genus",  "gen",
  "B",    "species", "s"
)

incertae_tree <- ape::read.tree(text = "((A,B),C);")

incertae_target <- tibble::tribble(
  ~node, ~rank, ~taxon,
  4,     "order", "ord",
  3,     "family", "fam",
  5,     "genus",  "gen",
  2,     "species", "s"
)

test_that("assignments work for incertae sedis taxa", {
  expect_true(dplyr::all_equal(phylotax::phylotax(incertae_tree, incertae_taxon)$node_taxa,
                               incertae_target))
})

