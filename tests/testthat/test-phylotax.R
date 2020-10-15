phylotax_out <- phylotax(tree = example_tree(), taxa = example_taxa())

test_that("phylotax does not revert", {
  expect_known_value(phylotax_out$rejected, "rejected.rds", update = FALSE)
  expect_known_value(phylotax_out$retained, "retained.rds", update = FALSE)
  expect_known_value(phylotax_out$tip_taxa, "tip_taxa.rds", update = FALSE)
  expect_known_value(phylotax_out$node_taxa, "node_taxa.rds", update = FALSE)
})
