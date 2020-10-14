test_that("phylotax does not revert", {
  expect_known_value(phylotax(tree = example_tree(), taxa = example_taxa()), "phylotax.rds", update = FALSE)
})
