o <- phyloseq::otu_table(
  matrix(
    1L,
    nrow = 5, ncol = 6,
    dimnames = list(letters[1:5], LETTERS[1:6])
  ),
  taxa_are_rows = FALSE
)

test_that("phyloseq conversion with tree works", {
  p <- phylotax(tree = example_tree(), taxa = example_taxa())
  pseq <- phylotax_to_phyloseq(p, o)
  checkmate::expect_class(pseq, "phyloseq")
  expect_equal(phyloseq::ntaxa(pseq), 6)
})

test_that("phyloseq conversion without tree works", {
  p <- phylotax(taxa = example_taxa())
  pseq <- phylotax_to_phyloseq(p, o)
  checkmate::expect_class(pseq, "phyloseq")
  expect_equal(phyloseq::ntaxa(pseq), 6)
})

o2 <- phyloseq::otu_table(
  matrix(
    1L,
    nrow = 5, ncol = 7,
    dimnames = list(letters[1:5], LETTERS[1:7])
  ),
  taxa_are_rows = FALSE
)

test_that("phyloseq conversion with tree and extra taxa works", {
  p <- phylotax(tree = example_tree(), taxa = example_taxa())
  pseq <- phylotax_to_phyloseq(p, o2)
  checkmate::expect_class(pseq, "phyloseq")
})
