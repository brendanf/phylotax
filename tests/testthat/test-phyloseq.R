check_phyloseq <- function() if (!requireNamespace("phyloseq")) skip()

o1 <- matrix( 1L, nrow = 5, ncol = 6,
              dimnames = list(letters[1:5], LETTERS[1:6]))

test_that("phyloseq conversion with tree works", {
  check_phyloseq()
  p <- phylotax(tree = example_tree(), taxa = example_taxa())
  o <- phyloseq::otu_table(o1, taxa_are_rows = FALSE)
  pseq <- phylotax_to_phyloseq(p, o)
  checkmate::expect_class(pseq, "phyloseq")
  expect_equal(phyloseq::ntaxa(pseq), 6)
})

test_that("phyloseq conversion without tree works", {
  check_phyloseq()
  p <- phylotax(taxa = example_taxa())
  o <- phyloseq::otu_table(o1, taxa_are_rows = FALSE)
  pseq <- phylotax_to_phyloseq(p, o)
  checkmate::expect_class(pseq, "phyloseq")
  expect_equal(phyloseq::ntaxa(pseq), 6)
})

o2 <- matrix(1L, nrow = 5, ncol = 7,
    dimnames = list(letters[1:5], LETTERS[1:7]))

test_that("phyloseq conversion with tree and extra taxa works", {
  check_phyloseq()
  p <- phylotax(tree = example_tree(), taxa = example_taxa())
  o <- phyloseq::otu_table(o2, taxa_are_rows = FALSE)
  pseq <- phylotax_to_phyloseq(p, o)
  checkmate::expect_class(pseq, "phyloseq")
  expect_equal(phyloseq::ntaxa(pseq), 6)
})
