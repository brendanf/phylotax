test_that("keep_tips works with mrca=FALSE", {
  phylotax_test <- phylotax(example_tree(), example_taxa())
  phylotax_keep <- keep_tips(phylotax_test, LETTERS[1:3], mrca = FALSE)
  checkmate::expect_set_equal(phylotax_keep$tree$tip.label, LETTERS[1:3])
  checkmate::expect_subset(phylotax_keep$assigned$label, LETTERS[1:3])
  checkmate::expect_subset(phylotax_keep$missing$label, LETTERS[1:3])
  checkmate::expect_subset(phylotax_keep$retained$label, LETTERS[1:3])
  checkmate::expect_subset(phylotax_keep$rejected$label, LETTERS[1:3])
})

test_that("keep_tips works with mrca=TRUE", {
  phylotax_test <- phylotax(example_tree(), example_taxa())
  phylotax_keep <- keep_tips(phylotax_test, LETTERS[4:5], mrca = TRUE)
  checkmate::expect_set_equal(phylotax_keep$tree$tip.label, LETTERS[4:6])
  checkmate::expect_subset(phylotax_keep$assigned$label, LETTERS[4:6])
  checkmate::expect_subset(phylotax_keep$missing$label, LETTERS[4:6])
  checkmate::expect_subset(phylotax_keep$retained$label, LETTERS[4:6])
  checkmate::expect_subset(phylotax_keep$rejected$label, LETTERS[4:6])
})

test_that("keep_tips works with invert=TRUE", {
  phylotax_test <- phylotax(example_tree(), example_taxa())
  phylotax_keep <- keep_tips(phylotax_test, LETTERS[4:5], mrca = TRUE, invert = TRUE)
  checkmate::expect_set_equal(phylotax_keep$tree$tip.label, LETTERS[1:3])
  checkmate::expect_subset(phylotax_keep$assigned$label, LETTERS[1:3])
  checkmate::expect_subset(phylotax_keep$missing$label, LETTERS[1:3])
  checkmate::expect_subset(phylotax_keep$retained$label, LETTERS[1:3])
  checkmate::expect_subset(phylotax_keep$rejected$label, LETTERS[1:3])
})
