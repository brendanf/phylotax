example_tree()
test_tree <- example_tree()
test_tree$tip.label <- tolower(test_tree$tip.label)

test_that("relabel_tree works", {
  expect_equal(
    relabel_tree(example_tree(), LETTERS[1:6], letters[1:6]),
    test_tree
  )
})