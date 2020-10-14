test_that("check_ranks accepts example", {
  expect_equal(
    check_ranks(example_taxa()),
    check_ranks(example_taxa(), ranks = NULL)
  )
  expect_equal(
    check_ranks(example_taxa()),
    check_ranks(example_taxa(), ranks = default_ranks)
  )
})

test_that("check_ranks makes the ranks an ordered factor", {
  expect_true(is.ordered(check_ranks(example_taxa())$rank))
})

test_that("check_ranks fails when the ranks don't match", {
  expect_error(check_ranks(example_taxa(), ranks = toupper(default_ranks)))
})