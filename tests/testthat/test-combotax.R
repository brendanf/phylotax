extra_taxa <- tibble::tribble(
  ~label, ~method, ~rank, ~taxon,
  "G",   "XTAX",  "genus", "Tax2",
  "G",   "YTAX",  "genus", "Tax1",
  "H",   "XTAX",  "genus", "Tax1",
  "H",   "YTAX",  "genus", "Tax1",
  "I",   "XTAX",  "genus", "Tax2",
  "J",   "YTAX",  "genus", "Tax1"
)

test_that("combined PHYLOTAX and LCA consensus work", {
  p <- phylotax(
    tree = example_tree(),
    taxa = dplyr::bind_rows(example_taxa(), extra_taxa),
    cons_method = "consensus",
    fallback = TRUE
  )
  checkmate::expect_tibble(p$missing, nrows = 0)
  checkmate::expect_tibble(p$retained, nrows = 9)
  checkmate::expect_tibble(p$assigned, nrows = 8)
  expect_equal(nrow(p$missing) + nrow(p$retained) + nrow(p$rejected),
               nrow(extra_taxa) + nrow(example_taxa()))
})
