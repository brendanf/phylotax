check_executable <- function(exec) if (nchar(Sys.which(exec)) == 0) skip(
  paste("skipping because executable: ", exec, "not present")
)

unknown_file <- system.file("extdata/unknowns.fasta.gz", package = "phylotax")
unknowns <- Biostrings::readDNAStringSet(unknown_file)
dada_ref <- system.file("extdata/sebacinales.dada2.fasta.gz", package = "phylotax")
sintax_ref <- system.file("extdata/sebacinales.sintax.fasta.gz", package = "phylotax")

test_that("dada2 taxonomy works", {
  skip_if_not_installed("dada2")
  set.seed(1)
  dada2_result <- taxonomy(unknowns, reference = dada_ref, method = "dada2")
  checkmate::expect_list(dada2_result, len = 2)
  checkmate::expect_set_equal(names(dada2_result), c("tax", "boot"))
  dada2_taxonomy <- taxtable(dada2_result)
  checkmate::expect_tibble(dada2_taxonomy, min.rows = 50)
  checkmate::expect_subset(c("label", "taxon", "rank"), names(dada2_taxonomy))
  checkmate::expect_subset(dada2_taxonomy$label, names(unknowns))
})

test_that("sintax taxonomy works", {
  check_executable("vsearch")
  set.seed(1)
  sintax_result <- taxonomy(unknowns, reference = sintax_ref, method = "sintax")
  checkmate::expect_tibble(sintax_result, min.rows = 10)
  checkmate::expect_set_equal(names(sintax_result), c("label", "seq", "hit", "strand", "c12n"))
  checkmate::expect_subset(sintax_result$label, names(unknowns))
  sintax_taxonomy <- taxtable(sintax_result)
  checkmate::expect_tibble(sintax_taxonomy, min.rows = 50)
  checkmate::expect_subset(c("label", "taxon", "rank"), names(sintax_taxonomy))
  checkmate::expect_subset(sintax_taxonomy$label, names(unknowns))
})



test_that("idtaxa taxonomy works", {
  skip_if_not_installed("DECIPHER")
  set.seed(1)
  idtaxa_model <- train_idtaxa(sintax_ref)
  checkmate::expect_class(idtaxa_model, c("Taxa", "Train"))
  idtaxa_result <- taxonomy(unknowns, reference = idtaxa_model, method = "idtaxa")
  checkmate::expect_class(idtaxa_result, "Taxa")
  expect_length(idtaxa_result, length(unknowns))
  checkmate::expect_set_equal(names(idtaxa_result), names(unknowns))
  idtaxa_taxonomy <- taxtable(idtaxa_result)
  checkmate::expect_tibble(idtaxa_taxonomy, min.rows = 30)
  checkmate::expect_subset(c("label", "taxon", "rank"), names(idtaxa_taxonomy))
  checkmate::expect_subset(idtaxa_taxonomy$label, names(unknowns))
})
