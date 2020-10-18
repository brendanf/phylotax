
test_that("seqhash makes hashes", {
   expect_equal(seqhash("AGCT"), seqhash(Biostrings::DNAStringSet("AGCT")))
   expect_equivalent(seqhash("AGCT"), digest::digest("AGCT", algo = "xxhash32"))
   expect_equivalent(
      seqhash("AGCT", algo = "sha256"),
      digest::digest("AGCT", algo = "sha256")
   )
})

test_that("preserve_na works", {
   expect_equivalent(seqhash(NA_character_), NA_character_)
   expect_equivalent(
      seqhash(NA_character_, preserve_na = FALSE),
      digest::digest(NA_character_, algo = "xxhash32")
   )
   expect_equivalent(
      seqhash(c("AGCT", NA_character_)),
      c(digest::digest("AGCT", algo = "xxhash32"),
        NA_character_)
   )
})

test_that("len works", {
   expect_equivalent(nchar(seqhash("AGCT", len = 4)), 4)
})

testseqs <- c(A = "AGCT", B = "GGTCTAN")
test_that("seqhash preserves names", {
   expect_mapequal(
      seqhash(testseqs),
      seqhash(Biostrings::DNAStringSet(testseqs))
   )
})
