# defined as a private function to avoid dependence on tzara
seqhash <- function(seq, algo = "xxhash32", len = NA, preserve_na = TRUE) {
  UseMethod("seqhash")
}

seqhash.character <- function(
  seq,
  algo = "xxhash32",
  len = NA,
  preserve_na = TRUE
) {
  h <- vapply(seq, digest::digest, "", algo = algo)
  if (preserve_na) h[is.na(seq)] <- NA_character_
  if (is.na(len)) {
    return(h)
  } else {
    substring(h, 1, len)
  }
}

seqhash.XStringSet <- function(
  seq,
  algo = "xxhash32",
  len = NA,
  preserve_na = TRUE
) {
  seqhash.character(as.character(seq), algo = algo, len = len, preserve_na)
}