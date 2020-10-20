#' Assign taxonomy to nucleotide sequences
#'
#' This method uses a common interface to call primary taxonomic assignment
#' algorithms (i.e., those which assign taxonomy based on a taxonomically
#' classified reference sequence database, but not based on the results of other
#' algorithms) from other R packages or external programs.
#'
#'
#'
#' @param seq (`character`` vector or something that can be coerced to
#' one, or a matrix with sequences as the column names ) Sequences to
#' assign taxonomy
#' @param reference (`character`` string giving a path to a file or the
#' result from [DECIPHER::LearnTaxa()]/[train_idtaxa()]) An appropriately
#' formatted reference database (see Details).
#' @param method (`character`` string) taxonomy assignment method.
#' Currently accepted values are "dada2", "sintax", and "idtaxa".
#' @param min_confidence (`integer`` between 0 and 100) The minimum
#' confidence to report results.
#' @param multithread (`integer` scalar) the number of processors to use
#' for assignment.
#' @param ... additional arguments for methods
#'
#' @return raw results of the taxonomy assignment, of various types depending on
#' `method`.
#'
#' @details # Return types
#'    * `taxonomy_dada2` and `taxonomy(..., method = "dada2")` return a `list`
#'      with elements "tax" and "boot", as [dada2::assignTaxonomy].
#'    * `taxonomy_sintax` and `taxonomy(..., method = "sintax")` return a
#'      `data.frame` with columns "label", "hit", "strand", and, if
#'      `min_confidence` is not given, "c12n" (short for "classification").
#'    * `taxonomy_idtaxa` and `taxonomy(..., method = "idtaxa")` gives
#'        an S4 object of classes "Taxa" and "Train".
#' Any of these can be passed to [taxtable] to get a uniform format suitable
#' for use in [phylotax].
#' @export
taxonomy <- function(seq, reference, method, min_confidence = 50, multithread = FALSE, ...) {
  # we can take a community matrix (in which case the sequences are the column
  # names) or the sequences
  if (is.matrix(seq)) {
    nreads <- Matrix::colSums(seq)
    seq <- colnames(seq)
  } else {
    nreads <- 1
    seqnames <- names(seq)
    seq <- as.character(seq)
    names(seq) <- seqnames
  }
  is.RNA <- grepl("[Uu]", seq[1])
  # seq <- if (is.RNA) Biostrings::RNAStringSet(seq) else Biostrings::DNAStringSet(seq)

  f <- switch(method,
         dada2 = taxonomy_dada2,
         sintax = taxonomy_sintax,
         idtaxa = taxonomy_idtaxa,
         stop("unknown taxonomy assignment method: ", method)
  )
  f(seq = seq, reference = reference, multithread = multithread, min_confidence = min_confidence, ...)
}


#' @param tryRC (`logical` scalar) passed on to `dada2::assignTaxonomy()`
#' @param outputBootstraps (`logical` scalar) passed on to `dada2::assignTaxonomy()`
#' @param verbose (`logical` scalar) passed on to `dada2::assignTaxonomy()`
#' @rdname taxonomy
#' @export
taxonomy_dada2 <- function(seq, reference, multithread = FALSE, min_confidence,
                           tryRC = FALSE,
                           outputBootstraps = TRUE,
                           verbose = TRUE, ...) {
  assertthat::assert_that(
    requireNamespace("dada2"),
    msg = "'dada2' package is required to assign taxonomy using DADA2."
  )
  dada2::assignTaxonomy(seqs = chartr("Uu", "Tt", seq),
                        refFasta = reference,
                        tryRC = tryRC,
                        outputBootstraps = outputBootstraps,
                        multithread = multithread,
                        verbose = verbose,
                        minBoot = min_confidence,
                        ...)
}

#' @param exec (`character` string) name of the executable to use for
#' SINTAX search.  The default is "vsearch", but "usearch" should also
#' work.  In either case, the executable should be installed and on the
#' system path.
#' @rdname taxonomy
#' @export
taxonomy_sintax <- function(seq, reference, min_confidence = NULL, multithread = FALSE, exec = NULL, ...) {
  if (is.null(exec)) {
    exec <- Sys.which(c("usearch", "vsearch"))
    exec <- exec[nchar(exec) > 0]
  }
  assertthat::assert_that(
    length(exec) > 0,
    msg = paste("External software USEARCH or VSEARCH are required to assign",
                "taxonomy with SINTAX. If they are installed and you are still",
                "getting this message, try supplying",
                "the full path to the executable with argument 'exec='")
  )
  assertthat::assert_that(
    suppressWarnings(system2(exec, "--version")) == 0,
    msg = paste("External software USEARCH or VSEARCH are required to assign",
    "taxonomy with SINTAX. Executable at", exec, "failed.  Try supplying a",
    "different executable with argument 'exec='")
  )
  args <- character()
  if (assertthat::is.string(seq) && file.exists(seq)) {
    seqfile <- seq
    seq <- seqinr::read.fasta(
      file = seqfile,
      as.string = TRUE) %>%
      {tibble::tibble(seq = unlist(seqinr::getSequence(., as.string = TRUE)),
              label = seqinr::getName(.))}
  } else {
    seqfile <- tempfile("seq", fileext = ".fasta")
    on.exit(file.remove(seqfile))
    if (methods::is(seq, "XStringSet")) {
      Biostrings::writeXStringSet(seq, seqfile)
      seq <- tibble::tibble(label = names(seq),
                            seq = as.character(seq, use.names = FALSE))
    } else if (methods::is(seq, "ShortRead")) {
      ShortRead::writeFasta(seq, seqfile)
      seq <- tibble::tibble(label = as.character(seq@id, use.names = FALSE),
                            seq = as.character(seq@sread, use.names = FALSE))
    } else if (is.character(seq)) {
      if (is.null(names(seq))) names(seq) <- seqhash(seq)
      is.RNA <- any(grepl( "[Uu]", seq))
      if (is.RNA) {
        seq <- Biostrings::RNAStringSet(chartr("Tt", "Uu", seq))
      } else {
        seq <- Biostrings::DNAStringSet(seq)
      }
      Biostrings::writeXStringSet(seq, seqfile)
      seq <- tibble::tibble(label = names(seq),
                            seq = as.character(seq, use.names = FALSE))
    }
  }
  args <- c("--sintax", seqfile)

  if (assertthat::is.string(reference) && file.exists(reference)) {
    dbfile <- reference
  } else {
    dbfile <- tempfile("db", fileext = ".fasta")
    on.exit(file.remove(dbfile))

    if (methods::is(reference, "XStringSet")) {
      Biostrings::writeXStringSet(reference, dbfile)
    } else if (methods::is(reference, "ShortRead")) {
      ShortRead::writeFasta(reference, dbfile)
    } else if (is.character(reference)) {
      if (is.null(names(reference))) stop("Taxonomy database given as character vector must be named.")
      is.RNA <- grepl("[Uu]", reference[1])
      if (is.RNA) {
        reference <- Biostrings::RNAStringSet(reference)
      } else {
        reference <- Biostrings::DNAStringSet(reference)
      }
      Biostrings::writeXStringSet(reference, dbfile)
    }
  }
  args <- c(args, "--db", dbfile)
  tablefile <- tempfile("table", fileext = ".tsv")
  args <- c(args, "--tabbedout", tablefile)
  on.exit(file.remove(tablefile))

  if (!missing(min_confidence)) {
    args <- c(args, "--sintax_cutoff", min_confidence/100)
  }
  if (!missing(multithread)) {
    if (isFALSE(multithread)) multithread <- 1
    args <- c(args, "--threads", multithread)
  }
  system2(exec, args = args)
  if (missing(min_confidence)) {
    # vsearch outputs an extra tab when it cannot place the sequence
    system2("sed", args = c("--in-place", "'s/\\t\\t\\t/\\t\\t/'", tablefile))
    readr::read_tsv(tablefile,
                    col_names = c("label", "hit", "strand"),
                    col_types = "ccc") %>%
      dplyr::left_join(seq, ., by = "label")
  } else {
    # vsearch outputs an extra tab when it cannot place the sequence
    system2("sed", args = c("--in-place", "'s/\\t\\t\\t\\t/\\t\\t\\t/'", tablefile))
    readr::read_tsv(tablefile,
                    col_names = c("label", "hit", "strand", "c12n"),
                    col_types = "cccc") %>%
      dplyr::left_join(seq, ., by = "label")
  }
}

#' @param strand (`character` string) passed on to `DECIPHER::IdTaxa()`
#' @rdname taxonomy
#' @export
taxonomy_idtaxa <- function(seq, reference, multithread = FALSE, strand = "top", min_confidence = 40, ...) {
  assertthat::assert_that(
    requireNamespace("DECIPHER"),
    msg = "'DECIPHER' package is required to assign taxonomy using IDTAXA."
  )
  if (isTRUE(multithread)) multithread <- NULL
  if (isFALSE(multithread)) multithread <- 1
  if (!methods::is(seq, "XStringSet")) {
    seq <- Biostrings::DNAStringSet(chartr("Uu", "Tt", seq))
  }
  DECIPHER::IdTaxa(test = seq,
                   trainingSet = reference,
                   strand = strand,
                   processors = multithread,
                   threshold = min_confidence,
                   ...)
}

sintax_format <- function(tax) {
  tax <- dplyr::mutate_at(
    tax,
    "c12n",
    sub,
    pattern = "(,?[dkpcofgs]:unidentified)+$",
    replacement = ""
  ) %>%
    dplyr::mutate(
      name = sub(.data$c12n, pattern = ".*,?[dkpcofgs]:", replacement = "")
    ) %>%
    tidyr::separate_rows("c12n", sep = ",") %>%
    dplyr::mutate_at("c12n", dplyr::na_if, "") %>%
    tidyr::separate("c12n", into = c("rank", "taxon"), sep = ":") %>%
    dplyr::mutate_at("taxon", dplyr::na_if, "unidentified") %>%
    dplyr::mutate_at("rank", rank_factor, abbrev = TRUE) %>%
    tidyr::spread(key = "rank", value = "taxon") %>%
    dplyr::select(-"<NA>") %>%
    dplyr::mutate_at("name", rlang::`%|%`, "unknown") %>%
    tidyr::unite(col = "Taxonomy",
                 dplyr::one_of("domain", "kingdom", "phylum", "class", "order",
                               "family", "genus", "species"),
                 sep = ";",
                 remove = FALSE) %>%
    dplyr::mutate_at("Taxonomy", sub, pattern = "(;?NA)+", replacement = "") %>%
    dplyr::mutate_at("Taxonomy", sub, pattern = "[|].*", replacement = "") %>%
    dplyr::mutate_at("Taxonomy", sub, pattern = "_sp^", replacement = "") %>%
    dplyr::mutate_at("Taxonomy", sub, pattern = "([^;]+);\\1",
                     replacement = "\\1")

  if ("species" %in% names(tax)) {
    tax <- dplyr::mutate(
      tax,
      name = paste0(.data$name, ifelse(is.na(.data$species), "_sp", ""))
    )
  } else {
    tax <- dplyr::mutate(tax, name = paste0(.data$name, "_sp"))
  }
  tax
}


#' Convert results from different taxonomic assignment algorithms to a uniform
#' format
#'
#' @param tax Results from `taxonomy()`
#' @param ... passed to methods
#'
#' @return a `tibble::tibble()` with columns: \describe{
#' \item{`label`}{sequence identifier}
#' \item{`rank`}{the rank of the assignment}
#' \item{`taxon`}{the taxon which was assigned}
#' \item{`confidence`}{the confidence of the assignment}
#' } Each query sequence will typically occupy several rows of the output, one
#' for each rank which was assigned.
#' @export
taxtable <- function(tax, ...) {
  if (methods::is(tax, "Taxa") && methods::is(tax, "Test")) {
    taxtable_idtaxa(tax, ...)
  } else if (is.list(tax) && "tax" %in% names(tax) && "boot" %in% names(tax)) {
    taxtable_dada2(tax, ...)
  } else if (is.data.frame(tax) && all(rlang::has_name(tax, c("label", "hit")))) {
    taxtable_sintax(tax, ...)
  } else {
    stop("Unknown taxonomy table format.")
  }
}


#' @param min_confidence (`integer`) The minimum confidence to include in
#' results. May be higher than the value given in `taxonomy()`,
#' but will have no effect if it is lower.
#' @rdname taxtable
#' @export
taxtable_sintax <- function(tax, min_confidence = 0, ...) {
  tidyr::separate_rows(tax, "hit", sep = ",") %>%
    dplyr::mutate_at("hit", dplyr::na_if, "") %>%
    tidyr::extract("hit", into = c("rank", "taxon", "confidence"),
                   regex = "([dkpcofgs]):([^(]+)\\(([01]\\.\\d+)\\)") %>%
    dplyr::mutate_at("taxon", dplyr::na_if, "unidentified") %>%
    dplyr::mutate_at("taxon", gsub, pattern = "Incertae",
                     replacement = "incertae") %>%
    dplyr::mutate_at("confidence", as.numeric) %>%
    dplyr::mutate_at("rank", rank_factor, abbrev = TRUE) %>%
    dplyr::select("label", "rank", "taxon", "confidence") %>%
    dplyr::filter(.data$confidence >= min_confidence, !is.na(.data$taxon))
}

#' @param seq_id (`character` vector) names for the sequences; these will be the
#' values that end up in the "`label`" column. If not given explicitly, they are
#' taken from the taxonomy results if this is possible.
#' @rdname taxtable
#' @export
taxtable_idtaxa <- function(tax, min_confidence = 0, seq_id = NULL, ...) {
  if (!missing(seq_id) && !is.null(seq_id)) names(tax) <- seq_id
  purrr::imap_dfr(tax, ~tibble::tibble(label = .y,
                                       rank = rank_factor(.x$rank),
                                       taxon = gsub(" ", "_", .x$taxon),
                                       confidence = .x$confidence / 100)) %>%
    dplyr::mutate_at("taxon", gsub, pattern = "Incertae", replacement = "incertae") %>%
    dplyr::filter(.data$rank != "rootrank", !startsWith(.data$taxon, "unclassified_"),
                  .data$confidence >= min_confidence)
}


#' @rdname taxtable
#' @export
taxtable_dada2 <- function(tax,
                           min_confidence = 0, ...) {
  taxa <- tax$tax %>%
    magrittr::set_rownames(names(rownames(.))) %>%
    tibble::as_tibble(rownames = "label") %>%
    tidyr::gather(key = "rank", value = "taxon", -1)
  conf <- (tax$boot / 100) %>%
    magrittr::set_rownames(names(rownames(.))) %>%
    tibble::as_tibble(rownames = "label") %>%
    tidyr::gather(key = "rank", value = "confidence", -1)
  dplyr::full_join(taxa, conf, by = c("label", "rank")) %>%
    dplyr::arrange(.data$label) %>%
    dplyr::mutate_at("taxon", gsub, pattern = "Incertae", replacement = "incertae") %>%
    dplyr::filter(!is.na(.data$taxon), .data$confidence >= min_confidence) %>%
    dplyr::mutate_at("rank", tolower) %>%
    dplyr::mutate_at("rank", rank_factor)
}

#' Train an IDTAXA model to a fasta file with \\[UV\\]SEARCH/SINTAX-style taxonomy
#'
#' IDTAXA needs to train a model on a taxonomic reference database before it can
#' be used to classify sequences.  This step can be time consuming, so it is
#' performed separately from the actual classification, and the result can be
#' saved for future analyses. This is a convenience function that makes it easy
#' to fit a model using `DECIPHER::LearnTaxa()` on the same reference
#' file that would be used for `taxonomy_sintax()`.
#'
#' @param fasta (`character` string) path to a fasta file containing the
#' reference sequences, with headers formatted as required for SINTAX.
#'
#' @return an object of classes `Taxa` and `Train`, as required for
#' `taxonomy_idtaxa()` or `DECIPHER::IdTaxa()`
#' @export
train_idtaxa <- function(fasta) {
  assertthat::assert_that(
    requireNamespace("DECIPHER"),
    msg = "'DECIPHER' package is required to train an IDTAXA model."
  )
  seqdata <- Biostrings::readDNAStringSet(fasta)
  taxonomy <- names(seqdata)

  #IDTAXA does not need intermediate/trailing placeholder nodes
  taxonomy <- gsub("[dkpcofgs]:unidentified[^,]+,?", "", taxonomy)
  taxonomy <- gsub("[dkpcofgs]:[^,]+incertae_sedis[^,;]*,?", "", taxonomy)

  #it does need a root node
  taxonomy <- gsub(".*tax=", "r:Root,", taxonomy)
  taxonomy <- gsub(";.*", "", taxonomy)

  taxdata <- taxa::parse_tax_data(taxonomy,
                                  class_sep = ",",
                                  class_regex = "([rdkpcofgs]):(.+)",
                                  class_key = c("taxon_rank", "taxon_name"))
  edgelist <- taxdata$edge_list
  ranklist <- c("rootrank", "domain", "kingdom", "phylum", "class", "order",
                "family", "genus", "species")
  ranks <- data.frame(Index = taxdata$taxon_indexes()[edgelist$to],
                      Name = taxdata$taxon_names()[edgelist$to],
                      Parent = taxdata$taxon_indexes()[edgelist$from],
                      Level = taxdata$n_supertaxa()[edgelist$to],
                      Rank = match.arg(taxdata$taxon_ranks()[edgelist$to],
                                       ranklist,
                                       several.ok = TRUE),
                      stringsAsFactors = FALSE)
  ranks$Parent[is.na(ranks$Parent)] <- 0

  taxonomy <- gsub("[rdkpcofgs]:", "", taxonomy)
  taxonomy <- gsub(",", ";", taxonomy)

  trainSet <- DECIPHER::LearnTaxa(train = seqdata,
                                  taxonomy = taxonomy,
                                  rank = ranks,
                                  verbose = TRUE)

}

abbrev_myco_taxa <- function() {
  tibble::tribble(
    ~pattern, ~replacement,
    "(mycota|mycetes|ales|aceae)", "",
    "incertae_sedis", "i_s",
    "Fungi\\b", "F",
    "Basidio\\b", "B",
    "Asco\\b", "A",
    "Chytridio\\b", "Chy",
    "Mucoro\\b", "Muc",
    "Glomero\\b", "Glo"
  )
}

condense_taxa <- function(taxa, abbrev = FALSE) {
  checkmate::assert_character(taxa)
  checkmate::assert(
    checkmate::check_flag(abbrev),
    checkmate::check_data_frame(abbrev)
  )
  taxa <- table(taxa)
  labels <- paste0(taxa, names(taxa), collapse = "/")
  labels <- gsub(pattern = "(.+/.+)", replacement = "<\\1>", labels)
  if (isFALSE(abbrev)) return(labels)
  if (isTRUE(abbrev)) abbrev <- abbrev_myco_taxa()
  checkmate::check_subset(c("pattern", "replacement"), colnames(abbrev))
  for (i in seq_len(nrow(abbrev))) {
    labels <- gsub(abbrev$pattern[i], abbrev$replacement[i], labels)
  }
  labels
}

collapse_non_na <- function(x) {
  x <- x[!is.na(x)]
  if (dplyr::n_distinct(x) == 1) unique(x) else NA
}


#' Make labels summarizing the taxonomy of each sequence
#'
#' @param t (`data.frame`) Taxonomy assignments to summarize. At a minimum,
#' should include columns "`label`", "`rank`", and "`taxon`".
#' @param cols (`character`) (optional) Additional columns from `t` which will
#' be included in the output. They should be identical for each value of
#' "`label`".
#' @param abbrev (`logical` or `data.frame`) If `TRUE`, use a standard set of
#' abbreviations for mycological taxon names.  Alternatively, define your own
#' abbreviations using a `data.frame` with columns "`pattern`" and
#' "`replacement`", which are passed on to [gsub()].
#'
#' @details The standard abbreviations are:
#'
#' `r gsub("\\", "&#92;&#92;", knitr::kable(abbrev_myco_taxa()))`
#'
#' @return A `data.frame` giving the old and new taxon taxon labels.
#' @export
make_taxon_labels <- function(t, cols = character(), abbrev = FALSE) {
  checkmate::assert_data_frame(t)
  checkmate::assert_subset(c("label", "rank", "taxon"), colnames(t))
  checkmate::assert_character(cols)
  checkmate::assert_subset(cols, colnames(t))
  dplyr::group_by_at(t, c("label", "rank")) %>%
    dplyr::summarize_at(
      c("taxon", cols),
      c(
        list(~ condense_taxa(., abbrev = abbrev)),
        rep(list(collapse_non_na), length(cols))
      )
    )%>%
    dplyr::group_by_at("label") %>%
    dplyr::arrange(rank) %>%
    dplyr::summarize(
      new = do.call(
        paste,
        c(
          list(.data$label[1]),
          purrr::map(cols, ~ .data[[.]]),
          list(paste0(.data$taxon, collapse = "-"))
        )
      )
    ) %>%
    dplyr::rename(old = "label")
}

#' Replace tree tip labels
#'
#' @param tree ([phylo][ape::read.tree()] object) phylogenetic tree with tip labels
#' @param old (`character`) Tip labels from `tree` which should be changed.
#' This does not have to be all the labels in the tree.
#' @param new (`character`) Replacement tip labels.  Should have the same length
#' as `old`.
#' @param quote (`logical`) If TRUE, double quotes are added around the values
#' in `new` before they are used for replacement.  This is helpful if they (may)
#' contain characters which will must be quoted in Newick format.
#' @param drop (`character`) (optional) Tip labels which should be dropped
#' before replacement.
#'
#' @return The tree with modified tip labels.
#' @export
relabel_tree <- function(tree, old, new, quote = FALSE, drop = character(0)) {
  tree <- ape::drop.tip(tree, intersect(drop, tree$tip.label))
  if (isTRUE(quote)) new <- paste0('"', new, '"')
  tree$tip.label <-
    plyr::mapvalues(tree$tip.label, old, new, warn_missing = FALSE)
  tree
}
