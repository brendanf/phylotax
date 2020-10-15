default_ranks <- c("rootrank", "domain", "kingdom", "phylum",
                   "class", "order", "family", "genus",
                   "species")

rank_factor <- function(r,
                        ranks = c("rootrank", "domain", "kingdom", "phylum",
                                  "class", "order", "family", "genus",
                                  "species"),
                        abbrev = FALSE) {
  if (abbrev) {
    factor(r, levels = substr(ranks, 1, 1), labels = ranks, ordered = TRUE)
  } else {
    factor(r, levels = ranks, ordered = TRUE)
  }
}

# check that a taxonomic assignment table has ranks we can use
check_ranks <- function(taxa, ranks = NULL) {
  # there need to be some ranks
  assertthat::assert_that(assertthat::has_name(taxa, "rank"))
  # NA doesn't make sense
  assertthat::assert_that(assertthat::noNA(taxa$rank))
  # an ordered factor or integer is usable as-is
  if (is.ordered(taxa$rank) | is.numeric(taxa$rank)) return(taxa)
  # if it's a character, investigate
  if (is.character(taxa$rank)) {
    # if we weren't given ranks, then they should be from the default ranks
    if (is.null(ranks)) {
      ranks <- sort(unique(taxa$rank))
      unknown <- setdiff(ranks, default_ranks)
      if (length(unknown)) {
        stop("No explicit list of ranks given, and given rank(s) ", unknown,
             "are not known.")
      }
      ranks <- default_ranks
    }
    assertthat::assert_that(all(taxa$rank %in% ranks))
    taxa$rank <- rank_factor(taxa$rank, ranks = ranks)
    return(taxa)
  }
  # if it's a factor, investigate
  if (is.factor(taxa$rank)) {
    # if we weren't given ranks, then check if this matches the default ranks
    if (is.null(ranks)) {
      if (all(levels(taxa$rank) %in% default_ranks)) {
        taxa$rank <- ordered(taxa$rank, default_ranks)
      } else {
        futile.logger::flog.warn(
          paste("Non-default rank factor given without explicit ranks.",
                "Ordering based on existing factor levels.")
        )
        taxa$rank <- ordered(taxa$rank)
      }
    }
  }
  stop("taxa$rank should be an (ordered) factor, character, or integer.")
}

check_method <- function(taxa, method) {
  if (is.null(method)) {
    method <- character()
    names(method) <- character()
    invisible(character())
  } else if (is.character(method)) {
    if (!is.null(names(method))) {
      assertthat::assert_that(assertthat::has_name(taxa, names(method)))
      invisible(method)
    } else {
      assertthat::assert_that(length(method) == 1)
      assertthat::assert_that(assertthat::has_name(taxa, "method"))
      names(method) <- "method"
      invisible(method)
    }
  } else if (is.list(method)) {
    assertthat::assert_that(vapply(method, is.character, TRUE))
    method <- unlist(method)
    assertthat::assert_that(assertthat::has_name(taxa, names(method)))
    invisible(method)
  } else {
    stop("argument 'method' should be a single string, a named character",
         " vector, or a named list of characters.")
  }
}

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
  is.RNA <- stringr::str_detect(seq[1], "[Uu]")
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
taxonomy_sintax <- function(seq, reference, min_confidence = NULL, multithread = FALSE, exec = "vsearch", ...) {
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
      is.RNA <- any(stringr::str_detect(seq, "[Uu]"))
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
      is.RNA <- stringr::str_detect(reference[1], "[Uu]")
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

#' @param names (`character` vector) names for the sequences; these will be the
#' values that end up in the "`label`" column. If not given explicitly, they are
#' taken from the taxonomy results if this is possible.
#' @rdname taxtable
#' @export
taxtable_idtaxa <- function(tax, min_confidence = 0, names = NULL, ...) {
  if (!missing(names) && !is.null(names)) names(tax) <- names
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
taxtable_dada2 <- function(tax, names = rownames(tax$tax),
                           min_confidence = 0, ...) {
  taxa <- tax$tax %>% magrittr::set_rownames(names) %>%
    tibble::as_tibble(rownames = "label") %>%
    tidyr::gather(key = "rank", value = "taxon", -1)
  conf <- (tax$boot / 100) %>%
    magrittr::set_rownames(names) %>%
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


#### taxon_labels ####
# make labels summarizing the taxonomy of each sequence
make_taxon_labels <- function(t) {
    dplyr::group_by_at(t, c("label", "rank", "n_reads")) %>%
    dplyr::summarize(
      taxon =
        table(.data$taxon) %>%
        paste0(names(.), collapse = "/") %>%
        gsub(pattern = "(.+/.+)", replacement = "<\\1>") %>%
        gsub(pattern = "(mycota|mycetes|ales|aceae)", replacement = "") %>%
        gsub(pattern = "incertae_sedis", replacement = "i_s") %>%
        gsub(pattern = "Fungi\\b", replacement = "F") %>%
        gsub(pattern = "Basidio\\b", replacement = "B") %>%
        gsub(pattern = "Asco\\b", replacement = "A") %>%
        gsub(pattern = "Chytridio\\b", replacement = "Chy") %>%
        gsub(pattern = "Zygo\\b", replacement = "Z")
    ) %>%
    dplyr::group_by_at(c("label", "n_reads")) %>%
    dplyr::arrange(rank) %>%
    dplyr::summarize(
      tip_label = paste(
        .data$label[1],
        format(.data$n_reads[1], width = 5),
        paste0(.data$taxon, collapse = "-")
      )
    )
}

#### relabel_tree ####
# replaces tree tip labels from old with labels from new
relabel_tree <- function(tree, old, new, chimeras = character(0)) {
  tree <- ape::drop.tip(tree, intersect(chimeras, tree$tip.label))
  tree$tip.label <-
    plyr::mapvalues(tree$tip.label, old, paste0('"', new, '"'), warn_missing = FALSE)
  tree
}

count_assignments <- function(taxa) {
  dplyr::group_by_at(taxa, c("label", "rank")) %>%
    dplyr::mutate(
      n_diff = dplyr::n_distinct(.data$taxon),
      n_tot = dplyr::n()
    ) %>%
    dplyr::ungroup()
}

# If all members of the clade are assigned uniquely to one taxon, then assign
# that taxon
# Otherwise, if there exists one taxon that is at least one of the possibilities
# for all of the members, then assign that one.
clade_taxon <- function(tree, tax, node, rank) {
  # Check the direct children
  # If one of them is completely unidentified, then we don't want to assign a
  # taxon here, because there's no way to know if the unassigned child is in the
  # group or not.
  children <- phangorn::Children(tree, node)
  if (length(children) == 0) {
    tips <- tree$tip.label[node]
  } else {
    for (child in children) {
      tips <- tree$tip.label[phangorn::Descendants(tree, child, type = "tips")[[1]]]
      taxa <- dplyr::filter(tax, .data$label %in% tips, rank == !!rank)
      if (nrow(taxa) == 0) return(NA_character_)
    }
    tips <- tree$tip.label[phangorn::Descendants(tree, node, type = "tips")[[1]]]
  }
  taxa <- dplyr::filter(tax, .data$label %in% tips, .data$rank == !!rank) %>%
    dplyr::group_by_at("label")

  # If there are no relevant taxon assignments, then we can't do anything.
  if (nrow(taxa) == 0) return(NA_character_)
  # If only one thing is assigned, then assign that.
  if (dplyr::n_distinct(taxa$taxon) == 1) return(unique(taxa$taxon))
  best_taxon <- unique(taxa$taxon[taxa$n_diff == 1])
  if (length(best_taxon) > 1) return(NA_character_)
  consensus_taxon <-
    dplyr::group_map(taxa, ~ unique(.$taxon)) %>%
    purrr::reduce(intersect)
  if (length(consensus_taxon) != 1) return(NA_character_)
  consensus_taxon
}


# create an environment to hold the current state of the taxonomic assignment
# between branching recursive calls to phylotax_
new_phylotax_env <- function(tree, taxa, parent = parent.frame()) {
  rlang::child_env(
    .parent = parent,
    node_taxa = tibble::tibble(
      node = integer(),
      label = NULL,
      rank = taxa$rank[FALSE],
      taxon = character()
    ),
    tip_taxa = dplyr::filter(taxa, FALSE),
    retained = dplyr::filter(taxa, .data$label %in% tree$tip.label),
    rejected = dplyr::filter(taxa, FALSE),
    missing = dplyr::filter(taxa, !.data$label %in% tree$tip.label),
    tree = tree
  )
}


phylotax_ <- function(tree, taxa, node, ranks, method, e) {
  if (length(ranks) == 0) return()
  nodelabel <- if (!is.null(tree$node.label)) {
    tree$node.label[node - ape::Ntip(tree)]
  } else {
    as.character(node)
  }
  parents <- phangorn::Ancestors(tree, node, type = "all")
  for (r in ranks) {
    if (is.ordered(ranks)) r <- ordered(r, levels = levels(ranks))
    if (any(e$node_taxa$node %in% parents & e$node_taxa$rank == r)) next
    taxon <- clade_taxon(tree, e$retained, node, r)
    if (is.na(taxon)) {
      futile.logger::flog.debug("Could not assign a %s to node %s.", r, nodelabel)
      for (n in phangorn::Children(tree, node)) {
          phylotax_(tree, e$retained, n, ranks, method, e)
      }
      break
    } else {
      children <- phangorn::Children(tree, node)
      if (length(children) > 0) {
        futile.logger::flog.info(
          "Assigned node %s and its %d children to %s %s.",
          nodelabel, length(children), as.character(r), taxon)
      } else {
        futile.logger::flog.info("Assigned node %s to %s %s.", nodelabel,
                                 as.character(r), taxon)
      }
      ranks <- ranks[-1]
      e$node_taxa <- dplyr::bind_rows(
        e$node_taxa,
        tibble::tibble(node = node, label = nodelabel, rank = r, taxon = taxon)
      )
      tips <- tree$tip.label[phangorn::Descendants(tree, node, type = "tips")[[1]]]
      wrongTaxa <- e$retained %>%
        dplyr::filter(
          .data$label %in% tips,
          .data$rank == r,
          .data$taxon != !!taxon
        ) %>%
        dplyr::select("label", dplyr::one_of(names(method)))
      newAssign <- tibble::tibble(
        label = tips,
        rank = r,
        taxon
      )
      for (n in names(method)) {
        newAssign[[n]] <- unname(method[n])
      }
      # remove assignments which are not consistent with the one we just chose
      e$tip_taxa <- dplyr::bind_rows(e$tip_taxa, newAssign)
      e$rejected <- dplyr::bind_rows(
        e$rejected,
        dplyr::filter(e$retained, .data$rank >= r) %>%
          dplyr::semi_join(wrongTaxa, by = names(wrongTaxa))
      )
      e$retained <- dplyr::bind_rows(
        dplyr::filter(e$retained, .data$rank < r),
        dplyr::filter(e$retained, .data$rank >= r) %>%
          dplyr::anti_join(wrongTaxa, by = names(wrongTaxa))
      )
    }
  }
}

#' Assign taxon labels to nodes in a tree when there is a consensus of IDs on descendent tips.
#'
#' @param tree (`ape::phylo()` object) A tree including the taxa to be
#' classified.  The tip labels should match the "label" column in `taxa`. If no
#' tree is given, a star tree will be used, which results in each tip being
#' assigned to the strict consensus of its primary assignments, if any.
#' @param taxa (`data.frame`) Taxon assignments for the taxa on the tree,
#' as returned by `taxtable()`.  Should have columns "label",
#' "rank", and "taxon", but may have other columns as well; see the `method`
#' argument.  Multiple assignments for each label are allowed, and
#' can be generated by using `rbind()` or `dplyr::bind_rows()`
#' on results from multiple calls to `taxonomy()` and
#' `taxtable()`.
#' @param ranks (`character`) names of ranks
#' used in the taxon assignments, from most to least inclusive. This can be
#' omitted if the "`rank`" column of "`taxa`" is an `integer` or ordered `factor`,
#' in which case the *smallest* value is the *most inclusive* rank, or
#' alternatively if the "`rank`" column of "`taxa`" is a `character` or `factor`
#' with values/levels from the set `r paste(default_ranks, collapse = ", ")`.
#' @param method (a single `character` string, or a named `character` vector)
#' how to identify different methods. See details.
#'
#' @details # Distinguishing different primary methods
#' Primary methods can be distinguished in three ways:
#' 1. Not at all.  For this option, `taxa` should not have a column named
#'   "`method`", and the the `method` argument to `phylotax()` should be
#'   `NULL`.
#' 2. A single column named "`method`" in `taxa`. If `method=NULL` but `taxa`
#'   has a "`method`" column, then this column is assumed to uniquely
#'   identify the methods. Assignments made by `phylotax()` will have
#'   "`PHYLOTAX`" in the method column. This value can be changed by setting
#'   the `method` argument to an unnamed `character` string, e.g.,
#'   "`method = 'consensus'`".
#' 3. Custom columns. If the `method` argument is a named character vector,
#'   then the names are taken to be columns in `taxa` (which must exist) and
#'   the values are taken to be the values for each column which should be
#'   used for PHYLOTAX annotations, e.g.,
#'   "`method = c(algorithm = "PHYLOTAX", region = "ITS2")`". PHYLOTAX will
#'   treat each unique combination of values in these columns as a distinct
#'   method.
#'
#' @return an S3 object with class "`phylotax`", with five elements:
#' * "`tip_taxa` a `tibble::tibble()` with the same format as `taxa`, containing
#'   taxonomy assignments made by PHYLOTAX to tips.
#' * "`node_taxa`" a `tibble::tibble()` with columns "`node`", "`label`",
#'   "`rank`" and "`taxon`" giving taxonomy assignments made by PHYLOTAX to
#'   internal nodes.
#' * "`rejected`" a `tibble::tibble()` with the same format as `taxa` giving
#'   primary assignments which have been rejected by PHYLOTAX.
#' * "`retained`" a `tibble::tibble()` with the same format as `taxa` giving
#'   primary assignments which have not been rehected by PHYLOTAX. These may
#'   contain inconsistencies that PHYLOTAX was unable to resolve.
#' * "`missing`" a `tibble::tibble()` with the same format as `taxa`, giving the
#'   primary assignments which have not been assessed by PHULOTAX because they
#'   have labels which are not present on the tree.
#'
#' @export
phylotax <- function(
  tree = NULL, taxa,
  ranks = NULL,
  method = if (utils::hasName(taxa, "method")) "PHYLOTAX" else NULL
) {
  if (is.null(tree)) tree <- ape::read.tree(
    text = paste0("(", paste(unique(taxa$label), collapse = ","), ");")
  )
  method <- check_method(taxa, method)
  taxa <- check_ranks(taxa, ranks)
  e <- new_phylotax_env(tree, count_assignments(taxa), ranks)
  ranks <- sort(unique(taxa$rank))
  phylotax_(tree, taxa, phangorn::getRoot(tree), ranks, method, e)
  for (member in c("missing", "retained", "rejected", "tip_taxa"))
    for (n in c("n_tot", "n_diff"))
      e[[member]][[n]] <- NULL
  structure(
    as.list(e),
    class = "phylotax"
  )
}

#' Simple phylogenetic tree for use in examples
#'
#' @return a [`phylo`][ape::read.tree()] object giving a very simple tree;
#' tip labels are the same as those in [example_taxa()].
#' @export
#'
#' @examples example_tree()
example_tree <- function() {
  tree <- ape::read.tree(text = "(A:1,((B:1,C:1):1,((E:1,F:1):1,D:1):1):1);")
  tree$node.label <- 1:tree$Nnode
  tree
}

#' Example taxonomy assignments
#'
#' @return a [tibble::tibble()] in the format used by [taxtable()] giving a
#' simple set of taxonomy assignments suitable for use in [phylotax()].
#' Tip labels are the same as those used in [example_tree()].
#' @export
#'
#' @examples example_taxa()
example_taxa <- function() {
  tibble::tibble(
    label = rep(LETTERS[1:6], 2),
    method = rep(c("XTAX", "YTAX"), each = 6),
    rank = "genus",
    taxon = c(NA, "Tax1", "Tax2", "Tax2", NA, NA,
              NA, "Tax2", "Tax2", "Tax1", NA, "Tax1")
  ) %>%
    dplyr::filter(stats::complete.cases(.))
}