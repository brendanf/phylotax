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

#' Assign taxonomy to nucleotide sequences
#'
#' This method dispatches to other R packages or external programs.
#'
#' @param seq.table (\code{character} vector or something that can be coerced to
#'        one, or a matrix with sequences as the column names) Sequences to
#'        assign taxonomy
#' @param reference (\code{character} string giving a path to a file or the
#'        result from \code{\link[DECIPHER]{LearnTaxa}}) An appropriately
#'        formatted reference database
#' @param method (\code{character} string) taxonomy assignment method.
#'        Currently accepted values are "dada2", "sintax", and "idtaxa".
#' @param min_confidence (\code{integer} between 0 and 100) The minimum
#'        confidence to report results.
#' @param multithread (\code{integer} scalar) the number of processors to use
#'        for assignment.
#' @param ... additional arguments for methods
#'
#' @return raw results of the taxonomy assignment, of various types depending on
#'        \code{method}.  (\code{method = "dada2"}) returns a \code{list} with
#'        elements "tax" and "boot", as \code{\link[dada2]{assignTaxonomy}}.
#'        (\code{method = "sintax"}) returns a \code{data.frame} with columns
#'        "label", "hit", "strand", and, if \code{min_confidence} is not given,
#'        "c12n" (short for "classification"). (\code{method = "idtaxa"}) gives
#'        an S4 object of classes "Taxa" and "Train".
#'        Any of these can be passed to \code{\link{taxtable}} to get a
#'        uniform format.
#' @export
taxonomy <- function(seq.table, reference, method, min_confidence = 50, multithread = FALSE, ...) {
  # we can take a community matrix (in which case the sequences are the column
  # names) or the sequences
  if (is.matrix(seq.table)) {
    nreads <- Matrix::colSums(seq.table)
    seq <- colnames(seq.table)
  } else {
    nreads <- 1
    seq <- as.character(seq.table)
    names(seq) <- names(seq.table)
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


#' @param tryRC (\code{logical} scalar) passed on to \code{\link[dada2]{assignTaxonomy}}
#' @param outputBootstraps (\code{logical} scalar) passed on to \code{\link[dada2]{assignTaxonomy}}
#' @param verbose (\code{logical} scalar) passed on to \code{\link[dada2]{assignTaxonomy}}
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

#' @param exec (\code{character} string) name of the executable to use for
#'       SINTAX search.  The default is "vsearch", but "usearch" should also
#'       work.  In either case, the executable should be installed and on the
#'       system path.
#' @rdname taxonomy
#' @export
taxonomy_sintax <- function(seq, reference, min_confidence = NULL, multithread = FALSE, exec = "vsearch", ...) {
  args <- character()
  if (assertthat::is.string(seq) && file.exists(seq)) {
    seqfile <- seq
    seq <- seqinr::read.fasta(
      file = seqfile,
      as.string = TRUE) %>%
      {tibble(seq = unlist(seqinr::getSequence(., as.string = TRUE)),
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
      if (is.null(names(seq))) names(seq) <- tzara::seqhash(seq)
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

  if (!missing(min_confidence)) args <- c(args, "--sintax_cutoff", min_confidence)
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

#' @param strand (\code{character} string) passed on to \code{\link[DECIPHER]{IdTaxa}}
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
  tax <- dplyr::mutate_at(tax, "c12n", sub, pattern = "(,?[dkpcofgs]:unidentified)+$", replacement = "") %>%
    dplyr::mutate(name = sub(c12n, pattern = ".*,?[dkpcofgs]:", replacement = "")) %>%
    tidyr::separate_rows(c12n, sep = ",") %>%
    dplyr::mutate_at("c12n", dplyr::na_if, "") %>%
    tidyr::separate(c12n, into = c("rank", "taxon"), sep = ":") %>%
    dplyr::mutate_at("taxon", dplyr::na_if, "unidentified") %>%
    dplyr::mutate_at("rank", rank_factor, abbrev = TRUE) %>%
    tidyr::spread(key = "rank", value = "taxon") %>%
    dplyr::select(-`<NA>`) %>%
    dplyr::mutate_at("name", rlang::`%|%`, "unknown") %>%
    tidyr::unite(col = "Taxonomy",
                 dplyr::one_of("domain", "kingdom", "phylum", "class", "order",
                               "family", "genus", "species"),
                 sep = ";",
                 remove = FALSE) %>%
    dplyr::mutate_at("Taxonomy", sub, pattern = "(;?NA)+", replacement = "") %>%
    dplyr::mutate_at("Taxonomy", sub, pattern = "[|].*", replacement = "") %>%
    dplyr::mutate_at("Taxonomy", sub, pattern = "_sp^", replacement = "") %>%
    dplyr::mutate_at("Taxonomy", sub, pattern = "([^;]+);\\1", replacement = "\\1")

  if ("species" %in% names(tax)) {
    tax <- dplyr::mutate(tax, name = paste0(name, ifelse(is.na(species),
                                                         "_sp", "")))
  } else {
    tax <- dplyr::mutate(tax, name = paste0(name, "_sp"))
  }
  tax
}




#' Convert results from different taxonomic assignment algorithms to a uniform
#' format
#'
#' @param tax Results from \code{\link{taxonomy}}
#' @param min_confidence (\code{integer}) The minimum confidence to include in
#'   results. May be higher than the value given in \code{\link{taxonomy}},
#'   but will have no effect if it is lower.
#' @param ... passed to methods
#'
#' @return a \code{\link[tibble]{tibble}} with columns: \describe{
#'   \item{\code{label}}{sequence identifier}
#'   \item{\code{rank}}{the rank of the assignment}
#'   \item{\code{taxon}}{the taxon which was assigned}
#'   \item{\code{confidence}}{the confidence of the assignment}
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

taxtable_sintax <- function(tax, min_confidence = 0, ...) {
  tidyr::separate_rows(tax, hit, sep = ",") %>%
    dplyr::mutate_at("hit", dplyr::na_if, "") %>%
    tidyr::extract(hit, into = c("rank", "taxon", "confidence"),
                   regex = "([dkpcofgs]):([^(]+)\\(([01]\\.\\d+)\\)") %>%
    dplyr::mutate_at("taxon", dplyr::na_if, "unidentified") %>%
    dplyr::mutate_at("taxon", gsub, pattern = "Incertae",
                     replacement = "incertae") %>%
    dplyr::mutate_at("confidence", as.numeric) %>%
    dplyr::mutate_at("rank", rank_factor, abbrev = TRUE) %>%
    dplyr::select(label, rank, taxon, confidence) %>%
    dplyr::filter(confidence >= min_confidence, !is.na(taxon))
}

taxtable_idtaxa <- function(tax, min_confidence = 0, names = NULL, ...) {
  if (!missing(names) && !is.null(names)) names(tax) <- names
  purrr::imap_dfr(tax, ~tibble::tibble(label = .y,
                                       rank = rank_factor(.x$rank),
                                       taxon = gsub(" ", "_", .x$taxon),
                                       confidence = .x$confidence / 100)) %>%
    dplyr::mutate_at("taxon", gsub, pattern = "Incertae", replacement = "incertae") %>%
    dplyr::filter(rank != "rootrank", !startsWith(taxon, "unclassified_"),
                  confidence >= min_confidence)
}

taxtable_dada2 <- function(tax, names = rownames(tax$tax),
                           min_confidence = 0, ...) {
  taxa <- tax$tax %>% magrittr::set_rownames(names) %>%
    tibble::as_tibble(rownames = "label") %>%
    tidyr::gather(key = "rank", value = "taxon", -1)
  conf <- (tax$boot / 100) %>% magrittr::set_rownames(names) %>%
    tibble::as_tibble(rownames = "label") %>%
    tidyr::gather(key = "rank", value = "confidence", -1)
  dplyr::full_join(taxa, conf, by = c("label", "rank")) %>%
    dplyr::arrange(label) %>%
    dplyr::mutate_at("taxon", gsub, pattern = "Incertae", replacement = "incertae") %>%
    dplyr::filter(!is.na(taxon), confidence >= min_confidence) %>%
    dplyr::mutate_at("rank", tolower) %>%
    dplyr::mutate_at("rank", rank_factor)
}

#' Train an IDTAXA model to a fasta file with [UV]SEARCH/SINTAX-style taxonomy
#'
#' IDTAXA needs to train a model on a taxonomic reference database before it can
#' be used to classify sequences.  This step can be time consuming, so it is
#' performed separately from the actual classification, and the result can be
#' saved for future analyses. This is a convenience function that makes it easy
#' to fit a model using \code{\link[DECIPHER]{LearnTaxa}} on the same reference
#' file that would be used for \code{\link{taxonomy_sintax}}.
#'
#' @param fasta (\code{character} string) path to a fasta file containing the
#'   reference sequences, with headers formatted as required for SINTAX.
#'
#' @return an object of classes \code{Taxa} and \code{Train}, as required for
#'   \code{\link{taxonomy_idtaxa}} or \code{\link[DECIPHER]{IdTaxa}}
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
    dplyr::group_by(t, label, rank, n_reads) %>%
    dplyr::summarize(
      taxon =
        table(taxon) %>%
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
    dplyr::group_by(label, n_reads) %>%
    dplyr::arrange(rank) %>%
    dplyr::summarize(
      tip_label = paste(
        label[1],
        format(n_reads[1], width = 5),
        paste0(taxon, collapse = "-")
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
      taxa <- dplyr::filter(tax, label %in% tips, rank == !!rank)
      if (nrow(taxa) == 0) return(NA_character_)
    }
    tips <- tree$tip.label[phangorn::Descendants(tree, node, type = "tips")[[1]]]
  }
  taxa <- dplyr::filter(tax, label %in% tips, rank == !!rank) %>%
    dplyr::group_by(label)

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
      rank = taxa$rank[FALSE],
      taxon = character()
    ),
    tip_taxa = dplyr::filter(taxa, label %in% tree$tip.label),
    tree = tree
  )
}


phylotax_ <- function(tree, taxa, node, ranks, e) {
  if (length(ranks) == 0) return()

  parents <- phangorn::Ancestors(tree, node, type = "all")
  for (rank in ranks) {
    r <- rank_factor(rank)
    if (any(e$node_taxa$node %in% parents & e$node_taxa$rank == rank)) next
    taxon <- clade_taxon(tree, e$tip_taxa, node, r)
    if (is.na(taxon)) {
      futile.logger::flog.debug("Could not assign a %s to node %d.", rank, node)
      for (n in phangorn::Children(tree, node)) {
          phylotax_(tree, e$tip_taxa, n, ranks, e)
      }
      break
    } else {
      futile.logger::flog.info("Assigned node %d to %s %s.", node, rank, taxon)
      ranks <- ranks[-1]
      e$node_taxa <- dplyr::bind_rows(
        e$node_taxa,
        tibble::tibble(node = node, rank = r, taxon = taxon)
      )
      tips <- tree$tip.label[phangorn::Descendants(tree, node, type = "tips")[[1]]]
      wrongTaxa <- e$tip_taxa %>%
        dplyr::filter(
          label %in% tips,
          rank == r,
          taxon != !!taxon
        ) %>%
        dplyr::select(-rank, -taxon, -n_tot, -n_diff, -n_reads, -confidence, -n_method, -n_reference)
      newAssign <- tibble::tibble(
        label = tips,
        rank = r,
        taxon,
        method = "phylotax"
      )
      # remove assignments which are not consistent with the one we just chose
      e$tip_taxa <- dplyr::bind_rows(
        dplyr::filter(e$tip_taxa, rank < r),
        dplyr::filter(e$tip_taxa, rank >= r) %>%
          dplyr::anti_join(wrongTaxa, by = names(wrongTaxa)),
        newAssign
      )
    }
  }
}

#' Assign taxon labels to nodes in a tree when there is a consensus of IDs on descendent tips.
#'
#' @param tree (\code{\link[ape]{phylo}} object) A tree including the taxa to be
#'   classified.  The tip labels should match the "label" column in \code{taxa}.
#' @param taxa (\code{data.frame}) Taxon assignments for the taxa on the tree,
#'   as returned by \code{\link{taxtable}}.  Should have columns "label",
#'   "rank", and "taxon", but may have other columns as well.  Multiple
#'   assignments for each lable are allowed, and
#'   can be generated by using \code{\link{rbind}} or \code{\link[dplyr]{bind_rows}}
#'   on results from multiple calls to \code{\link{taxonomy}} and
#'   \code{\link{taxtable}}.
#'
#' @return a list with two elements, "tip_taxa" and "node_taxa".  "tip_taxa" is
#'     a \code{\link[tibble]{tibble}} with the same format as \code{taxa}, in
#'     which assignments which are inconsistent with the phylogeny have been
#'     removed, and new assignments deduced or confirmed from the phylogeny.
#'     These are identified by the value "phylotax" in the "method" column,
#'     which is created if it does not already exist.  "node_taxa" has columns
#'     "node", "rank" and "taxon", giving taxonomic assignments for the nodes of
#'     the tree.
#'
#' @export
phylotax <- function(tree = ape::read.tree(text = paste0("(", paste(unique(taxa$label), collapse = ","), ");")), taxa) {
  e <- new_phylotax_env(tree, taxa)
  ranks <- sort(unique(taxa$rank))
  phylotax_(tree, taxa, phangorn::getRoot(tree), ranks, e)
  as.list(e)
}