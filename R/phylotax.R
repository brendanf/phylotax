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

# if the incoming taxonomy has _incertae sedis_ taxa (i.e., the taxonomic
# classfication skips ranks) then add dummy taxa to fill in those ranks
interpolate_ranks <- function(taxa, method) {
  out <- taxa
  taxa <- dplyr::group_by_at(taxa, c("label", names(method)))
  ranks <- sort(unique(taxa$rank))
  for (r in ranks) {
    # for converts factors to characters
    if (is.factor(ranks)) {
      r <- factor(r, levels = levels(ranks), ordered = is.ordered(ranks))
    }
    incertae_taxa <- dplyr::filter(
      taxa,
      !r %in% .data$rank,
      r > min(.data$rank),
      r < max(.data$rank)
    )
    incertae_taxa <- dplyr::summarize(
      incertae_taxa,
      rank = r,
      taxon = "..phylotax_placeholder.."
    )
    out <- dplyr::bind_rows(
      out,
      incertae_taxa
    )
  }
  out
}

# remove the dummy taxa added by interpolate_ranks
deinterpolate_ranks <- function(taxa) {
  dplyr::filter(taxa, taxon != "..phylotax_placeholder..")
}



count_assignments <- function(taxa) {
  dplyr::group_by_at(taxa, c("label", "rank")) %>%
    dplyr::mutate(
      ..phylotax_n_diff = dplyr::n_distinct(.data$taxon),
      ..phylotax_n_tot = dplyr::n()
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
  best_taxon <- unique(taxa$taxon[taxa$..phylotax_n_diff == 1])
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
      futile.logger::flog.debug("Could not assign a %s to node %d.", r, node)
      for (n in phangorn::Children(tree, node)) {
        phylotax_(tree, e$retained, n, ranks, method, e)
      }
      break
    } else {
      children <- phangorn::Descendants(tree, node, "tips")
      if (taxon != "..phylotax_placeholder.."){
        if (length(children) > 0) {
          futile.logger::flog.info(
            "Assigned node %d and its %d descendant(s) to %s %s.",
            node, length(children), as.character(r), taxon)
        } else {
          futile.logger::flog.info("Assigned node %d to %s %s.", node,
                                   as.character(r), taxon)
        }
      }
      ranks <- ranks[-1]
      e$node_taxa <- dplyr::bind_rows(
        e$node_taxa,
        tibble::tibble(
          node = node,
          label = nodelabel,
          rank = r,
          taxon = taxon
        )
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

#' Assign taxon labels to nodes in a tree when there is a consensus of IDs on
#' descendant tips.
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
  taxa <- interpolate_ranks(taxa, method)
  e <- new_phylotax_env(tree, count_assignments(taxa), ranks)
  ranks <- sort(unique(taxa$rank))
  phylotax_(tree, taxa, phangorn::getRoot(tree), ranks, method, e)
  for (member in c("missing", "retained", "rejected", "tip_taxa", "node_taxa")) {
    e[[member]] <- deinterpolate_ranks(e[[member]])
    for (n in c("..phylotax_n_tot", "..phylotax_n_diff")) {
      e[[member]][[n]] <- NULL
    }
  }
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