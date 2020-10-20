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
  maxranks <- dplyr::group_by_at(taxa, c("label", names(method))) %>%
    dplyr::summarize(..max_rank.. = max(.data$rank)) %>%
    dplyr::ungroup()
  ranks <- sort(unique(taxa$rank))
  for (r in ranks) {
    # for converts factors to characters
    if (is.factor(ranks)) {
      r <- factor(r, levels = levels(ranks), ordered = is.ordered(ranks))
    }
    # it will be faster on each loop if we do this filter destructively now.
    maxranks <- dplyr::filter(
      maxranks,
      r < .data$..max_rank..
    )
    incertae_taxa <- maxranks
    incertae_taxa$rank <- r
    incertae_taxa <- dplyr::anti_join(
      incertae_taxa,
      taxa,
      by = c("label", "rank", names(method))
    )
    incertae_taxa$taxon <- "..phylotax_placeholder.."
    incertae_taxa$..max_rank.. <- NULL
    out <- dplyr::bind_rows(
      out,
      incertae_taxa
    )
  }
  out
}

# remove the dummy taxa added by interpolate_ranks
deinterpolate_ranks <- function(taxa) {
  dplyr::filter(taxa, .data$taxon != "..phylotax_placeholder..")
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
    node_assigned = tibble::tibble(
      node = integer(),
      label = NULL,
      rank = taxa$rank[FALSE],
      taxon = character()
    ),
    assigned = dplyr::filter(taxa, FALSE),
    retained = dplyr::filter(taxa, .data$label %in% tree$tip.label),
    rejected = dplyr::filter(taxa, FALSE),
    missing = dplyr::filter(taxa, !.data$label %in% tree$tip.label),
    tree = tree
  )
}


phylotax_ <- function(tree, taxa, node, ranks, method, e) {
  if (length(ranks) == 0) return()
  if (!is.null(tree$node.label)) {
    nodelabel <- tree$node.label[node - ape::Ntip(tree)]
    nodename <- sprintf("%d (label: %s)", node, nodelabel)
  } else {
    nodelabel <- as.character(node)
    nodename <- nodelabel
  }
  parents <- phangorn::Ancestors(tree, node, type = "all")
  for (r in ranks) {
    if (is.ordered(ranks)) r <- ordered(r, levels = levels(ranks))
    if (any(e$node_assigned$node %in% parents & e$node_assigned$rank == r)) next
    taxon <- clade_taxon(tree, e$retained, node, r)
    if (is.na(taxon)) {
      futile.logger::flog.debug("Could not assign a %s to node %s.", r, nodename)
      for (n in phangorn::Children(tree, node)) {
        phylotax_(tree, e$retained, n, ranks, method, e)
      }
      break
    } else {
      children <- phangorn::Descendants(tree, node, "tips")[[1]]
      if (taxon != "..phylotax_placeholder.."){
        if (length(children) > 0) {
          futile.logger::flog.info(
            "Assigned node %s and its %d descendant(s) to %s %s.",
            nodename, length(children), as.character(r), taxon)
        } else {
          futile.logger::flog.info("Assigned node %s to %s %s.", node,
                                   as.character(r), nodename)
        }
      }
      ranks <- ranks[-1]
      e$node_assigned <- dplyr::bind_rows(
        e$node_assigned,
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
      e$assigned <- dplyr::bind_rows(e$assigned, newAssign)
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
#' @param cons_method (a single `character` string, or a named `character` vector)
#' how to identify different methods for the fallback consensus. See details.
#' @param fallback (`logical`) If `TRUE`, use [lca_consensus()] for tips which
#' do not have a tree (which may be all of the tips, if no tree is given).
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
#' * "`assigned` a `tibble::tibble()` with the same format as `taxa`, containing
#'   taxonomy assignments made by PHYLOTAX to tips.
#' * "`node_assigned`" a `tibble::tibble()` with columns "`node`", "`label`",
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
  method = if (utils::hasName(taxa, "method")) "PHYLOTAX" else NULL,
  cons_method = method,
  fallback = TRUE
) {
  checkmate::check_flag(fallback)
  if (is.null(tree)) {
    stopifnot(isTRUE(fallback) || !is.null(tree))
    futile.logger::flog.info("'tree' is NULL, falling back on LCA consensus...")
    return(lca_consensus(taxa, ranks, cons_method))
  }
  checkmate::assert_data_frame(taxa)
  checkmate::assert_subset(c("label", "rank", "taxon"), colnames(taxa))
  method <- check_method(taxa, method)
  taxa <- check_ranks(taxa, ranks)
  taxa <- interpolate_ranks(taxa, method)
  e <- new_phylotax_env(tree, count_assignments(taxa), ranks)
  ranks <- sort(unique(taxa$rank))
  phylotax_(tree, taxa, phangorn::getRoot(tree), ranks, method, e)
  for (member in c("missing", "retained", "rejected", "assigned", "node_assigned")) {
    e[[member]] <- deinterpolate_ranks(e[[member]])
    for (n in c("..phylotax_n_tot", "..phylotax_n_diff")) {
      e[[member]][[n]] <- NULL
    }
  }
  out <- structure(
    as.list(e),
    class = "phylotax"
  )
  if (fallback && nrow(out$missing)) {
    lca <- lca_consensus(out$missing, ranks = ranks, method = cons_method)
    out <- combotax(out, lca, cons_method)
  }
  out
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


# ape::keep.tip, but also return a map from new tip numbers
# to old tip numbers
keep_tips_carefully <- function(tree, tips, keep_fun = ape::keep.tip) {
  old_node_labels <- tree$node.label
  old_ntip <- ape::Ntip(tree)
  tree$node.label <- seq_len(tree$Nnode) + old_ntip
  tree <- keep_fun(tree, tips)
  key <- tibble::tibble(
    new_node = seq_len(tree$Nnode) + ape::Ntip(tree),
    old_node = tree$node.label
  )
  if (is.null(old_node_labels)) {
    tree$node.label <- NULL
  } else {
    tree$node.label <- old_node_labels[key$old_node - old_ntip]
  }
  list(
    tree = tree,
    key = key
  )
}

#' Keep only certain tips from a phylotax object
#'
#' This is analogous to [ape::keep.tip], but also takes care of the taxonomic
#' annotations associated with nodes and tips. It also works for
#' [phylotax][phylotax()] objects without an associated tree (if `mrca=FALSE`)
#' and for "tip" labels which are not missing from the tree, but present in
#' the taxonomic annotations.
#'
#' @param phylotax ([phylotax][phylotax()] object)
#' @param tips (`character` vector) Tip labels to keep
#' @param mrca (`logical`) If `TRUE`, also keep all tips descended from the most
#' recent common ancestor (according to the tree) of `tips`. If `FALSE`,
#' just keep `tips`.
#' @param invert (`logical`) If `TRUE`, then the named `tips` (and all
#' descendants of their MRCA if `mrca=TRUE`) are removed, and all other tips are
#' kept.
#'
#' @return A [phylotax][phylotax()] object with the tree and taxonomic
#' assignments pruned to only include the specified tips.
#' @export
keep_tips <- function(phylotax, tips, mrca = (!is.null(phylotax$tree)),
                      invert = FALSE) {
  checkmate::assert_class(phylotax, "phylotax")
  has_tree <- !is.null(phylotax$tree)
  checkmate::assert_flag(mrca)
  checkmate::assert_flag(invert)
  tree_tips <-
    if (has_tree) intersect(tips, phylotax$tree$tip.label) else character()
  if (isTRUE(mrca)) {
    assertthat::assert_that(
      has_tree,
      msg = "For mrca=TRUE, the supplied phylotax object must include a tree."
    )
    mrca_node <- ape::getMRCA(phylotax$tree, tree_tips)
    mrca_tips <- phangorn::Descendants(phylotax$tree, mrca_node)[[1]]
    mrca_tips <- phylotax$tree$tip.label[mrca_tips]
    tips <- union(tips, mrca_tips)
    tree_tips <- union(tree_tips, mrca_tips)
  }
  if (has_tree) {
    phylofun <- if(invert) ape::drop.tip else ape::keep.tip
    newtree <- keep_tips_carefully(phylotax$tree, tree_tips, phylofun)
    phylotax$tree <- newtree$tree
    phylotax$node_assigned <- dplyr::filter(
      phylotax$node_assigned,
      .data$node %in% newtree$key$old_node
    )
    phylotax$node_assigned$node <- plyr::mapvalues(
      phylotax$node_assigned$node,
      newtree$key$old_node,
      newtree$key$new_node,
      warn_missing = FALSE
    )
  }
  filterfun <- if (invert) (function(x, y) ! x %in% y) else magrittr::is_in
  purrr::modify_if(
    phylotax,
    ~ all(utils::hasName(., c("rank", "taxon", "label"))) &
      !utils::hasName(., "node"),
    ~ dplyr::filter(., filterfun(.data$label, tips))
  )
}

#' Last common ancestor consensus
#'
#' This is an alternative to [phylotax()] which does not require a phylogenetic
#' tree. Instead, multiple assignments are resolved by strict consensus at each
#' rank in turn.  For eaxample, if for a particular sequence all available
#' assignments at the kingdom, phylum, and class agree, but there is
#' disagreement at the order level, then the consensus assignment takes the
#' kingdom, phylum, and class assignments to be correct, but does not include
#' any assignment at the order level or below.
#'
#'
#' @param taxa (`data.frame`) Taxon assignments for the taxa on the tree,
#' as returned by `taxtable()`.  Should have columns "`label`",
#' "`rank`", and "`taxon`", but may have other columns as well; see the `method`
#' argument.  Multiple assignments for each label are allowed, and
#' can be generated by using `rbind()` or `dplyr::bind_rows()`
#' on results from multiple calls to `taxonomy()` and
#' `taxtable()`.
#' @param ranks (`character`) names of ranks
#' used in the taxon assignments, from most to least inclusive. This can be
#' omitted if the "`rank`" column of "`taxa`" is an `integer` or ordered
#' `factor`, in which case the *smallest* value is the *most inclusive* rank, or
#' alternatively if the "`rank`" column of "`taxa`" is a `character` or `factor`
#' with values/levels from the set `r paste(default_ranks, collapse = ", ")`.
#' @param method (a single `character` string, or a named `character` vector)
#' how to identify different methods. See details.
#'
#' @details # Distinguishing different primary methods
#' Primary methods can be distinguished in three ways:
#' 1. Not at all.  For this option, `taxa` should not have a column named
#'   "`method`", and the the `method` argument to `lca_consensus()` should be
#'   `NULL`.
#' 2. A single column named "`method`" in `taxa`. If `method=NULL` but `taxa`
#'   has a "`method`" column, then this column is assumed to uniquely
#'   identify the methods. Assignments made by `lca_consensus()` will have
#'   "`consensus`" in the method column. This value can be changed by setting
#'   the `method` argument to an unnamed `character` string, e.g.,
#'   "`method = 'cons'`".
#' 3. Custom columns. If the `method` argument is a named character vector,
#'   then the names are taken to be columns in `taxa` (which must exist) and
#'   the values are taken to be the values for each column which should be
#'   used for consensus annotations, e.g.,
#'   "`method = c(algorithm = "consensus", region = "ITS2")`". Each unique
#'   combination of values in these columns is treated as a distinct method.
#'
#' @return A [phylotax][phylotax()] object.
#' @export
lca_consensus <- function(
  taxa, ranks = NULL,
  method = if (utils::hasName(taxa, "method")) "consensus" else NULL
) {
  checkmate::assert_data_frame(taxa)
  checkmate::assert_subset(c("label", "rank", "taxon"), colnames(taxa))
  method <- check_method(taxa, method)
  taxa <- dplyr::select(taxa, "label", dplyr::one_of(names(method)), "rank",
                        "taxon")
  taxa <- check_ranks(taxa, ranks)
  taxa <- interpolate_ranks(taxa, method)
  taxa <- count_assignments(taxa)
  assigned <- dplyr::group_by_at(taxa, c("label", names(method))) %>%
    dplyr::arrange(.data$rank) %>%
    dplyr::filter(dplyr::cumall(.data$..phylotax_n_diff == 1)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-dplyr::one_of(c(names(method), "..phylotax_n_diff", "..phylotax_n_tot"))) %>%
    unique()
  for (n in names(method)) {
    assigned[[n]] <- unname(method[n])
  }
  taxa <- dplyr::select(taxa, -"..phylotax_n_diff", -"..phylotax_n_tot")
  taxa <- deinterpolate_ranks(taxa)
  assigned <- deinterpolate_ranks(assigned)
  structure(
    list(
      node_assigned = NULL,
      tree = NULL,
      assigned = assigned,
      retained = taxa, # LCA retains everything
      rejected = dplyr::filter(taxa, FALSE), # LCA never rejects anything
      missing = dplyr::filter(taxa, FALSE) # nothing is missing from LCA
    ),
    class = "phylotax"
  )
}

combotax <- function(
  phylotax,
  lca = NULL,
  method = if (utils::hasName(phylotax$assigned, "method")) "PHYLOTAX" else NULL
) {
  method <- check_method(phylotax$assigned, method)
  assertthat::assert_that(methods::is(phylotax, "phylotax"))
  if (is.null(lca)) {
    lca <- lca_consensus(phylotax$missing, method = method)
  } else {
    for (n in names(method))
      lca$assigned[[n]] <- unname(method[n])
  }
  phylotax$assigned <- dplyr::bind_rows(
    phylotax$assigned,
    dplyr::anti_join(lca$assigned, phylotax$assigned, by = "label")
  )
  phylotax$rejected <- unique(dplyr::bind_rows(phylotax$rejected, lca$rejected))
  phylotax$retained <- unique(dplyr::bind_rows(phylotax$retained, lca$retained))
  phylotax$missing <- purrr::reduce(
    list(phylotax$missing, phylotax$rejected, phylotax$retained),
    dplyr::anti_join, by = names(method)
  )
  phylotax
}

#' Extract a taxon by name from a `phylotax` object
#'
#' For a [phylotax][phylotax()] object with a tree, the taxon is extracted based
#' on node annotations.  If more than one node is annotated as belonging to the
#' given taxon (i.e., it seems to be para/polyphyletic in the tree) then the
#' "correct" clade can be chosen by providing "true members"; the result will
#' only include the clade(s) which contain the "true members".
#'
#' If no tree is present (meaning the `phylotax` object was generated using
#' [lca_consensus()]), and also for sequences which are not present in the tree,
#' then the chosen taxon is simply extracted using the assigned taxonomy.
#'
#' @param phylotax ([phylotax][phylotax()] object) The object to extract from.
#' @param taxon (`character` string) The name of the taxon to extract.
#' @param true_members (`character`) (optional) Labels of sequences known to
#' belong in the taxon.
#'
#' @return A [phylotax][phylotax()] object containing only the specified taxon.
#' @export
extract_taxon <- function(phylotax, taxon, true_members = NULL) {
  checkmate::assert_class(phylotax, "phylotax")
  checkmate::assert_string(taxon)
  tips <-
    dplyr::group_by_at(phylotax$assigned, "label") %>%
    dplyr::filter(!!taxon %in% .data$taxon) %>%
    dplyr::pull("label") %>%
    unique()
  if (!is.null(phylotax$tree)) {
    tips <- setdiff(tips, phylotax$tree$tip.label)
    nodes <- dplyr::filter(phylotax$node_assigned, !!taxon == .data$taxon)$node
    nodes <- unique(nodes)
    nodetips <- phangorn::Descendants(phylotax$tree, nodes, "tips")
    nodetips <- purrr::map(nodetips, ~ phylotax$tree$tip.label[.])
    if (!is.null(true_members)) {
      checkmate::assert_character(true_members)
      nodetips <- purrr::keep(nodetips, ~any(. %in% true_members))
    }
    tips <- union(
      unlist(nodetips),
      tips
    )
  }
  keep_tips(phylotax, tips, mrca = FALSE)
}
