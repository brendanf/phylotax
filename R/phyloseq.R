#' Convert assigned taxonomy and tree to a phyloseq object
#'
#' The output will have [`tax_table`][phyloseq::tax_table()] defined based on
#' the *assigned* taxa in `phylotax`, and [`phy_tree`][phyloseq::phy_tree()]
#' defined based on the tree in `phylotax`, if one is present.
#'
#' The resulting object will contain all taxon labels found in *both*
#' `otu_table` and the tree associated with `phylotax`, if any. Taxon labels
#' with no assignments in `phylotax` will have `NA` for all ranks in the
#' output's taxonomy table.
#'
#' @param phylotax ([`phylotax`][phylotax()] object)
#' @param otu_table ([`phyloseq::otu_table-class`] object)
#' @param ...  Additional arguments to pass to [phyloseq::phyloseq()]. These
#' should generally be the result of [phyloseq::sample_data()], but could also
#' include a [`phylo`][ape::read.tree()] object if `phylotax` does not already
#' have one or if `use_tree=FALSE`.
#' @param use_tree (`logical`) If `FALSE`, do not use the tree from `phylotax`
#' even if one exists.
#'
#' @return A [`phyloseq::phyloseq-class`]
#' object.
#'
#' @export
phylotax_to_phyloseq <- function(phylotax, otu_table, ..., use_tree = TRUE) {
  checkmate::assert_class(phylotax, "phylotax")
  checkmate::check_flag(use_tree)
  if (!requireNamespace("phyloseq"))
    stop("'phyloseq' package is required for as_phyloseq.")
  checkmate::assert_class(otu_table, "otu_table")
  tax <- phylotax$assigned %>%
    dplyr::select("label", "rank", "taxon") %>%
    dplyr::mutate_all(as.character) %>%
    tidyr::pivot_wider(names_from = "rank", values_from = "taxon") %>%
    tidyr::complete(label = phyloseq::taxa_names(otu_table)) %>%
    tibble::column_to_rownames("label") %>%
    as.matrix() %>%
    phyloseq::tax_table()
  if (methods::is(phylotax$tree, "phylo") && use_tree){
    phyloseq::phyloseq(tax, otu_table, phylotax$tree, ...)
  } else {
    phyloseq::phyloseq(tax, otu_table, ...)
  }
}
