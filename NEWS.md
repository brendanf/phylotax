# phylotax 0.0.4

* Argument "`names`" is no longer required for `taxonomy(method = 'dada2')`,
  these are not taken directly from the dada2 taxonomy result (which will have
  names as long as the input sequences had names).
* New function `phylotax_to_phyloseq()` aids in followup analysis using the
  [phyloseq](https://joey711.github.io/phyloseq/index.html) package.
* New function `extract_taxon()` can be used to extract a particular taxon of
  interest from the results of `phylotax()` or `lca_consensus()`.
* `phylotax()` automatically calls `lca_consensus()` for any labels which are
  missing from the tree.  This can be disabled with `fallback=FALSE`.
* Node labels given in addition to node numbers, it was very confusing to use
  only the node labels when they were actually bootstrap values (a common
  situation).
* New function `lca_consensus()` calculates rankwise strict consensus.
  If no tree is give, `phylotax()` now dispatched to `lca_consensus()`, which is
  MUCH faster than calculating `phylotax()` on a star tree.
* New function `make_taxon_labels()` summarizes taxonomic assignments so that
  they can fit as the tip labels on a tree.
* **BREAKING CHANGE** Renamed the `phylotax` class elements "`tip_taxa`" ->
  "`assigned`" and "`node_taxa`" -> "`node_assigned`", to be more consistent
  with the other class elements (`accepted`, `rejected`, `missing`)
* Removed dependence on `stringr`.
* Added function `keep_tips()` to make subsets of the tips in a `phylotax`
  object, including the tree and the taxonomic annotations.
* Additional utility function `relabel_tree()` maps the tip labels of a tree
  to new values.
* Fixed a bug where PHYLOTAX would not make assignments inside an implicit
  *incertae sedis* taxon; e.g. assignments exist for a order and genus, but
  not for family, because there is genuine taxonomic uncertainty about which
  family the genus belongs in.
* `phylotax()` now counts all descendents, not just direct children, in its
  logged output. This does not change the actual results, just explains them
  correctly.

# phylotax 0.0.3

* **BREAKING CHANGE** `phylotax` returns taxonomic tables in four categories;
  * "`tip_taxa`" is the assignments which PHYLOTAX has made itself.
  * "`rejected`" are primary assignments which PHYLOTAX has rejected.
  * "`retained`" are primary assignments which PHYLOTAX has not rejected;
    however some of them may still be ambiguous.
  * "`missing`" are primary assignments whose labels are not present in the
    tree, so PHYLOTAX has not done anything with them. (But note that this will
    be empty if no tree was given).
* `phylotax()` now returns an S3 object of class "`phylotax`".  This should not
  break anything, and it allows the possibility of nice improvements in the
  future.
* `phylotax()$node_taxa` now includes a "`label`" column, and populates it with
  node labels if they exist, or just the numbers if they don't.  Node labels are
  also used in trace output.

# phylotax 0.0.2.1

* Two quick bugfix, applying to errors in `taxonomy_sintax()` and
  `taxtable_sintax()`.

# phylotax 0.0.2

* The undocumented requirement for the "`taxa`" argument to `phylotax()` to
  already include columns `n_tot` and `n_diff` is removed. `phylotax()` now
  generated these columns internally and deletes them when it is done, which
  will clobber these columns if they are present in the input.
* `phylotax()` gains a "`method`" argument, used to specify which of the
  columns in the input taxonomic assignment table are used to distinguish
  different primary methods, and what values should be used for assignments
  made by PHYLOTAX.
* `phylotax()` gains a "`ranks`" argument, in case the incoming data does not
  use the default rank names. This can be omitted if the `rank` column is
  already an ordered factor.
* Algorithm-specific implementations for `taxonomy()` are now exported.
* **BREAKING CHANGE** The first argument of `taxonomy()` is now called "seq"
  instead of "seq.table".
* Add `example_tree()` and `example_taxa()` for use in examples and tests.
* Fixed a bug due to a missing argument in `taxonomy_dada2`.
* **BREAKING CHANGE** Renamed `fit_idtaxa` to `train_idtaxa`.
* Added a `NEWS.md` file to track changes to the package.
