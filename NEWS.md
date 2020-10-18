# phylotax (development version)

* Additional utility function `relabel_tree()` maps the tip labels of a tree
  to new values.
* Fixed a bug where PHYLOTAX would not make assignments inside an implicit
  *incertae sedis* taxon; e.g. assignments exist for a order and genus, but
  not for family, because there is genuine taxonomic uncertainty about which
  family the genus belongs in.
* Node labels are no longer used in trace output, this was very confusing when
  the node labels were actually bootstrap values (a common situation).
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
