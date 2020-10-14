# phylotax 0.0.2

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
