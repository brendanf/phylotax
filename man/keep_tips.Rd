% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phylotax.R
\name{keep_tips}
\alias{keep_tips}
\title{Keep only certain tips from a phylotax object}
\usage{
keep_tips(phylotax, tips, mrca = (!is.null(phylotax$tree)), invert = FALSE)
}
\arguments{
\item{phylotax}{(\link[=phylotax]{phylotax} object)}

\item{tips}{(\code{character} vector) Tip labels to keep}

\item{mrca}{(\code{logical}) If \code{TRUE}, also keep all tips descended from the most
recent common ancestor (according to the tree) of \code{tips}. If \code{FALSE},
just keep \code{tips}.}

\item{invert}{(\code{logical}) If \code{TRUE}, then the named \code{tips} (and all
descendants of their MRCA if \code{mrca=TRUE}) are removed, and all other tips are
kept.}
}
\value{
A \link[=phylotax]{phylotax} object with the tree and taxonomic
assignments pruned to only include the specified tips.
}
\description{
This is analogous to \link[ape:keep.tip]{ape::keep.tip}, but also takes care of the taxonomic
annotations associated with nodes and tips. It also works for
\link[=phylotax]{phylotax} objects without an associated tree (if \code{mrca=FALSE})
and for "tip" labels which are not missing from the tree, but present in
the taxonomic annotations.
}
