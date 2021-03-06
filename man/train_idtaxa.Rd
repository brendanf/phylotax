% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxonomy.R
\name{train_idtaxa}
\alias{train_idtaxa}
\title{Train an IDTAXA model to a fasta file with \[UV\]SEARCH/SINTAX-style taxonomy}
\usage{
train_idtaxa(fasta)
}
\arguments{
\item{fasta}{(\code{character} string) path to a fasta file containing the
reference sequences, with headers formatted as required for SINTAX.}
}
\value{
an object of classes \code{Taxa} and \code{Train}, as required for
\code{taxonomy_idtaxa()} or \code{DECIPHER::IdTaxa()}
}
\description{
IDTAXA needs to train a model on a taxonomic reference database before it can
be used to classify sequences.  This step can be time consuming, so it is
performed separately from the actual classification, and the result can be
saved for future analyses. This is a convenience function that makes it easy
to fit a model using \code{DECIPHER::LearnTaxa()} on the same reference
file that would be used for \code{taxonomy_sintax()}.
}
