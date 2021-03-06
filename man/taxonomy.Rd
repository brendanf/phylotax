% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxonomy.R
\name{taxonomy}
\alias{taxonomy}
\alias{taxonomy_dada2}
\alias{taxonomy_sintax}
\alias{taxonomy_idtaxa}
\title{Assign taxonomy to nucleotide sequences}
\usage{
taxonomy(seq, reference, method, min_confidence = 50, multithread = FALSE, ...)

taxonomy_dada2(
  seq,
  reference,
  multithread = FALSE,
  min_confidence,
  tryRC = FALSE,
  outputBootstraps = TRUE,
  verbose = TRUE,
  ...
)

taxonomy_sintax(
  seq,
  reference,
  min_confidence = NULL,
  multithread = FALSE,
  exec = NULL,
  ...
)

taxonomy_idtaxa(
  seq,
  reference,
  multithread = FALSE,
  strand = "top",
  min_confidence = 40,
  ...
)
}
\arguments{
\item{seq}{(`character`` vector or something that can be coerced to
one, or a matrix with sequences as the column names ) Sequences to
assign taxonomy}

\item{reference}{(`character`` string giving a path to a file or the
result from \code{\link[DECIPHER:LearnTaxa]{DECIPHER::LearnTaxa()}}/\code{\link[=train_idtaxa]{train_idtaxa()}}) An appropriately
formatted reference database (see Details).}

\item{method}{(`character`` string) taxonomy assignment method.
Currently accepted values are "dada2", "sintax", and "idtaxa".}

\item{min_confidence}{(`integer`` between 0 and 100) The minimum
confidence to report results.}

\item{multithread}{(\code{integer} scalar) the number of processors to use
for assignment.}

\item{...}{additional arguments for methods}

\item{tryRC}{(\code{logical} scalar) passed on to \code{dada2::assignTaxonomy()}}

\item{outputBootstraps}{(\code{logical} scalar) passed on to \code{dada2::assignTaxonomy()}}

\item{verbose}{(\code{logical} scalar) passed on to \code{dada2::assignTaxonomy()}}

\item{exec}{(\code{character} string) name of the executable to use for
SINTAX search.  The default is "vsearch", but "usearch" should also
work.  In either case, the executable should be installed and on the
system path.}

\item{strand}{(\code{character} string) passed on to \code{DECIPHER::IdTaxa()}}
}
\value{
raw results of the taxonomy assignment, of various types depending on
\code{method}.
}
\description{
This method uses a common interface to call primary taxonomic assignment
algorithms (i.e., those which assign taxonomy based on a taxonomically
classified reference sequence database, but not based on the results of other
algorithms) from other R packages or external programs.
}
\section{Return types}{
\itemize{
\item \code{taxonomy_dada2} and \code{taxonomy(..., method = "dada2")} return a \code{list}
with elements "tax" and "boot", as \link[dada2:assignTaxonomy]{dada2::assignTaxonomy}.
\item \code{taxonomy_sintax} and \code{taxonomy(..., method = "sintax")} return a
\code{data.frame} with columns "label", "hit", "strand", and, if
\code{min_confidence} is not given, "c12n" (short for "classification").
\item \code{taxonomy_idtaxa} and \code{taxonomy(..., method = "idtaxa")} gives
an S4 object of classes "Taxa" and "Train".
Any of these can be passed to \link{taxtable} to get a uniform format suitable
for use in \link{phylotax}.
}
}

