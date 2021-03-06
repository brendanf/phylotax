% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxonomy.R
\name{taxtable}
\alias{taxtable}
\alias{taxtable_sintax}
\alias{taxtable_idtaxa}
\alias{taxtable_dada2}
\title{Convert results from different taxonomic assignment algorithms to a uniform
format}
\usage{
taxtable(tax, ...)

taxtable_sintax(tax, min_confidence = 0, ...)

taxtable_idtaxa(tax, min_confidence = 0, seq_id = NULL, ...)

taxtable_dada2(tax, min_confidence = 0, ...)
}
\arguments{
\item{tax}{Results from \code{taxonomy()}}

\item{...}{passed to methods}

\item{min_confidence}{(\code{integer}) The minimum confidence to include in
results. May be higher than the value given in \code{taxonomy()},
but will have no effect if it is lower.}

\item{seq_id}{(\code{character} vector) names for the sequences; these will be the
values that end up in the "\code{label}" column. If not given explicitly, they are
taken from the taxonomy results if this is possible.}
}
\value{
a \code{tibble::tibble()} with columns: \describe{
\item{\code{label}}{sequence identifier}
\item{\code{rank}}{the rank of the assignment}
\item{\code{taxon}}{the taxon which was assigned}
\item{\code{confidence}}{the confidence of the assignment}
} Each query sequence will typically occupy several rows of the output, one
for each rank which was assigned.
}
\description{
Convert results from different taxonomic assignment algorithms to a uniform
format
}
