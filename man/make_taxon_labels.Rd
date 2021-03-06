% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxonomy.R
\name{make_taxon_labels}
\alias{make_taxon_labels}
\title{Make labels summarizing the taxonomy of each sequence}
\usage{
make_taxon_labels(t, cols = character(), abbrev = FALSE)
}
\arguments{
\item{t}{(\code{data.frame}) Taxonomy assignments to summarize. At a minimum,
should include columns "\code{label}", "\code{rank}", and "\code{taxon}".}

\item{cols}{(\code{character}) (optional) Additional columns from \code{t} which will
be included in the output. They should be identical for each value of
"\code{label}".}

\item{abbrev}{(\code{logical} or \code{data.frame}) If \code{TRUE}, use a standard set of
abbreviations for mycological taxon names.  Alternatively, define your own
abbreviations using a \code{data.frame} with columns "\code{pattern}" and
"\code{replacement}", which are passed on to \code{\link[=gsub]{gsub()}}.}
}
\value{
A \code{data.frame} giving the old and new taxon taxon labels.
}
\description{
Make labels summarizing the taxonomy of each sequence
}
\details{
The standard abbreviations are:\tabular{ll}{
   pattern \tab replacement \cr
   (mycota|mycetes|ales|aceae) \tab  \cr
   incertae_sedis \tab i_s \cr
   Fungi\\b \tab F \cr
   Basidio\\b \tab B \cr
   Asco\\b \tab A \cr
   Chytridio\\b \tab Chy \cr
   Mucoro\\b \tab Muc \cr
   Glomero\\b \tab Glo \cr
}
}
