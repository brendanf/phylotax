\docType{package}
\name{phylotax-package}
\alias{phylotax-package}
\title{phylotax: Refine taxonomic assignment of environmental sequences using a taxonomic tree}
\description{
Combines taxonomic assignments made with different methods or
    databases using a phylogenetic tree.  Clades within the tree are assigned
    to a taxon if the assignment is consistent with at least one result for all
    of the members of the clade, except the members that do not have any
    assignments at all.  This resolves conflict between different methods, while
    also adding assignments to unknown sequences which are nested in an
    identified clade.
    The package also includes helper functions to assign taxonomy using several
    methods, and compile the results into a uniform format.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://www.github.com/brendanf/phylotax}
  \item Report bugs at \url{https://www.github.com/brendanf/phylotax/issues}
}

}
\author{
\strong{Maintainer}: Brendan Furneaux \email{brendan.furneaux@gmail.com}

}
\keyword{internal}
