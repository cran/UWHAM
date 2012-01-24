\name{UWHAM-package}
\alias{UWHAM-package}
\alias{UWHAM}
\docType{package}

\title{
A R package for the unbinned weighted histogram analysis method (UWHAM)
}

\description{
UWHAM for estimating log-normalizing constants (or free energies) and
expectations from multiple distributions (such as
multiple generalized ensembles).
}

\details{

The R package \code{UWHAM} -- version 1.0 can be used for two main tasks:
\itemize{
\item to estimate log-normalizing constants (or free energies) from multiple distributions,

\item to estimate expectations from multiple distributions.
}
There are 3 main functions:
\itemize{
\item \code{uwham()}: to estimate free energies and associated variances,

\item \code{uwham.phi()}: to estimate expectations and associated variances,

\item \code{uwham.boot()}: to estimate variances based on block bootstrap.


}

The package also provides a function, \code{histw()}, for plotting weighted histograms. 
}

\keyword{UWHAM}