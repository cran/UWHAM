\name{uwham.phi}
\alias{uwham.phi}

\title{
Unbinned weighted histogram analysis method (UWHAM) for estimating expectations
}

\description{
This function implements UWHAM for estimating expectations for multiple distributions (such as
multiple generalized ensembles) as described in Tan et al. (2012).
}

\usage{
uwham.phi(phi, state, out.uwham, fisher = TRUE)
}

\arguments{
  \item{phi}{A vector of function (or observable) values on the pooled sample.}
  \item{state}{A vector of indices between \code{1} to \code{M} for the thermodynamic states under which expectations are to be computed.} 
  \item{out.uwham}{The output of \code{\link{uwham}}.}
  \item{fisher}{Logical; if \code{NULL}, no variance estimation; if \code{TRUE}, variance estimation is based
on Fisher information; if \code{FALSE} variance estimation is based on the Sandwich variance formula (see the details).}
}

\value{
  \item{phi}{The vector of estimated expectations.}
  \item{phi.v}{The vector of estimated variances for \code{phi}, if \code{fisher!=NULL}.}
  \item{phi.V}{The variance-covariance matrix for \code{phi}, if \code{fisher!=NULL}.}
}

\details{
The implementation is directly based on Tan et al. (2012). See the details for \code{\link{uwham}}.
}

\references{
Tan, Z., Gallicchio, E., Lapelosa, M., and Levy, R.M. (2012) "Theory of binless multi-state free energy estimation with applications to protein-ligand binding," \emph{Journal of Chemical Physics}, 136, 144102.
}

\examples{

#See the examples for uwham().

}


\keyword{expectation}

