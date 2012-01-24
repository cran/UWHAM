\name{obj.fcn}
\alias{obj.fcn}

\title{
The objective function for UWHAM
}

\description{
This function computes the objective function, its gradient and its hessian matrix for UWHAM.
}

\usage{
obj.fcn(ze, logQ, size, base)
}

\arguments{
  \item{ze}{A vector of log-normalizing constants (or free energies) for the sampled thermodynamic states.}
  \item{logQ}{The matrix of log of unnormalized density ratios for the sampled thermodynamic states over the baseline.}
  \item{size}{A vector giving the individual sample sizes for the sampled thermodynamic states.}
  \item{base}{The baseline index.}
}

\value{
  \item{value}{The value of the objective function.}
  \item{gradient}{The gradient of the objective function.}
  \item{hessian}{The hessian matrix of objective function.}
}

\details{
The objective function is convex as discussed in Tan et al. (2012). See also Gill et al. (1998) for related results on biased sampling models. 
}

\references{
Gill, R., Vardi, Y., and Wellner, J. (1998) "Large sample theory of empirical distributions in biased sampling models," \emph{Annals of Statistics}, 16, 1069--1112.

Tan, Z., Gallicchio, E., Lapelosa, M., and Levy, R.M. (2012) "Theory of binless multi-state free energy estimation with applications to protein-ligand binding," \emph{Journal of Chemical Physics}, 136, 144102.
}


\keyword{free energy}