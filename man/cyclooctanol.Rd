\name{cyclooctanol}
\docType{data}
\alias{cyclooctanol}

\title{A simulated tempering dataset of interaction parameters and binding energies}
\description{Time series of interaction parameters and corresponding
  binding energies from a Hamiltonian hopping simulation of molecular binding.}

\usage{data(cyclooctanol)}

\format{A data frame containing 2 columns and 35000 rows.}

\details{ 

The dataset contains ligand-receptor interaction
parameters (lambda) and binding energies for the equivalent of an
Hamiltonian hopping simulation of cyclooctanol binding to
beta-cyclodextrin. (Hamiltonian hopping is the same as simulated
tempering except that potential parameters are varied rather than
temperature.) This illustrative dataset was constructed
by concatenating the output of some of the replicas of a 
Hamiltonian replica exchange simulation of the same complex (Wickstrom
et al. 2013). Harmless discontinuities may be observed at the joining points.

The dataset contains 2 columns and 3500 rows. Rows correspond to
consecutive time points. The first column \code{cyclooctanol$V1}
contains lambda values (24 discrete values between 0 and 1) and the second column
\code{cyclooctanol$V2} records binding energy values in kcal/mol.
}

\references{

  Wickstrom, L., Gallicchio, E., He, P., and Levy, R.M. (2013)
  "Large scale affinity calculations of cyclodextrin host-guest complexes:
  Understanding the role of reorganization in the molecular recognition
  process," \emph{Journal of Chemical Theory and Computation}, submitted for publication.
}

\keyword{datasets}
