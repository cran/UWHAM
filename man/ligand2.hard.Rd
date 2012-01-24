\name{ligand2.hard}
\docType{data}
\alias{ligand2.hard}

\title{A dataset of binding energies obtained with a conventional interaction potential}
\description{A dataset for protein-ligand binding with an unmodified energy function.}

\usage{data(ligand2.hard)}

\format{A data frame containing 1 column and 18000 rows.}

\details{ 
The dataset contains binding energies (in kcal/mol) for the ligand
2/FKBP complex with an unmodified
energy function, simulated using replica exchange with 18 lambda values
(Tan et al. 2012).

Suppose that the dataset is converted to a \code{1000 x 18} matrix as
follows.
\preformatted{
  data(ligand2.hard)
  lig.data <- matrix(ligand2.hard$V1, nrow=1000)
}
The 18 columns of \code{lig.data}
correspond to the raw output of the 18 replicas. A replica is defined as a continuous
replica exchange molecular dynamics thread with variable lambda (exchanges are
performed by swapping lambda values rather than coordinates). The 1000
rows represent as many sequential time points during the replica
exchange trajectory. Each row contains binding energies from each
replica taken at the same simulation time (synchronous sampling).
}

\references{

Tan, Z., Gallicchio, E., Lapelosa, M., and Levy, R.M. (2012) "Theory of binless multi-state free energy
estimation with applications to protein-ligand binding," \emph{Journal of Chemical Physics}, 136, 144102.
}

\keyword{datasets}
