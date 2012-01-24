\name{ligand2.soft}
\docType{data}
\alias{ligand2.soft}

\title{A dataset of binding energies with a soft-core interaction potential}
\description{A dataset for protein-ligand binding with a soft energy function.}

\usage{data(ligand2.soft)}

\format{A data frame containing 1 column and 15000 rows.}

\details{ 
The dataset contains binding energies (in kcal/mol) for the ligand
2/FKBP complex with a soft-core interaction energy function, 
simulated using replica exchange with 15 lambda values
(Tan et al. 2012).
  
Suppose that the dataset is converted to a \code{1000 x 15} matrix as follows.
\preformatted{
data(ligand2.soft)
lig.data <- matrix(ligand2.soft$V1, nrow=1000)
}
The 15 columns of \code{lig.data} correspond to the 15 lambda states
from smallest (lambda=0) to largest (lambda=1).  The data were obtained
by reordering by lambda the raw output of the molecular
dynamics replicas at each observation time. A replica is defined as a continuous
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