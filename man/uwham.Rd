\name{uwham}
\alias{uwham}

\title{
Unbinned weighted histogram analysis method (UWHAM) for estimating free energies
}

\description{
This function implements UWHAM for estimating log-normalizing constants
(or free energies) and expectations from multiple distributions (such as
multiple generalized ensembles) as described in Tan et al. (2012).
}

\usage{
uwham(label=NULL, logQ, size=NULL, base = NULL, init = NULL, 
      fisher = TRUE)
}

\arguments{
  \item{label}{A vector of length \code{N} of labels between \code{1} to \code{M}
          such that \code{label[i]=j} when \code{i}th observation is obtained
          from \code{j}th thermodynamic state; 
          either \code{label} or \code{size} must be provided;
          if \code{fisher=FALSE} and if \code{label=NULL}, then \code{label} 
          is set such that the first \code{size[1]} observations are assumed to be from state \code{1}, the
          next \code{size[2]} observations from state \code{2}, etc.}
  \item{logQ}{\code{N x M} matrix of log unnormalized densities (such as \eqn{1/kT}
  times negative potential energies), where
         \code{N} is the total sample size, i.e., \code{sum(size)}, and
         \code{M} is the number of thermodynamic states for which free energies are to be computed;
         the \code{i}th row of \code{logQ} correspond to \code{i}th observations and 
         the \code{j}th column correspond to \code{j}th thermodynamic state.}
  \item{size}{A vector of length \code{M}, giving the individual sample
  sizes for the \code{M} thermodynamic states, ordered as the columns of \code{logQ};
  if \code{NULL}, then \code{label} is required and used to compute \code{size}.}
  \item{base}{The baseline index, between \code{1} to \code{M}, for the
  thermodynamic state (with sample size \code{>0}) whose free energy is set to 0; if \code{NULL}, 
  then \code{base} is set to the first index \code{j} such that \code{size[j]>0}.}
  \item{init}{A vector of length \code{M}, giving the initial values of
  the log-normalizing constants (or log of the partition
  functions); if \code{NULL}, then \code{init} is set to the zero vector.}
  \item{fisher}{Logical; if \code{NULL}, no variance estimation; if \code{TRUE}, variance estimation is based
on Fisher information; if \code{FALSE}, variance estimation is based on the Sandwich variance formula (see the details).}
}

\value{
  \item{ze}{The vector of estimated log-normalizing constants (or log of the partition
  functions).}
  \item{ve}{The vector of estimated variances for \code{ze}, if \code{fisher!=NULL}.}
  \item{Ve}{The estimated variance-covariance matrix for \code{ze}, if \code{fisher!=NULL}.} 
  \item{W}{The \code{N x M} matrix of UWHAM weights for each of the
           \code{N} observations at each of the \code{M} thermodynamic states.} 
  \item{check}{The column averages of \code{W[,sampled]}; 
               the elements of \code{check} should be equal to \code{1} to indicate a valid convergence of \emph{trust}.}
  \item{out}{The output of \emph{trust} used to minimize the objective function; see help(trust).}
  \item{label}{Same as argument \code{label}.}
  \item{size}{Same as argument \code{size}.}
  \item{base}{Same as argument \code{base}.}
}

\details{
The UWHAM method results from a number of interesting, sometimes independent, developments in physics and
statistics.
See Kong et al. (2003) for a formal statistical treatment along with earlier references, and Tan et al. (2012)
for a
more accessible account, presenting the method as a binless extension of the Weighted Histogram Analysis
Method (WHAM)
widely known in physics and chemistry (e.g., Ferrenberg and Swendsen
1989). The possibility of obtaining free energies from a binless extension of
WHAM was noticed by various authors (e.g., Kumar et al. 1992; Newman & Barkema 1999). 
The binless method was later reintroduced
by Shirts & Chodera (2008) to the physics literature and was called the Multi-state
Bennet Acceptance Ratio method (MBAR) because it can be interpreted as a multi-state extension of the Bennet
Acceptance Ratio (BAR) method. An implementation of MBAR in the Python language 
developed by Shirts & Chodera is freely available (see
https://simtk.org/home/pymbar).
The UWHAM package, while adopting an alternative numerical approach,
provides identical point estimates compared to
the MBAR Python package. In addition, the UWHAM package provides 
variance estimation based on variance formulas without using generalized inverses or, for
correlated data, by block bootstrap.

A typical application of UWHAM involves the computation of the relative
free energies of a series of thermodynamic states differing in
environmental conditions (temperature for example) and/or Hamiltonian
parameters (such as the strength of biasing potentials) from
data collected at these thermodynamic states. The
method takes as input the reduced energies (such as the inverse temperature
times the negative of the potential energy for canonical ensembles
differing in temperature) of the observations at all
thermodynamic states of interest.
In addition to the free energies, the
output includes estimates of the thermodynamic weights of the
observations for all states
to be used for thermodynamic reweighting
calculations.

The UWHAM method is statistically optimal in yielding the smallest asymptotic variances, provided that the
individual samples 
are independent and the observations in each sample are also
independent (Tan 2004).

To compute point estimates, the method is implemented here by minimizing
a convex objective function, as described in Tan et al. (2012).
This approach can be more effective than solving the nonlinear equations by the self-consistency or the
Newton-Raphson algorithm.
Currently, the optimization is done by using the R package \emph{trust}. 

Variance estimation provided here is based on the Fisher information or the Sandwich variance formula, 
as presented in Tan et al. (2012).  
In contrast with the Sandwich formula, the Fisher information based formula does not require labels indicating   
which thermodynamic state each observation is obtained from (see the dataset \code{\link{ligand2.hard}}). 
The analytical variance formulas are consistent when the observations are considered independent.
Alternatively, variance estimation can be done by block bootstrap, implemented in \code{\link{uwham.boot}}.
}

\references{
Ferrenberg, A.M. and Swendsen, R.H. (1989) "Optimized Monte Carlo data analysis,"
\emph{Physics Review Letters}, 63, 1195-1198.

Kong, A., McCullagh, P., Meng, X.-L., Nicolae, D., and Tan, Z. (2003) "A theory of statistical models for
Monte Carlo integration" (with discussion), \emph{Journal of the Royal Statistical Society}, Ser. B, 65, 
585-618. 

Kumar, S., Bouzida, D., Swendsen, R.H., Kollman, P.A. and Rosenberg, J.M. (1992) 
"The Weighted Histogram Analysis Method for free-energy
calculations on biomolecules. I. The method," \emph{Journal of Computational Chemistry}, 13, 1011-1021.

Newman, M.E.J. and Barkema, G.T. (1999) \emph{Monte Carlo Methods in Statistical
Physics}, Oxford University Press, New York.

Shirts, M.R. and J. D. Chodera, J.D. (2008) "Statistically optimal analysis of samples from multiple
equilibrium states,"
\emph{Journal of Chemical Physics}, 129, 124105.

Tan, Z. (2004) "On a likelihood approach for Monte Carlo integration," \emph{Journal of the American
Statistical Association}, 99, 1027-1036.

Tan, Z., Gallicchio, E., Lapelosa, M., and Levy, R.M. (2012) "Theory of binless multi-state free energy
estimation with applications to protein-ligand binding," \emph{Journal of Chemical Physics}, 136, 144102.
}


\examples{
#######################################
############## example 1 ##############  
#######################################

# This example illustrates the calculation of the standard free energy
# of binding of a ligand to a protein receptor by means of an alchemical
# perturbation potential of the form lambda*binding_energy(r), where
# lambda is a scaling parameter (lambda=0 corresponds to the
# protein-ligand uncoupled state, and lambda=1 corresponds to the
# coupled state) and binding_energy(r) is the (solvent averaged)
# potential energy of conformation r of the complex relative to one in
# which the receptor and the ligand are rigidly displaced at infinite
# separation. See Gallicchio et al., J. Chem. Theory Comput. 6,
# 2961-2977 (2010), and Tan et al. J. Chem. Phys., 136, 144102 (2012).

# Inverse temperature beta, in kcal/mol^-1
bet <- 1.0/(0.001986209*300.0)

# negative reduced potential function -beta*lambda*binding_energy.
#  "x" is the binding energy of a structure of the complex and "lam" the
#  value of lambda at which to compute the reduced energy
npot.fcn <- function(x, lam)
   -bet*lam*x

# values of lambda for the calculation with a "hard-core" potential 
lam <- c(0.0, 0.000000001, 0.00000001, 0.0000001, 0.000001, 0.00001, 0.0001,
         0.001, 0.01, 0.1, 0.15, 0.25, 0.35, 0.5, 0.6, 0.75, 0.9, 1.0)

#number of alchemical states
m <- length(lam)

# load binding energies
data(ligand2.hard)
lig.data <- ligand2.hard$V1

# 1000 observations at each lambda state
size <- rep(1000, m)

# total sample size
N <- sum(size)

# compute negative potential of each observation at each lambda state
neg.pot <- matrix(0, N,m)
for (j in 1:length(lam))
  neg.pot[,j] <- npot.fcn(x=lig.data, lam=lam[j])

# estimate free energies using UWHAM
out <- uwham(logQ=neg.pot, size=size, fisher=TRUE)

# convergence diagnosis: the elements should be equal to 1 
out$check

# the "ze" values are dimensionless free energies, that is the
# log of the partition functions (log Z). 
# To obtain thermodynamic free energies, multiply by -kT.
# 0.71 kcal/mol is the standard state correction; see
# Lapelosa et al. J Chem Theory Comput, 8, 44-60 (2012).
-out$ze/bet + 0.71

# variances of free energies from Fisher information
sqrt(out$ve)/bet

# perform block bootstrap for free energies to take into account
# time correlations in the binding energy data
# To save time for package checking, this is not run.
#out.boot <- uwham.boot(proc.type="parallel", block.size=50, boot.size=100, 
#                       logQ=neg.pot, size=size)
#
#-out.boot$ze/bet + 0.71
#sqrt(out.boot$ve)/bet

# estimation of average binding energies and variances 
# at lambda = 0.6, 0.75, 0.9, 1.0
state <- 15:18
out.phi <- uwham.phi(phi=lig.data, state=state, out.uwham=out, fisher=TRUE)

out.phi$phi
out.phi$phi.v 

# block bootstrap for both free energies and expectations
# To save time for package checking, this is not run.
#out.boot <- uwham.boot(proc.type="parallel", block.size=50, boot.size=100, 
#                       logQ=neg.pot, size=size, 
#                       phi=lig.data, state=state)
#
#out.boot$phi
#out.boot$phi.v

################################################
## example 2 (unequal and zero sample sizes) ###
################################################

# same calculation as above but with a "soft-core" potential and
# illustrating the ability to compute free energies and expectations
# for states with unequal sample sizes including those with
# zero sample sizes (or states that have not been sampled).
# See Tan et al. J. Chem. Phys., 136, 144102 (2012).

rm(list=ls())

# inverse temperature
bet <- 1.0/(0.001986209*300.0)

# negative potential function 
npot.fcn <- function(x, lam)
   -bet*lam*x

# read data (soft core)
lam <- c(0.0, 0.001, 0.002, 0.004, 0.006, 0.008, 0.01, 0.02, 0.06, 0.1,
         0.25, 0.5, 0.75, 0.9, 1.0)
m <- length(lam)

data(ligand2.soft)
lig.data <- ligand2.soft$V1

### unequal and zero sample sizes
size <- c( rep(1000, 5), rep(500, 3), rep(0, 2), rep(1000,5))
subs <- c(rep(TRUE, 5000), rep(c(rep(TRUE,500),rep(FALSE,500)), 3),
          rep(FALSE, 2000), rep(TRUE,5000))
lig.data <- lig.data[subs]

N <- sum(size)

# compute negative potential
neg.pot <- matrix(0, N,m)
for (j in 1:length(lam))
  neg.pot[,j] <- npot.fcn(x=lig.data, lam=lam[j])

# estimate free energies 
out <- uwham(logQ=neg.pot, size=size, fisher=TRUE)

-out$ze/bet + 0.71
sqrt(out$ve)/bet

# block bootstrap for free energies, 
# pretending that the data are generated from independent chains.
# To save time for package checking, this is not run.
#out.boot <- uwham.boot(proc.type="indep", 
#                       block.size=rep(50,m-2), boot.size=100, 
#                       logQ=neg.pot, size=size)
#
#-out.boot$ze/bet + 0.71
#sqrt(out.boot$ve)/bet

#######################################
## example 3 (serial tempering data) ##
#######################################

rm(list=ls())

# inverse temperature
bet <- 1.0/(0.001986209*300.0)

# negative potential function 
npot.fcn <- function(x, lam)
   -bet*lam*x

# lambda states
lam <- c(0.0,0.001,0.002,0.004,0.005,0.006,0.008,0.01,0.02,0.04,
         0.07,0.1,0.25,0.5,0.55,0.6,0.65,0.7,0.75,0.8,
         0.85,0.9,0.95,1.0)
m <- length(lam)

# loads cyclooctanol dataset
data(cyclooctanol)
lig.data <- cyclooctanol$V2

# sample size
N <- length(lig.data)

# state labels based on lambda values
# note that labels=1:m, not 0:(m-1)
state.labels <- factor(cyclooctanol$V1, labels=1:m)

# compute negative potential
neg.pot <- matrix(0, N,m)
for (j in 1:m)
  neg.pot[,j] <- npot.fcn(x=lig.data, lam=lam[j])

# estimate free energies, note that size=NULL because label is given
out <- uwham(label=state.labels, logQ=neg.pot, fisher=TRUE)

# free energies as a function of lambda, 0.36 kcal/mol is a standard
# state correction
-out$ze/bet + 0.36
sqrt(out$ve)/bet

# block bootstrap for free energies, note that proc.type="serial"
# for simulated tempering data.
# To save time for package checking, this is not run.
#out.boot <- uwham.boot(proc.type="serial", block.size=10, boot.size=100, 
#                       label=state.labels, logQ=neg.pot)
#
#-out.boot$ze/bet + 0.36
#sqrt(out.boot$ve)/bet
}


\keyword{free energy}