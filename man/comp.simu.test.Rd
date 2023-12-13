\name{comp.simu.test}
\alias{comp.simu.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Selection of Nonnormal Invariant Components Using Simulations
}
\description{
Identifies invariant coordinates that are nonnormal using simulations under a standard multivariate normal model for a specific data setup and scatter combination.
}
\usage{
comp.simu.test(object, m = 10000, type = "smallprop", level = 0.05, 
  adjust = TRUE, ncores = NULL, iseed = NULL, pkg = "ICSOutlier", 
  qtype = 7, ...)
}

\arguments{
  \item{object}{object of class \code{ics2} where both \code{S1} and \code{S2} are specified as functions. The sample size and the dimension of interest
                are also obtained from the object.}
  \item{m}{number of simulations. Note that since extreme quantiles are of interest \code{m} should be large. }
  \item{type}{currently the only type option is \code{"smallprop"}. See details.}
  \item{level}{the initial level used to make a decision. The cut-off values are the (1-\code{level})th quantile of the eigenvalues obtained from simulations. See details.}
  \item{adjust}{logical. If \code{TRUE}, the quantiles levels are adjusted. Default is \code{TRUE}. See details. }
  \item{ncores}{number of cores to be used. If \code{NULL} or 1, no parallel computing is used. Otherwise \link[parallel]{makeCluster} with \code{type = "PSOCK"} is used.}
  \item{iseed}{If parallel computation is used the seed passed on to \code{\link[parallel]{clusterSetRNGStream}}. Default is NULL which means no fixed seed is used.}
  \item{pkg}{When using parallel computing, a character vector listing all the packages which need to be loaded on the different cores via \code{\link{require}}. Must be at least "ICSOutlier"
  and must contain the packages needed to compute the scatter matrices.}
   \item{qtype}{specifies the quantile algorithm used in \code{\link{quantile}}.}
  \item{\dots}{further arguments passed on to the function \code{\link{quantile}}.}
}
\details{
Based on simulations it detects which of the components follow a univariately normal distribution. More precisely it identifies the observed eigenvalues larger than the ones coming
from normal distributed data. \code{m} standard normal data sets are simulated using the same data size and scatters as specified in the \code{ics2} object.
The cut-off values are determined based on a quantile of these simulated eigenvalues. 

As the eigenvalues, aka generalized kurtosis values, of ICS are ordered it is natural to perform the comparison in a specific order depending on the purpose. 
Currently the only available \code{type} is \code{"smallprop"} so starting with the first component, the observed eigenvalues are successively compared to
these cut-off values. The precedure stops when an eigenvalue is below the corresponding cut-off, so when a normal component is detected.
 
If \code{adjust = FALSE} all eigenvalues are compared to the same (\code{1-level})th level of the quantile. This leads however often to too many selected components. 
Therefore some multiple testing adjustment might be useful. The current default adjusts the quantile for the jth component as \code{1-level}/j. 

Note that depending on the data size and scatters used this can take a while and so it is more efficient to parallelize computations.
Note also that the function is seldomly called directly by the user but internally by \code{\link{ics.outlier}}.

}
\value{
 A list containing:
  \item{index }{integer vector indicating the indices of the selected components.}
  \item{test }{string \code{"simulation"}.}
  \item{criterion }{vector of the cut-off values for all the eigenvalues.}
  \item{levels}{vector of the levels used to derive the cut-offs for each component.}
  \item{adjust}{logical. \code{TRUE} if adjusted.}
  \item{type}{\code{type} used.}
  \item{m}{number of iterations \code{m} used in the simulations.}
}
\note{
Function \code{\link{comp.simu.test}} reached the end of its lifecycle, please use \code{comp_simu_test()} instead. In future versions, \code{\link{comp.simu.test}} will be deprecated and eventually removed.
}

\references{
 \cite{Archimbaud, A., Nordhausen, K. and Ruiz-Gazen, A. (2018), ICS for multivariate outlier detection with application to quality control. Computational Statistics & Data Analysis, 128:184-199. ISSN 0167-9473.  <https://doi.org/10.1016/j.csda.2018.06.011>.}   
}
\author{
Aurore Archimbaud and Klaus Nordhausen
}

\seealso{
\code{\link[ICS]{ics2}}, \code{\link{comp.norm.test}}
}
\examples{
# For a real analysis use larger values for m and more cores if available

set.seed(123)
Z <- rmvnorm(1000, rep(0, 6))
# Add 20 outliers on the first component
Z[1:20, 1] <- Z[1:20, 1] + 10
pairs(Z)
icsZ <- ics2(Z)
# For demo purpose only small m value, should select the first component
comp.simu.test(icsZ, m = 400, ncores = 1)

\dontrun{
# For using two cores
# For demo purpose only small m value, should select the first component
comp.simu.test(icsZ, m = 500, ncores = 2, iseed = 123)
  
# For using several cores and for using a scatter function from a different package
# Using the parallel package to detect automatically the number of cores
library(parallel)
# ICS with MCD estimates and the usual estimates
# Need to create a wrapper for the CovMcd function to return first the location estimate
# and the scatter estimate secondly.
library(rrcov)
myMCD <- function(x,...){
  mcd <- CovMcd(x,...)
  return(list(location = mcd@center, scatter = mcd@cov))
}
icsZmcd <- ics2(Z, S1 = myMCD, S2 = MeanCov, S1args = list(alpha = 0.75))
# For demo purpose only small m value, should select the first component
comp.simu.test(icsZmcd, m = 500, ncores = detectCores()-1, 
               pkg = c("ICSOutlier", "rrcov"), iseed = 123)
}

# Example with no outlier
Z0 <- rmvnorm(1000, rep(0, 6))
pairs(Z0)
icsZ0 <- ics2(Z0)
#Should select no component
comp.simu.test(icsZ0, m = 400, level = 0.01, ncores = 1)
}
\keyword{multivariate}
