\name{dist.simu.test}
\alias{dist.simu.test}
\title{
Cut-Off Values Using Simulations for the Detection of Extreme ICS Distances
}
\description{
Computes the cut-off values for the identification of the outliers based on the squared ICS distances.
It uses simulations under a multivariate standard normal model for a specific data setup and scatters combination.
}
\usage{
dist.simu.test(object, index, m = 10000, level = 0.025, ncores = NULL, 
               iseed = NULL, pkg = "ICSOutlier", qtype = 7, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object}{object of class \code{ics2} where both \code{S1} and \code{S2} are specified as functions. The sample size and the dimension of interest
                are also obtained from the object.}
  \item{index}{integer vector specifiying which components are used to compute the 
  \code{\link{ics.distances}}.}
  \item{m}{number of simulations. Note that extreme quantiles are of interest and hence \code{m} should be large. }
  \item{level}{the (1-\code{level}(s))th quantile(s) used to choose the cut-off value(s). Usually just one number between 0 and 1. However a vector is
  also possible.}
  \item{ncores}{number of cores to be used. If \code{NULL} or 1, no parallel computing is used. Otherwise \link[parallel]{makeCluster} with \code{type = "PSOCK"} is used.}
  \item{iseed}{If parallel computation is used the seed passed on to \code{\link[parallel]{clusterSetRNGStream}}. Default is NULL which means no fixed seed is used.}
  \item{pkg}{When using parallel computing, a character vector listing all the packages which need to be loaded on the different cores via \code{\link{require}}. Must be at least "ICSOutlier"
  and must contain the packages needed to compute the scatter matrices.}
  \item{qtype}{specifies the quantile algorithm used in \code{\link{quantile}}.}
  \item{\dots}{further arguments passed on to the function \code{\link{quantile}}.}
}
\details{
The function extracts basically the dimension of the data from the \code{ics2} object and simulates \code{m} times, from a multivariate standard normal distribution,
the squared ICS distances with the components specified in \code{index}. The resulting value is then the mean of the \code{m} correponding quantiles of these distances 
at level 1-\code{level}.

Note that depending on the data size and scatters used this can take a while and so it is more efficient to parallelize computations.

Note that the function is seldomly called directly by the user but internally by \code{\link{ics.outlier}}.
}
\value{
A vector with the values of the (1-\code{level})th quantile.

}

\note{
Function \code{\link{dist.simu.test}} reached the end of its lifecycle, please use \code{\link{dist_simu_test}} instead. In future versions, \code{\link{dist.simu.test}} will be deprecated and eventually removed.
}

\references{
 \cite{Archimbaud, A., Nordhausen, K. and Ruiz-Gazen, A. (2018), ICS for multivariate outlier detection with application to quality control. Computational Statistics & Data Analysis, 128:184-199. ISSN 0167-9473.  <https://doi.org/10.1016/j.csda.2018.06.011>.}   
}
\author{
Aurore Archimbaud and Klaus Nordhausen
}

\seealso{
\code{\link[ICS]{ics2}}, \code{\link{ics.distances}}
}
\examples{
# For a real analysis use larger values for m and more cores if available

Z <- rmvnorm(1000, rep(0, 6))
Z[1:20, 1] <- Z[1:20, 1] + 10
A <- matrix(rnorm(36), ncol = 6)
X <- tcrossprod(Z, A)

pairs(X)
icsX <- ics2(X)

icsX.dist.1 <- ics.distances(icsX, index = 1)
CutOff <- dist.simu.test(icsX, 1, m = 500, ncores = 1)

# check if outliers are above the cut-off value
plot(icsX.dist.1, col = rep(2:1, c(20, 980)))
abline(h = CutOff)


library(REPPlab)
data(ReliabilityData)
# The observations 414 and 512 are suspected to be outliers
icsReliability <- ics2(ReliabilityData, S1 = MeanCov, S2 = Mean3Cov4)
# Choice of the number of components with the screeplot: 2
screeplot(icsReliability)
# Computation of the distances with the first 2 components
ics.dist.scree <- ics.distances(icsReliability, index = 1:2)
# Computation of the cut-off of the distances 
CutOff <- dist.simu.test(icsReliability, 1:2, m = 50, level = 0.02, ncores = 1)
# Identification of the outliers based on the cut-off value
plot(ics.dist.scree)
abline(h = CutOff)
outliers <- which(ics.dist.scree >= CutOff)
text(outliers, ics.dist.scree[outliers], outliers, pos = 2, cex = 0.9)

\dontrun{
# For using three cores
# For demo purpose only small m value, should select the first component
dist.simu.test(icsReliability, 1:2, m = 500, level = 0.02, ncores = 3, iseed = 123)
  
# For using several cores and for using a scatter function from a different package
# Using the parallel package to detect automatically the number of cores
library(parallel)
# ICS with Multivariate Median and Tyler's Shape Matrix and the usual estimates
library(ICSNP)
icsReliabilityHRMest <- ics2(ReliabilityData, S1 = HR.Mest, S2 = MeanCov, 
                             S1args = list(maxiter = 1000))
# Computation of the cut-off of the distances. For demo purpose only small m value.
dist.simu.test(icsReliabilityHRMest,  1:2, m = 500, level = 0.02, ncores = detectCores()-1, 
               pkg = c("ICSOutlier","ICSNP"), iseed = 123)
}
}
\keyword{multivariate}
