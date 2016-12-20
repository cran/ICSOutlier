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
dist.simu.test(object, index, m = 10000, level = 0.025, ...)
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
  \item{\dots}{further arguments passed on to the function \code{\link{quantile}}.}
}
\details{
The function extracts basically the dimension of the data from the \code{ics2} object and simulates \code{m} times, from a multivariate standard normal distribution,
the squared ICS distances with the components specified in \code{index}. The resulting value is then the mean of the \code{m} correponding quantiles of these distances 
at level 1-\code{level}.

Note that the function is seldomly called directly by the user but internally by \code{\link{ics.outlier}}.
}
\value{
A vector with the values of the (1-\code{level})th quantile.

}

\references{
 \cite{Archimbaud, A., Nordhausen, K. and Ruiz-Gazen, A. (2016), Multivariate Outlier Detection With ICS, <https://arxiv.org/abs/1612.06118>.}   
}
\author{
Aurore Archimbaud and Klaus Nordhausen
}

\seealso{
\code{\link[ICS]{ics2}}, \code{\link{ics.distances}}
}
\examples{
Z <- rmvnorm(1000, rep(0, 6))
Z[1:20, 1] <- Z[1:20, 1] + 10
A <- matrix(rnorm(36), ncol = 6)
X <- tcrossprod(Z, A)

pairs(X)
icsX <- ics2(X)

icsX.dist.1 <- ics.distances(icsX, index = 1)
CutOff <- dist.simu.test(icsX, 1, m = 500)

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
CutOff <- dist.simu.test(icsReliability, 1:2, m = 50, level = 0.02)
# Identification of the outliers based on the cut-off value
plot(ics.dist.scree)
abline(h = CutOff)
outliers <- which(ics.dist.scree >= CutOff)
text(outliers, ics.dist.scree[outliers], outliers, pos = 2, cex = 0.9)

}
\keyword{multivariate}
