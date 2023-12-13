\name{ics.distances}
\alias{ics.distances}

\title{
Squared ICS Distances for Invariant Coordinates
}
\description{
Computes the squared ICS distances, defined as the Euclidian distances of the selected centered components.

}
\usage{
ics.distances(object, index = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{object of class \code{ics2} where both \code{S1} and \code{S2} are specified as functions.}
  \item{index}{vector of integers indicating the indices of the components to select.}
}
\details{
	For outlier detection, the squared ICS distances can be used as a measure of outlierness. Denote as \eqn{Z}
    the invariant coordinates centered with the location estimate specified in \code{S1} (for details see \code{\link{ics2}}).
    Let \eqn{Z_k} be the \eqn{k} components of \eqn{Z} selected by \code{index}, then the ICS distance of the observation \eqn{i} 
    is defined as: 
	\deqn{ICSD^2(x_i,k) = || Z_k||^2.}{ICSD^2(x_i,k) = || Z_k||^2.}
	
	Note that if all components are selected, the ICS distances are equivalent to the Mahlanobis distances computed with
	respect of the first scatter and associated location specified in \code{S1}. 
}
\value{
A numeric vector containing the squared ICS distances.
}

\note{
Function \code{ics.distances()} reached the end of its lifecycle, please use \code{\link{ics_distances}} instead. In future versions, \code{ics_distances()} will be deprecated and eventually removed.
}

\references{
  \cite{Archimbaud, A., Nordhausen, K. and Ruiz-Gazen, A. (2018), ICS for multivariate outlier detection with application to quality control. Computational Statistics & Data Analysis, 128:184-199. ISSN 0167-9473.  <https://doi.org/10.1016/j.csda.2018.06.011>.}   
}
\author{
Aurore Archimbaud and Klaus Nordhausen
}

\seealso{
\code{\link[ICS]{ics2}}, \code{\link{mahalanobis}}
}
\examples{
Z <- rmvnorm(1000, rep(0, 6))
Z[1:20, 1] <- Z[1:20, 1] + 5
A <- matrix(rnorm(36), ncol = 6)
X <- tcrossprod(Z, A)

pairs(X)
icsX <- ics2(X)

icsX.dist.all <- ics.distances(icsX, index = 1:6)
maha <- mahalanobis(X, center = colMeans(X), cov = cov(X))
# in this case the distances should be the same
plot(icsX.dist.all, maha)
all.equal(icsX.dist.all, maha)

icsX.dist.first <- ics.distances(icsX, index = 1)
plot(icsX.dist.first)



}
\keyword{multivariate}
