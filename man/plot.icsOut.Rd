\name{plot.icsOut}
\docType{methods}
\alias{plot-icsOut}
\alias{plot.icsOut}
\alias{plot-method}
\alias{plot,icsOut,missing-method}
\title{Distances Plot for an icsOut Object}

\description{
  Distances plot for an icsOut object visualizing the separation of the outliers from the good data points.
}
\usage{
\S4method{plot}{icsOut,missing}(x, pch.out = 16, pch.good = 4, col.out = 1, col.good = grey(0.5), 
     col.cut = 1, lwd.cut = 1, lty.cut = 1, xlab = "Observation Number", 
     ylab = "ICS distances", ...)
}

\arguments{
  \item{x}{object of class \code{icsOut}.}
  \item{pch.out}{ploting symbol for the outliers.}
  \item{pch.good}{plotting symbol for the `good' data points.}
  \item{col.out}{color for the outliers.}
  \item{col.good}{color for the `good' data points.}
  \item{col.cut}{color for cut-off line.}
  \item{lwd.cut}{lwd value for cut-off line.}
  \item{lty.cut}{lty value for cut-off line.}
  \item{xlab}{default x-axis label.}
  \item{ylab}{default y-axis label.}
  \item{\dots}{other arguments for \code{plot}}
}
\details{
 For the figure the IC distances are plotted versus their index. The cut-off value for distances is given as a horizontal line
 and all observations above the line are considered as outliers.
}


\author{Aurore Archimbaud and Klaus Nordhausen}

\seealso{\code{\link{icsOut-class}} and \code{\link{ics.outlier}}}

\examples{
    Z <- rmvnorm(1000, rep(0, 6))
    Z[1:20, 1] <- Z[1:20, 1] + 10
    A <- matrix(rnorm(36), ncol = 6)
    X <- tcrossprod(Z, A)
    icsX <- ics2(X)
    # For demonstation purposes mDist is small, should be larger for real data analysis
    icsXoutliers <- ics.outlier(icsX, mDist = 500)
    plot(icsXoutliers, col.out = 2)
}


\keyword{hplot}
\keyword{methods}
