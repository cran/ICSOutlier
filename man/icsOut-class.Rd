\name{icsOut-class}
\docType{class}
\alias{icsOut-class}

\title{Class icsOut}
\description{A S4 class to store results from performing outlier detection in an ICS context.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("icsOut", ...)}. But usually objects are created by the function \code{\link{ics.outlier}}.
     
}
\section{Slots}{
     \describe{
    \item{\code{outliers}:}{Object of class \code{"integer"}. A vector containing ones for outliers and zeros for non outliers.}
    \item{\code{ics.distances}:}{Object of class \code{"numeric"}. Vector giving the squared ICS distances of the observations from the invariant coordinates centered with the location estimate specified in \code{S1}.}
    \item{\code{ics.dist.cutoff}:}{Object of class \code{"numeric"}. The cut-off for the distances to decide if an observation is outlying or not.}
    \item{\code{level.dist}:}{Object of class \code{"numeric"}. The level for deciding upon the cut-off value for the ICS distances.}
    \item{\code{level.test}:}{Object of class \code{"numeric"}. The inital level for selecting the invariant coordinates.}
    \item{\code{method}:}{Object of class \code{"character"}. Name of the method used to decide upon the number of ICS components.}
    \item{\code{index}:}{Object of class \code{"numeric"}. Vector giving the indices of the ICS components selected.}
    \item{\code{test}:}{Object of class \code{"character"}. The name of the normality test as specified in the function call.}
    \item{\code{criterion}:}{Object of class \code{"numeric"}. Vector giving the marginal levels for the components selection.}
    \item{\code{adjust}:}{Object of class \code{"logical"}. Wether the initial level used to decide upon the number of components has been adjusted for multiple testing or not.}
    \item{\code{type}:}{Object of class \code{"character"}. Currently always the string \code{"smallprop"}.}
    \item{\code{mDist}:}{Object of class \code{"integer"}. Number of simulations performed to decide upon the cut-off for the ICS distances.}
    \item{\code{mEig}:}{Object of class \code{"integer"}. Number of simulations performed for selecting the ICS components based on simulations.}
    \item{\code{S1name}:}{Object of class \code{"character"}. Name of S1 in the original ics2 object.}
    \item{\code{S2name}:}{Object of class \code{"character"}. Name of S2 in the original ics2 object.}
  }
                                  
}
\section{Methods}{
For this class the following generic functions are available: \code{\link{print.icsOut}}, \code{\link{summary.icsOut}} and \code{\link{plot.ics}}
}

\section{Note}{
In case no extractor function for the slots exists, the component can be extracted the usual way using '@'. This S4 class is created by  \code{\link{ics.outlier}} that reached the end of its lifecycle, please use  \code{\link{ICS_outlier}} instead for which an object of class S3 is returned. In future versions,  \code{\link{ics.outlier}} will be deprecated and eventually removed.
}

\author{Aurore Archimbaud and Klaus Nordhausen}



\seealso{\code{\link{ics.outlier}}}

\keyword{classes}
