\name{ics.outlier}
\alias{ics.outlier}
\title{
Outlier Detection Using ICS
}
\description{
In a multivariate framework outlier(s) are detected using ICS. The function works on an object of class \code{\link[ICS]{ics2}}
 and decides automatically about the number of invariant components to use to search for the outliers and the number of outliers
detected on these components. Currently the function is restricted to the case of searching outliers only on the first components.

}
\usage{
ics.outlier(object, method = "norm.test", test = "agostino.test", mEig = 10000,
  level.test = 0.05, adjust = TRUE, level.dist = 0.025, mDist = 10000,
  type = "smallprop", ncores = NULL, iseed = NULL, pkg = "ICSOutlier", 
  qtype = 7, ...)
}

\arguments{
  \item{object}{object of class \code{ics2} where both \code{S1} and \code{S2} are specified as functions.}
  \item{method}{name of the method used to select the ICS components involved to compute ICS distances. Options are 
        \code{"norm.test"} and \code{"simulation"}. Depending on the method either \code{\link{comp.norm.test}} or \code{\link{comp.simu.test}} are used.}
  \item{test}{name of the marginal normality test to use if \code{method = "norm.test"}. Possibilities are 
              \code{"jarque.test"}, \code{"anscombe.test"}, \code{"bonett.test"}, \code{"agostino.test"}, \code{"shapiro.test"}.
              Default is \code{"agostino.test"}.}
  \item{mEig}{number of simulations performed to derive the cut-off values for selecting the ICS components. Only if \code{method = "simulation"}.
			See \code{\link{comp.simu.test}} for details.}
  \item{level.test}{\code{level} for the \code{\link{comp.norm.test}} or \code{\link{comp.simu.test}} functions. The inital level for selecting the invariant coordinates.}
  \item{adjust}{logical. For selecting the invariant coordinates, the level of the test can be adjusted for each component to deal with multiple testing. See \code{\link{comp.norm.test}} and \code{\link{comp.simu.test}}
                for details. Default is TRUE.}
  \item{level.dist}{\code{level} for the \code{\link{dist.simu.test}} function.  The (1-\code{level})th quantile used to determine the cut-off value for the ICS distances.}
  \item{mDist}{number of simulations performed to derive the cut-off value for the ICS distances.
               See \code{\link{dist.simu.test}} for details.}
  \item{type}{currently the only option is \code{"smallprop"} which means that only the first ICS components can be selected. See  \code{\link{comp.norm.test}}
				or \code{\link{comp.simu.test}} for details.}
  \item{ncores}{number of cores to be used in \code{\link{dist.simu.test}} and \code{\link{comp.simu.test}}. If \code{NULL} or 1, no parallel computing is used. Otherwise \link[parallel]{makeCluster} with \code{type = "PSOCK"} is used.}
  \item{iseed}{If parallel computation is used the seed passed on to \code{\link[parallel]{clusterSetRNGStream}}. Default is NULL which means no fixed seed is used.}
  \item{pkg}{When using parallel computing, a character vector listing all the packages which need to be loaded on the different cores via \code{\link{require}}. Must be at least "ICSOutlier"
  and must contain the packages needed to compute the scatter matrices.}
  \item{qtype}{specifies the quantile algorithm used in \code{\link{quantile}}.}
  \item{\dots}{passed on to other methods.}
}
\details{
The ICS method has attractive properties for outlier detection in the case of a small proportion of outliers. As for PCA three steps have to be performed:
 (i) select the components most useful for the detection, (ii) compute distances as outlierness measures for all observation and finally (iii) label outliers using some cut-off value.

This function performs these three steps automatically:
  \describe{
  \item{(i)}{For choosing the components of interest two methods are proposed: \code{"norm.test"} based on some marginal normality tests (see details in \code{\link{comp.norm.test}})
  or \code{"simulation"} based on a parallel analysis (see details in \code{\link{comp.simu.test}}). These two approaches lie on the intrinsic property of ICS in case of a small proportion 
  of outliers with the choice of S1 "more robust" than S2, which ensures to find outliers on the first components. Indeed  when using \code{S1 = MeanCov} and \code{S2 = Mean3Cov4}, 
  the Invariant Coordinates are ordered according to their classical Pearson kurtosis values in decreasing order. The information to find the outliers should be then contained in the first
  k nonnormal directions.}
  \item{(ii)}{Then the ICS distances are computed as the Euclidian distances on the selected k centered components \eqn{Z_k}. }
  \item{(iii)}{Finally the outliers are identified based on a cut-off derived from simulations. If the distance of an observation exceeds the expectation under the normal model, 
  this observation is labeled as outlier (see details in \code{\link{dist.simu.test}}).}
  }
  
As a rule of thumb, the percentage of contamination should be limited to 10\% in case of a mixture of gaussian distributions and using the default combination of locations and scatters for ICS.
}
\value{
an object of class \code{icsOut}
}

\note{
Function \code{\link{ics.outlier}} reached the end of its lifecycle, please use \code{\link{ICS_outlier}} instead. In future versions, \code{\link{ics.outlier}} will be deprecated and eventually removed.
}
\references{
  \cite{Archimbaud, A., Nordhausen, K. and Ruiz-Gazen, A. (2018), ICS for multivariate outlier detection with application to quality control. Computational Statistics & Data Analysis, 128:184-199. ISSN 0167-9473.  <https://doi.org/10.1016/j.csda.2018.06.011>.}   
  
\cite{Archimbaud, A., Nordhausen, K. and Ruiz-Gazen, A. (2018), ICSOutlier: Unsupervised Outlier Detection for Low-Dimensional Contamination Structure. The R Journal, 10:234-250.  <doi:10.32614/RJ-2018-034>.}   
}
\author{
Aurore Archimbaud and Klaus Nordhausen
}

\seealso{
\code{\link[ICS]{ics2}}, \code{\link{comp.norm.test}}, \code{\link{comp.simu.test}}, \code{\link{dist.simu.test}},
\code{\link{icsOut-class}}
}
\examples{
# ReliabilityData example: the observations 414 and 512 are suspected to be outliers  
library(REPPlab)
data(ReliabilityData)
icsReliabilityData <- ics2(ReliabilityData, S1 = tM, S2 = MeanCov)
# For demo purpose only small mDist value, but as extreme quantiles 
# are of interest mDist should be much larger. Also number of cores used
# should be larger if available
icsOutlierDA <- ics.outlier(icsReliabilityData, level.dist = 0.01, mDist = 50, ncores = 1)
icsOutlierDA
summary(icsOutlierDA)
plot(icsOutlierDA)

\dontrun{
# For using several cores and for using a scatter function from a different package
# Using the parallel package to detect automatically the number of cores
library(parallel)
# ICS with MCD estimates and the usual estimates
# Need to create a wrapper for the CovMcd function to return first the location estimate
# and the scatter estimate secondly.
data(HTP)
library(rrcov)
myMCD <- function(x,...){
  mcd <- CovMcd(x,...)
  return(list(location = mcd@center, scatter = mcd@cov))
}
icsHTP <- ics2(HTP, S1 = myMCD, S2 = MeanCov, S1args = list(alpha = 0.75))
# For demo purpose only small m value, should select the first seven components
icsOutlier <- ics.outlier(icsHTP, mEig = 50, level.test = 0.05, adjust = TRUE, 
                          level.dist = 0.025, mDist = 50,
                          ncores =  detectCores()-1, iseed = 123, 
                          pkg = c("ICSOutlier", "rrcov"))
icsOutlier
}
  
# Exemple of no direction and hence also no outlier
set.seed(123)
X = rmvnorm(500, rep(0, 2), diag(rep(0.1,2)))
icsX <- ics2(X)
icsOutlierJB <- ics.outlier(icsX, test = "jarque", level.dist = 0.01, 
				level.test = 0.01, mDist = 100, ncores = 1)
summary(icsOutlierJB)
plot(icsOutlierJB)
rm(.Random.seed)

     

# Example of no outlier
set.seed(123)
X = matrix(rweibull(1000, 4, 4), 500, 2)
X = apply(X,2, function(x){ifelse(x<5 & x>2, x, runif(sum(!(x<5 & x>2)), 5, 5.5))}) 
icsX <- ics2(X)
icsOutlierAG <- ics.outlier(icsX, test = "anscombe", level.dist = 0.01, 
				level.test = 0.05, mDist = 100, ncores = 1)
summary(icsOutlierAG)
plot(icsOutlierAG)
rm(.Random.seed)
  

}
\keyword{multivariate}
