\name{HTP}
\alias{HTP}
\docType{data}
\title{
Production Measurements of High-Tech Parts - Full Rank Case
}
\description{
The \code{HTP} data set contains 902 high-tech parts designed for consumer products characterized by 88 tests. These tests are performed to ensure a high quality 
of the production. All these 902 parts were considered functional and have been sold.  However the two parts 581 and 619 
showed defects in use and were returned to the manufacturer by the customer. Therefore these two can be considered as outliers.
}
\usage{data("HTP")}
\format{
  A data frame with 902 observations and 88 numeric variables V.1 - V.88.
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
Anonymized data from a nondisclosed manufacturer.
}

\examples{
# HTP data: the observations 581 and 619 are considered as outliers
data(HTP)
outliers <- c(581, 619)
boxplot(HTP)


# Outlier detection using ICS
icsHTP <- ics2(HTP)
# Selection of components based on a Normality Test, for demo purpose only small mDist value,
# but as extreme quantiles are of interest mDist should be much larger.
# Also more cores could be used if available.
icsOutlierDA <- ics.outlier(icsHTP, test = "agostino.test", level.test = 0.05,
                            level.dist = 0.02, mDist = 50, ncores = 1)
icsOutlierDA
summary(icsOutlierDA)
plot(icsOutlierDA)
text(outliers, icsOutlierDA@ics.distances[outliers], outliers, pos = 2, cex = 0.9, col = 2)

\dontrun{
# Selection of components based on simulations
# This might take a while to run (around 30 minutes)
icsOutlierPA <- ics.outlier(icsHTP, method = "simulation", level.dist = 0.02,
level.test = 0.05, mEig = 10000, mDist = 10000)
icsOutlierPA
summary(icsOutlierPA)
plot(icsOutlierPA)
text(outliers, icsOutlierPA@ics.distances[outliers], outliers, pos = 2, cex = 0.9, col = 2)
}

}
\keyword{datasets}
