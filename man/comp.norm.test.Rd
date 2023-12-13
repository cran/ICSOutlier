\name{comp.norm.test}
\alias{comp.norm.test}

\title{Selection of Nonnormal Invariant Components Using Marginal Normality Tests}

\description{
Identifies invariant coordinates that are non normal using univariate normality tests.
}
\usage{
comp.norm.test(object, test = "agostino.test", type = "smallprop", level = 0.05, 
  adjust = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{    
  \item{object}{object of class \code{ics2} where both \code{S1} and \code{S2} are specified as functions. The sample size and the dimension of interest
                are also obtained from the object.}
  \item{test}{name of the normality test to be used. Possibilites are \code{"jarque.test"}, 
      \code{"anscombe.test"}, \code{"bonett.test"}, \code{"agostino.test"}, \code{"shapiro.test"}.
  Default is \code{"agostino.test"}.}
  \item{type}{currently the only option is \code{"smallprop"}. See details.}
  \item{level}{the initial level used to make a decision based on the test p-values. See details.}
  \item{adjust}{logical. If \code{TRUE}, the quantiles levels are adjusted. Default is \code{TRUE}. See details.}
}
\details{
Currently the only available \code{type} is \code{"smallprop"} which detects which of the components follow a univariately normal distribution. It starts 
from the first component and stops when a component is detected as gaussian. Five tests for univariate normality are available.

If \code{adjust = FALSE} all tests are performed at the same \code{level}. This leads however often to too many components.
Therefore some multiple testing adjustments might be useful. The current default adjusts the level for the jth component as   
\code{level}/j.

Note that the function is seldomly called directly by the user but internally by \code{\link{ics.outlier}}.

}
\value{
 A list containing:
  \item{index }{integer vector indicating the indices of the selected components.}
  \item{test }{string with the name of the normality test used.}
  \item{criterion }{vector of the p-values from the marginal normality tests for each component.}
  \item{levels}{vector of the levels used for the decision for each component.}
  \item{adjust}{logical. \code{TRUE} if adjusted.}
  \item{type}{\code{type} used.}
  %\item{m}{number of iterations \code{m} used in the simulations.}
}
\note{
Function \code{\link{comp.norm.test}} reached the end of its lifecycle, please use \code{\link{comp_norm_test}} instead. In future versions, \code{\link{comp.norm.test}} will be deprecated and eventually removed.
}

\references{
 \cite{Archimbaud, A., Nordhausen, K. and Ruiz-Gazen, A. (2018), ICS for multivariate outlier detection with application to quality control. Computational Statistics & Data Analysis, 128:184-199. ISSN 0167-9473.  <https://doi.org/10.1016/j.csda.2018.06.011>.}   
}
\author{
Aurore Archimbaud and Klaus Nordhausen
}

\seealso{
\code{\link[ICS]{ics2}}, \code{\link{comp.simu.test}}, \code{\link[moments]{jarque.test}}, \code{\link[moments]{anscombe.test}},
 \code{\link[moments]{bonett.test}}, \code{\link[moments]{agostino.test}}, 
 \code{\link{shapiro.test}}
}
\examples{
Z <- rmvnorm(1000, rep(0, 6))
# Add 20 outliers on the first component
Z[1:20, 1] <- Z[1:20, 1] + 10
pairs(Z)
icsZ <- ics2(Z)
# The shift located outliers can be displayed in one dimension
comp.norm.test(icsZ)
# Only one invariant component is non normal and selected.
comp.norm.test(icsZ, test = "bo")

# Example with no outlier
Z0 <- rmvnorm(1000, rep(0, 6))
pairs(Z0)
icsZ0 <- ics2(Z0)
# Should select no component
comp.norm.test(icsZ0, level = 0.01)$index

}
\keyword{multivariate}
