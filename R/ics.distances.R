#' @export
ics.distances <-
function(object, index = NULL)
     {
    if (!inherits(object, "ics2")) stop("'object' must be of class 'ics2'" )
    if(is.null(index)) index <- 1:ncol(object@Scores)
    DIST <- rowSums(object@Scores[, index, drop = FALSE]^2)
    return(DIST)
}

#' Squared ICS Distances for Invariant Coordinates
#'
#' @param object object of class \code{"ICS"} where both \code{S1} and \code{S2} are specified as functions.
#' @param index vector of integers indicating the indices of the components to select.
#'
#'
#' @details
#' For outlier detection, the squared ICS distances can be used as a measure of outlierness. Denote as \eqn{Z} the invariant coordinates centered with the location estimate specified in \code{S1} (for details see [ICS()][ICS::ICS()]). Let \eqn{Z_k} be the \eqn{k} components of \eqn{Z} selected by \code{index}, then the ICS distance of the observation \eqn{i}  is defined as:     \deqn{ICSD^2(x_i,k) = || Z_k||^2.}{ICSD^2(x_i,k) = || Z_k||^2.}
#' 
#' Note that if all components are selected, the ICS distances are equivalent to the Mahalanobis distances computed with respect of the first scatter and associated location specified in \code{S1}. 
#' 
#' @return A numeric vector containing the squared ICS distances.
#' 
#' @references
#' Archimbaud, A., Nordhausen, K. and Ruiz-Gazen, A. (2018), ICS for multivariate outlier detection with application to quality control. Computational Statistics & Data Analysis, 128:184-199. ISSN 0167-9473.  \doi{10.1016/j.csda.2018.06.011}. 

#' @author Aurore Archimbaud and Klaus Nordhausen
#' 
#' @seealso [ICS()][ICS::ICS()], [mahalanobis()]
#' @export
#'
#' @examples
#' Z <- rmvnorm(1000, rep(0, 6))
#' Z[1:20, 1] <- Z[1:20, 1] + 5
#' A <- matrix(rnorm(36), ncol = 6)
#' X <- tcrossprod(Z, A)
#' 
#' pairs(X)
#' icsX <- ICS(X, center = TRUE)
#' 
#' icsX.dist.all <- ics_distances(icsX, index = 1:6)
#' maha <- mahalanobis(X, center = colMeans(X), cov = cov(X))
#' # in this case the distances should be the same
#' plot(icsX.dist.all, maha)
#' all.equal(icsX.dist.all, maha)
#' 
#' icsX.dist.first <- ics_distances(icsX, index = 1)
#' plot(icsX.dist.first)
ics_distances <-
    function(object, index = NULL)
    {
        if (!inherits(object, "ICS")) stop("'object' must be of class 'ICS'" )
        if(is.null(index)) index <- 1:ncol(object$scores)
        dists <- rowSums(object$scores[, index, drop = FALSE]^2)
        dists
    }