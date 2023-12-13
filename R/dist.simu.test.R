#' @export
dist.simu.test <-
    function(object,
             index,
             m = 10000,
             level = 0.025,
             ncores = NULL,
             iseed = NULL,
             pkg = "ICSOutlier",
             qtype = 7,
             ...)
    {
        LEVEL <- 1 - level
        nameObject <- deparse(substitute(object))
        if (!inherits(object, "ics2"))
            stop("'object' must be of class ics2")
        S1 <- get(object@S1name)
        S2 <- get(object@S2name)
        if (!is.function(S1))
            stop(paste("S1 in '", nameObject, ", must be a specified as a function"))
        if (!is.function(S2))
            stop(paste("S2 in '", nameObject, ", must be a specified as a function"))
        n <- nrow(object@Scores)
        p <- ncol(object@Scores)
        MEAN <- rep(0, p)
        
        if (!is.null(ncores) && ncores > 1) {
            DOTS <- list(...)
            if (is.null(DOTS$na.rm)) {
                Qna.rm <- FALSE
            } else {
                Qna.rm <- DOTS$na.rm
            }
            if (is.null(DOTS$names)) {
                Qnames <- 7
            } else {
                Qnames <- DOTS$names
            }
            
            if (is.null(iseed)) {
                if (exists(".Random.seed", envir = globalenv())) {
                    oldseed <- get(".Random.seed", envir = globalenv())
                    rm(.Random.seed, envir = globalenv())
                    on.exit(assign(".Random.seed", oldseed, envir = globalenv()))
                }
            }
            
            type <-  "PSOCK"
            cl <- makeCluster(ncores, type = type)
            clusterExport(
                cl,
                c(
                    "n",
                    "m",
                    "MEAN",
                    "S1",
                    "S2",
                    "object",
                    "index",
                    "LEVEL",
                    "qtype",
                    "Qna.rm",
                    "Qnames",
                    "pkg",
                    "iseed"
                ),
                envir = environment()
            )
            clusterEvalQ(cl, lapply(pkg, require, character.only = TRUE))
            clusterSetRNGStream(cl = cl, iseed = iseed)
            DISTS <- parSapply(cl, 1:m, function(i, ...) {
                quantile(
                    dist.simu.test.internal(
                        rmvnorm(n,
                                mean = MEAN),
                        S1 = S1,
                        S2 = S2,
                        S1args = object@S1args,
                        S2args = object@S2args,
                        index = index
                    ),
                    probs = LEVEL,
                    type = qtype,
                    na.rm = Qna.rm,
                    names = Qnames
                )
            })
            stopCluster(cl)
            
        } else {
            DISTS <- replicate(m,
                               quantile(
                                   dist.simu.test.internal(
                                       rmvnorm(n,
                                               mean = MEAN),
                                       S1 = S1,
                                       S2 = S2,
                                       S1args = object@S1args,
                                       S2args = object@S2args,
                                       index = index
                                   ),
                                   probs = LEVEL,
                                   ...
                               ))
        }
        
        if (length(LEVEL) > 1) {
            RES <- rowMeans(DISTS)
        } else {
            RES = mean(DISTS)
            names(RES) = paste0(LEVEL * 100, "%")
        }
        return(RES)
    }


#' Cut-Off Values Using Simulations for the Detection of Extreme ICS Distances
#' 
#' Computes the cut-off values for the identification of the outliers based on the squared ICS distances. It uses simulations under a multivariate standard normal model for a specific data setup and scatters combination.
#'
#' @param object object of class \code{"ICS"} where both \code{S1} and \code{S2} are specified as functions. 
#' The sample size and the dimension of interest are also obtained from the object.
#' The invariant coordinate are required to be centered.
#' @param S1 an object of class \code{"ICS_scatter"} or a function that 
#' contains the location vector and scatter matrix as \code{location} and \code{scatter} components.
#' @param S2 an object of class \code{"ICS_scatter"} or a function that 
#' contains the location vector and scatter matrix as \code{location} and \code{scatter} components.
#' @param S1_args a list containing additional arguments for \code{S1}. 
#' @param S2_args a list containing additional arguments for \code{S2}. 
#' @param index integer vector specifying which components are used to compute the  [ics_distances()].
#' @param m number of simulations. Note that extreme quantiles are of interest and hence \code{m} should be large. 
#' @param level the (1-\code{level}(s))th quantile(s) used to choose the cut-off value(s). Usually just one number between 0 and 1. However a vector is also possible.
#' @param n_cores number of cores to be used. If \code{NULL} or 1, no parallel computing is used. Otherwise \link[parallel]{makeCluster} with \code{type = "PSOCK"} is used.
#' @param iseed If parallel computation is used the seed passed on to \code{\link[parallel]{clusterSetRNGStream}}. Default is NULL which means no fixed seed is used.
#' @param pkg When using parallel computing, a character vector listing all the packages which need to be loaded on the different cores via \code{\link{require}}. Must be at least "ICSOutlier" and must contain the packages needed to compute the scatter matrices.
#' @param q_type specifies the quantile algorithm used in \code{\link{quantile}}.
#' @param ... further arguments passed on to the function \code{\link{quantile}}.
#'
#'
#' @details 
#' The function extracts basically the dimension of the data from the \code{"ICS"} object and simulates \code{m} times, from a multivariate standard normal distribution, the squared ICS distances with the components specified in \code{index}. The resulting value is then the mean of the \code{m} correponding quantiles of these distances  at level 1-\code{level}.
#' 
#' Note that depending on the data size and scatters used this can take a while and so it is more efficient to parallelize computations.
#' 
#' Note that the function is seldomly called directly by the user but internally by [ICS_outlier()].
#' @return A vector with the values of the (1-\code{level})th quantile.
#' @export
#' 
#' @references
#' Archimbaud, A., Nordhausen, K. and Ruiz-Gazen, A. (2018), ICS for multivariate outlier detection with application to quality control. Computational Statistics & Data Analysis, 128:184-199. ISSN 0167-9473.  \doi{10.1016/j.csda.2018.06.011}. 

#' @author Aurore Archimbaud and Klaus Nordhausen
#' 
#' @seealso [ICS()][ICS::ICS()], [ics_distances()]
#' 
#' @import parallel
#' @importFrom mvtnorm rmvnorm
#'
#' @examples
#' # For a real analysis use larger values for m and more cores if available
#'
#' Z <- rmvnorm(1000, rep(0, 6))
#' Z[1:20, 1] <- Z[1:20, 1] + 10
#' A <- matrix(rnorm(36), ncol = 6)
#' X <- tcrossprod(Z, A)
#'
#' pairs(X)
#' icsX <- ICS(X, center = TRUE)
#'
#' icsX.dist.1 <- ics_distances(icsX, index = 1)
#' CutOff <- dist_simu_test(icsX, S1 = ICS_cov, S2= ICS_cov4,
#'                         index = 1, m = 500, ncores = 1)
#'
#' # check if outliers are above the cut-off value
#' plot(icsX.dist.1, col = rep(2:1, c(20, 980)))
#' abline(h = CutOff)
#' 
#' 
#' library(REPPlab)
#' data(ReliabilityData)
#' # The observations 414 and 512 are suspected to be outliers
#' icsReliability <- ICS(ReliabilityData, center = TRUE)
#' # Choice of the number of components with the screeplot: 2
#' screeplot(icsReliability)
#' # Computation of the distances with the first 2 components
#' ics.dist.scree <- ics_distances(icsReliability, index = 1:2)
#' # Computation of the cut-off of the distances 
#' CutOff <- dist_simu_test(icsReliability, S1 = ICS_cov, S2= ICS_cov4,
#'                          index = 1:2, m = 50, level = 0.02, ncores = 1)
#' # Identification of the outliers based on the cut-off value
#' plot(ics.dist.scree)
#' abline(h = CutOff)
#' outliers <- which(ics.dist.scree >= CutOff)
#' text(outliers, ics.dist.scree[outliers], outliers, pos = 2, cex = 0.9)
#' 
#' \dontrun{
#'     # For using three cores
#'     # For demo purpose only small m value, should select the first #' component
#'     dist_simu_test(icsReliability, S1 = ICS_cov, S2= ICS_cov4,
#'                    index = 1:2, m = 500, level = 0.02, n_cores = 3, iseed #' = 123)
#'     
#'     # For using several cores and for using a scatter function from a different package
#'     # Using the parallel package to detect automatically the number of cores
#'     library(parallel)
#'     # ICS with Cauchy estimates
#'     library(ICSClust)
#'     icsReliabilityMLC <- ICS(ReliabilityData, S1 = ICS_mlc, 
#'                             S1_args = list(location = TRUE),
#'                                  S2 = ICS_cov, center = TRUE)
#'     # Computation of the cut-off of the distances. For demo purpose only small m value.
#'     dist_simu_test(icsReliabilityMLC,  S1 = ICS_mlc,  S1_args = list(location = TRUE), 
#'     S2 = ICS_cov, index = 1:2, m = 500, level = 0.02, 
#'     n_cores = detectCores()-1,  pkg = c("ICSOutlier","ICSClust"), iseed = 123)
#' }
dist_simu_test <- function(object,
                           S1 = NULL,
                           S2 = NULL,
                           S1_args = list(),
                           S2_args = list(),
                           index,
                           m = 10000,
                           level = 0.025,
                           n_cores = NULL,
                           iseed = NULL,
                           pkg = "ICSOutlier",
                           q_type = 7,
                           ...)
{
    level <- 1 - level
    name_object <- deparse(substitute(object))
    if (!inherits(object, "ICS"))
        stop("'object' must be of class 'ICS'")
    if (isFALSE(object$center)) 
        stop("'ICS' object must has been computed with `center` = TRUE.")
    
    # S1 <- get(object$S1_label)
    # S2 <- get(object$S2_label)
    if (!is.function(S1))
        stop(paste("S1 must be a specified as a function"))
    if (!is.function(S2))
        stop(paste("S2 must be a specified as a function"))
    n <- nrow(object$scores)
    p <- ncol(object$scores)
    mean_vec <- rep(0, p)
    
    if (!is.null(n_cores) && n_cores > 1) {
        list_dots <- list(...)
        if (is.null(list_dots$na.rm)) {
            Qna_rm <- FALSE
        } else {
            Qna_rm <- list_dots$na.rm
        }
        if (is.null(list_dots$names)) {
            Q_names <- 7
        } else {
            Q_names <- list_dots$names
        }
        
        if (is.null(iseed)) {
            if (exists(".Random.seed", envir = globalenv())) {
                oldseed <- get(".Random.seed", envir = globalenv())
                rm(.Random.seed, envir = globalenv())
                on.exit(assign(".Random.seed", oldseed, envir = globalenv()))
            }
        }
        
        type <-  "PSOCK"
        cl <- makeCluster(n_cores, type = type)
        clusterExport(
            cl,
            c(
                "n",
                "m",
                "mean_vec",
                "S1",
                "S2",
                "S1_args",
                "S2_args",
                "index",
                "level",
                "q_type",
                "Qna_rm",
                "Q_names",
                "pkg",
                "iseed"
            ),
            envir = environment()
        )
        clusterEvalQ(cl, lapply(pkg, require, character.only = TRUE))
        clusterSetRNGStream(cl = cl, iseed = iseed)
        dists <- parSapply(cl, 1:m, function(i, ...) {
            quantile(
                dist_simu_test_internal(
                    rmvnorm(n,
                            mean = mean_vec),
                    S1 = S1,
                    S2 = S2,
                    S1_args = S1_args,
                    S2_args = S2_args,
                    index = index
                ),
                probs = level,
                type = q_type,
                na_rm = Qna_rm,
                names = Q_names
            )
        })
        stopCluster(cl)
        
    } else {
        dists <- replicate(m,
                           quantile(
                               dist_simu_test_internal(
                                   rmvnorm(n,
                                           mean = mean_vec),
                                   S1 = S1,
                                   S2 = S2,
                                   S1_args = S1_args,
                                   S2_args = S2_args,
                                   index = index
                               ),
                               probs = level,
                               ...
                           ))
    }
    
    if (length(level) > 1) {
        res <- rowMeans(dists)
    } else {
        res = mean(dists)
        names(res) = paste0(level * 100, "%")
    }
    res
}
