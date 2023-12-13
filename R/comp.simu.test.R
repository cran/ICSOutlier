#' @export
comp.simu.test <- function(object, m = 10000, type = "smallprop",
                           level = 0.05, adjust = TRUE, ncores = NULL, iseed = NULL, pkg = "ICSOutlier", qtype = 7, ...)
{
  if (!inherits(object, "ics2"))
    stop("'object' must be of class ics2")
  S1 <- get(object@S1name)
  S2 <- get(object@S2name)
  if (!is.function(S1))
    stop(paste("S1 in '", S1, ", must be a specified as a function"))
  if (!is.function(S2))
    stop(paste("S2 in '", S2, ", must be a specified as a function"))
  type <- match.arg(type, c("smallprop"))
  n <- nrow(object@Scores)
  p <- ncol(object@Scores)
  MEAN <- rep(0, p)
  
  if(!is.null(ncores) && ncores > 1){
    
    
    if(is.null(iseed)){
      if (exists(".Random.seed", envir=globalenv())){
        oldseed <- get(".Random.seed", envir=globalenv())
        rm(.Random.seed, envir=globalenv())
        on.exit(assign(".Random.seed", oldseed, envir=globalenv()))
      }
    }
    
    ctype <-  "PSOCK"
    cl <- makeCluster(ncores, type=ctype)
    clusterExport(cl, c("n","m","MEAN","S1","S2","object", "pkg", "iseed"), envir = environment())
    clusterEvalQ(cl, lapply(pkg, require,character.only = TRUE))
    clusterSetRNGStream(cl = cl, iseed = iseed)
    
    EV <- parSapply(cl, 1:m, function(i,...) {
      ics2(rmvnorm(n, MEAN), S1 = S1, S2 = S2,
           S1args = object@S1args, S2args = object@S2args)@gKurt } )
    stopCluster(cl)
    
  } else {
    
    EV <- replicate(m, ics2(rmvnorm(n, MEAN), S1 = S1, S2 = S2,
                            S1args = object@S1args, S2args = object@S2args)@gKurt)
  }
  
  if (adjust == TRUE) {
    levels <- level/1:p
  }
  else {
    levels <- rep(level, p)
  }
  EV.quantile <- numeric(p)
  for (i in 1:p) {
    EV.quantile[i] <- quantile(EV[i, ], probs = 1 - levels[i], type = qtype,
                               ...)
  }
  decisions <- (object@gKurt > EV.quantile)
  k <- match(FALSE, decisions) - 1
  if (is.na(k)) {
    index <- 1:p
  }
  else {
    if (k == 0)
      index <- 0
    else index <- 1:k
  }
  RES <- list(index = index, test = "simulation", criterion = EV.quantile,
              levels = levels, adjust = adjust, type = type, m = m)
  RES
}

#' Selection of Nonnormal Invariant Components Using Simulations
#' 
#' Identifies invariant coordinates that are nonnormal using simulations under a standard multivariate normal model for a specific data setup and scatter combination.
#'
#' @param object object of class \code{"ICS"} where both \code{S1} and \code{S2} are specified as functions. 
#' The sample size and the dimension of interest are also obtained from the object.
#' It is also natural to expect that the invariant coordinate are centered.
#' @param S1 an object of class \code{"ICS_scatter"} or a function that 
#' contains the location vector and scatter matrix as \code{location} and \code{scatter} components.
#' @param S2 an object of class \code{"ICS_scatter"} or a function that 
#' contains the location vector and scatter matrix as \code{location} and \code{scatter} components.
#' @param S1_args a list containing additional arguments for \code{S1}. 
#' @param S2_args a list containing additional arguments for \code{S2}. 
#' @param m number of simulations. Note that since extreme quantiles are of interest \code{m} should be large.
#' @param type currently the only type option is \code{"smallprop"}. See details.
#' @param level the initial level used to make a decision. The cut-off values are the (1-\code{level})th quantile of the eigenvalues obtained from simulations. See details.
#' @param adjust logical. If \code{TRUE}, the quantiles levels are adjusted. Default is \code{TRUE}. See details.
#' @param n_cores number of cores to be used. If \code{NULL} or 1, no parallel computing is used. Otherwise \link[parallel]{makeCluster} with \code{type = "PSOCK"} is used.
#' @param iseed If parallel computation is used the seed passed on to \code{\link[parallel]{clusterSetRNGStream}}. Default is NULL which means no fixed seed is used.
#' @param pkg When using parallel computing, a character vector listing all the packages which need to be loaded on the different cores via \code{\link{require}}. Must be at least "ICSOutlier" and must contain the packages needed to compute the scatter matrices.
#' @param q_type specifies the quantile algorithm used in \code{\link{quantile}}.
#' @param ... further arguments passed on to the function \code{\link{quantile}}.
#'
#' @details 
#' Based on simulations it detects which of the components follow a univariately normal distribution. More precisely it identifies the observed eigenvalues larger than the ones coming
#' from normal distributed data. \code{m} standard normal data sets are simulated using the same data size and scatters as specified in the \code{"ICS"} object.
#' The cut-off values are determined based on a quantile of these simulated eigenvalues. 
#' 
#' 
#' As the eigenvalues, aka generalized kurtosis values, of ICS are ordered it is natural to perform the comparison in a specific order depending on the purpose. 
#' Currently the only available \code{type} is \code{"smallprop"} so starting with the first component, the observed eigenvalues are successively compared to these cut-off values. The precedure stops when an eigenvalue is below the corresponding cut-off, so when a normal component is detected. 
#' 
#' If \code{adjust = FALSE} all eigenvalues are compared to the same (\code{1-level})th level of the quantile. This leads however often to too many selected components. 
#' Therefore some multiple testing adjustment might be useful. The current default adjusts the quantile for the jth component as \code{1-level}/j.
#' 
#' Note that depending on the data size and scatters used this can take a while and so it is more efficient to parallelize computations. 
#' Note also that the function is seldomly called directly by the user but internally by [ICS_outlier()].

#' @return A list containing:
#' - `index`:  integer vector indicating the indices of the selected components.
#' - `test`: string \code{"simulation"}.
#' - `criterion`: vector of the cut-off values for all the eigenvalues.
#' - `levels`: vector of the levels used for the decision for each component.
#' - `adjust`: logical. \code{TRUE} if adjusted.
#' - `type`: \code{type} used
#' - `m`: number of iterations \code{m} used in the simulations.
#' 
#' 
#' @export
#' 
#' @author Aurore Archimbaud and Klaus Nordhausen
#' 
#' @references
#' Archimbaud, A., Nordhausen, K. and Ruiz-Gazen, A. (2018), ICS for multivariate outlier detection with application to quality control. Computational Statistics & Data Analysis, 128:184-199. ISSN 0167-9473.  \doi{10.1016/j.csda.2018.06.011}. 
#' 
#' @seealso [ICS()][ICS::ICS()], [comp_norm_test()]
#' @import parallel
#' @import ICS
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats quantile
#'
#' @examples
#' # For a real analysis use larger values for m and more cores if available
#' set.seed(123)
#' Z <- rmvnorm(1000, rep(0, 6))
#' # Add 20 outliers on the first component
#' Z[1:20, 1] <- Z[1:20, 1] + 10
#' pairs(Z)
#' icsZ <- ICS(Z)
#' # For demo purpose only small m value, should select the first component
#' comp_simu_test(icsZ, S1 = ICS_cov, S2= ICS_cov4, m = 400, n_cores = 1)
#' 
#' \dontrun{
#'  # For using two cores
#'   # For demo purpose only small m value, should select the first component
#'   comp_simu_test(icsZ, S1 = ICS_cov, S2 = ICS_cov4, m = 500, n_cores = 2, iseed = 123)
#'   # For using several cores and for using a scatter function from a different package
#'   # Using the parallel package to detect automatically the number of cores
#'   library(parallel)
#'   # ICS with MCD estimates and the usual estimates
#'   library(ICSClust)
#'         icsZmcd <- ICS(Z, S1 = ICS_mcd_raw, S2 = ICS_cov, S1_args = list(alpha = 0.75))
#'         # For demo purpose only small m value, should select the first component
#'         comp_simu_test(icsZmcd, S1 = ICS_mcd_raw, S2 = ICS_cov, 
#'         S1_args = list(alpha = 0.75, location = TRUE),
#'          m = 500, ncores = detectCores()-1, 
#'                     pkg = c("ICSOutlier", "ICSClust"), iseed = 123)
#'  }
#'  # Example with no outlier
#'  Z0 <- rmvnorm(1000, rep(0, 6))
#'  pairs(Z0)
#'  icsZ0 <- ICS(Z0)
#'  # Should select no component
#'  comp_simu_test(icsZ0,S1 = ICS_cov, S2 = ICS_cov4, m = 400, level = 0.01, n_cores = 1)
comp_simu_test <- function(object,
                           S1 = NULL,
                           S2 = NULL,
                           S1_args = list(),
                           S2_args = list(), m = 10000, type = "smallprop",
                           level = 0.05, adjust = TRUE, n_cores = NULL, 
                           iseed = NULL, pkg = "ICSOutlier", q_type = 7, ...)
{
  if (!inherits(object, "ICS"))
    stop("'object' must be of class 'ICS'")
  
  # S1 <- get(object$S1_label)
  # S2 <- get(object$S2_label)
  if (!is.function(S1))
    stop(paste("S1 must be specified as a function"))
  if (!is.function(S2))
    stop(paste("S2 must be specified as a function"))
  type <- match.arg(type, c("smallprop"))
  n <- nrow(object$scores)
  p <- ncol(object$scores)
  mean_vec <- rep(0, p)
  
  if (!is.null(n_cores) && n_cores > 1) {
    if (is.null(iseed)) {
      if (exists(".Random.seed", envir = globalenv())) {
        oldseed <- get(".Random.seed", envir = globalenv())
        rm(.Random.seed, envir = globalenv())
        on.exit(assign(".Random.seed", oldseed, envir = globalenv()))
      }
    }
    
    ctype <-  "PSOCK"
    cl <- makeCluster(n_cores, type = ctype)
    clusterExport(cl,
                  c("n", "m", "mean_vec", "S1", "S2", 
                    "S1_args", "S2_args",
                    "pkg", "iseed"),
                  envir = environment())
    clusterEvalQ(cl, lapply(pkg, require, character.only = TRUE))
    clusterSetRNGStream(cl = cl, iseed = iseed)
    
    EV <- parSapply(cl, 1:m, function(i, ...) {
      tryCatch({
        ICS(
          rmvnorm(n, mean_vec),
          S1 = S1,
          S2 = S2,
          S1_args = S1_args,
          S2_args = S2_args,
          center = TRUE,
          fix_signs = "scores"
        )$gen_kurtosis
      },
      warning = function(w) stop(w),
      error = function(e) stop(e))
    })
    stopCluster(cl)
    
  } else {
    EV <- replicate(
      m,
      tryCatch({
        ICS(
          rmvnorm(n, mean_vec),
          S1 = S1,
          S2 = S2,
          S1_args = S1_args,
          S2_args = S2_args,
          center = TRUE,
          fix_signs = "scores"
        )$gen_kurtosis
      },
      warning = function(w) stop(w),
      error = function(e) stop(e))
    )
  }
  
  if (adjust == TRUE) {
    levels <- level / 1:p
  }
  else {
    levels <- rep(level, p)
  }
  EV.quantile <- numeric(p)
  for (i in 1:p) {
    EV.quantile[i] <-
      quantile(EV[i,], probs = 1 - levels[i], type = q_type,
               ...)
  }
  decisions <- (object$gen_kurtosis > EV.quantile)
  k <- match(FALSE, decisions) - 1
  if (is.na(k)) {
    index <- 1:p
  }
  else {
    if (k == 0)
      index <- 0
    else
      index <- 1:k
  }
  res <-
    list(
      index = index,
      test = "simulation",
      criterion = EV.quantile,
      levels = levels,
      adjust = adjust,
      type = type,
      m = m
    )
  res
}

