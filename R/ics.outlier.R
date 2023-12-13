#' @export
ics.outlier <-
  function(object, method = "norm.test", test = "agostino.test", mEig = 10000, level.test = 0.05, adjust = TRUE, 
           level.dist = 0.025, mDist = 10000, type = "smallprop", ncores = NULL, iseed = NULL, pkg = "ICSOutlier", qtype = 7, ...)
  {
    # choose method - interpolation should be added when available
    METHOD <- match.arg(method, c("norm.test", "simulation"))
    
    if (!inherits(object, "ics2")) stop("'object' must be of class ics2")
    
    S1 <- get(object@S1name)
    S2 <- get(object@S2name)
    if (!is.function(S1)) stop(paste("S1 in '", S1, ", must be a specified as a function"))
    if (!is.function(S2)) stop(paste("S2 in '", S2, ", must be a specified as a function"))
    
    ROWNAMES <- rownames(object@Scores)
    
    #type <- match.arg(type, c("smallprop", "cluster", "all"))
    type <- match.arg(type, c("smallprop"))
    
    ResMethod <- switch(METHOD, norm.test = {comp.norm.test(object, test = test, level = level.test, adjust = adjust, type = type)},
                        simulation = {comp.simu.test (object, m = mEig, level = level.test, adjust = adjust, type = type, ncores = ncores, iseed = iseed, pkg = pkg, qtype = qtype, ...)}
    )
    
    n <- nrow(object@Scores)
    p <- ncol(object@Scores)
    
    if(sum(ResMethod$index<0.5)) {
      outliers <- rep(0L, n)
      names(outliers) <- ROWNAMES
      ICdistances <- rep(0L, n)
      names(ICdistances) <- ROWNAMES
      ICdistancesQuantile <- rep(0, n)
    } else {
      ICdistancesQuantile <- dist.simu.test(object, m = mDist, index = ResMethod$index, level = level.dist, ncores = ncores, iseed = iseed, pkg = pkg, qtype = qtype, ...)
      ICdistances <- ics.distances(object, index = ResMethod$index)
      outliers <- as.integer(ICdistances > ICdistancesQuantile)
      names(outliers) <- ROWNAMES
    }
    #RES <- list(outliers = outliers, ICdistances = ICdistances, ICdistancesQuantile = ICdistancesQuantile, quant = quant, ResMethod, ICS2 = object)
    
    
    RES <-new("icsOut", outliers = outliers, 
              ics.distances = ICdistances, 
              ics.dist.cutoff = ICdistancesQuantile, 
              level.dist = level.dist, 
              level.test = level.test, 
              method = METHOD,
              index = ResMethod$index,
              test =  ResMethod$test,
              criterion = ResMethod$criterion,
              adjust = ResMethod$adjust,
              type =  ResMethod$type,
              mDist = as.integer(mDist),
              mEig = as.integer(mEig),
              S1name = object@S1name, 
              S2name = object@S2name)
    
    RES
  }



#' Outlier Detection Using ICS
#' 
#' In a multivariate framework outlier(s) are detected using ICS. The function performs [ICS()][ICS::ICS()] and decides automatically about the number of invariant components to use to search for the outliers and the number of outliers detected on these components. Currently the function is restricted to the case of searching outliers only on the first components. 
#'
#'
#' @param X a numeric matrix or data frame containing the data to be
#' transformed.
#' @param S1 an object of class \code{"ICS_scatter"} or a function that 
#' contains the location vector and scatter matrix as \code{location} and \code{scatter} components.
#' @param S2 an object of class \code{"ICS_scatter"} or a function that 
#' contains the location vector and scatter matrix as \code{location} and \code{scatter} components.
#' @param S1_args a list containing additional arguments for \code{S1}. 
#' @param S2_args a list containing additional arguments for \code{S2}. 
#' @param ICS_algorithm a character string specifying with which algorithm
#' the invariant coordinate system is computed. Possible values are
#'  \code{"whiten"}, \code{"standard"} or \code{"QR"}.
#' @param method name of the method used to select the ICS components involved to compute ICS distances. Options are \code{"norm_test"} and \code{"simulation"}. Depending on the method either \code{\link{comp_norm_test}} or \code{\link{comp_simu_test}} are used.
#' @param test name of the marginal normality test to use if \code{method = "norm_test"}. Possibilities are \code{"jarque.test"}, \code{"anscombe.test"}, \code{"bonett.test"}, \code{"agostino.test"}, \code{"shapiro.test"}.Default is \code{"agostino.test"}.
#' @param n_eig number of simulations performed to derive the cut-off values for selecting the ICS components. Only if \code{method = "simulation"}. See \code{\link{comp_simu_test}} for details.
#' @param level_test  for the \code{\link{comp_norm_test}} or \code{\link{comp_simu_test}} functions. The initial level for selecting the invariant coordinates.
#' @param adjust logical. For selecting the invariant coordinates, the level of the test can be adjusted for each component to deal with multiple testing. See \code{\link{comp_norm_test}} and \code{\link{comp_simu_test}} for details. Default is TRUE.
#' @param level_dist \code{level} for the \code{\link{dist_simu_test}} function.  The (1-\code{level})th quantile used to determine the cut-off value for the ICS distances.
#' @param n_dist number of simulations performed to derive the cut-off value for the ICS distances. See \code{\link{dist_simu_test}} for details.
#' @param type currently the only option is \code{"smallprop"} which means that only the first ICS components can be selected. See  \code{\link{comp_norm_test}} or \code{\link{comp_simu_test}} for details.
#' @param n_cores number of cores to be used in \code{\link{dist_simu_test}} and \code{\link{comp_simu_test}}. If \code{NULL} or 1, no parallel computing is used. Otherwise \link[parallel]{makeCluster} with \code{type = "PSOCK"} is used.
#' @param iseed If parallel computation is used the seed passed on to \code{\link[parallel]{clusterSetRNGStream}}. Default is NULL which means no fixed seed is used.
#' @param pkg When using parallel computing, a character vector listing all the packages which need to be loaded on the different cores via \code{\link{require}}. Must be at least "ICSOutlier" and must contain the packages needed to compute the scatter matrices.
#' @param q_type specifies the quantile algorithm used in \code{\link{quantile}}.
#' @param ... passed on to other methods.
#' 
#' @details
#' The ICS method has attractive properties for outlier detection in the case of a small proportion of outliers. As for PCA three steps have to be performed:(i) select the components most useful for the detection, (ii) compute distances as outlierness measures for all observation and finally (iii) label outliers using some cut-off value.
#' 
#' This function performs these three steps automatically:
#' 
#' - For choosing the components of interest two methods are proposed: \code{"norm_test"} based on some marginal normality tests (see details in \code{\link{comp_norm_test}}) or \code{"simulation"} based on a parallel analysis (see details in \code{\link{comp_simu_test}}). These two approaches lie on the intrinsic property of ICS in case of a small proportion of outliers with the choice of S1 "more robust" than S2, which ensures to find outliers on the first components. Indeed  when using \code{S1 = ICS_cov} and \code{S2 = ICS_cov4}, the Invariant Coordinates are ordered according to their classical Pearson kurtosis values in decreasing order. The information to find the outliers should be then contained in the first k non-normal directions.
#' 
#' - Then the ICS distances are computed as the Euclidean distances on the selected k centered components \eqn{Z_k}. 
#' 
#' - Finally the outliers are identified based on a cut-off derived from simulations. If the distance of an observation exceeds the expectation under the normal model, this observation is labeled as outlier (see details in \code{\link{dist_simu_test}}).
#' 
#' As a rule of thumb, the percentage of contamination should be limited to 10% in case of a mixture of gaussian distributions and using the default combination of locations and scatters for ICS.
#' 
#'
#' @return 
#' An object of S3-class 'ICS_Out' which contains:
#' 
#' - `outliers`:  A vector containing ones for outliers and zeros for non outliers.
#' - `ics_distances`:  A numeric vector containing the squared ICS distances.
#' - `ics_dist_cutoff`:  The cut-off for the distances to decide if an observation is outlying or not.
#' - `level_dist`:  The level for deciding upon the cut-off value for the ICS distances.
#' - `level_test`:  The initial level for selecting the invariant coordinates.
#' - `method`: Name of the method used to decide upon the number of ICS components.
#' - `index`: Vector giving the indices of the ICS components selected.
#' - `test`:  The name of the normality test as specified in the function call.
#' - `criterion`:  Vector giving the marginal levels for the components selection.
#' - `adjust`:  Wether the initial level used to decide upon the number of components has been adjusted for multiple testing or not.
#' - `type`:  Currently always the string \code{"smallprop"}.
#' - `n_dist`: Number of simulations performed to decide upon the cut-off for the ICS distances.
#' - `n_eig`:  Number of simulations performed for selecting the ICS components based on simulations.
#' - `S1_label`:  Name of S1.
#' - `S2_label `: Name of S2.
#' @export
#' 
#' @import ICS
#' 
#' @references
#' Archimbaud, A., Nordhausen, K. and Ruiz-Gazen, A. (2018), ICS for multivariate outlier detection with application to quality control. Computational Statistics & Data Analysis, 128:184-199. ISSN 0167-9473.  \doi{10.1016/j.csda.2018.06.011}. 

#' @author Aurore Archimbaud and Klaus Nordhausen
#' 
#' @seealso [ICS()][ICS::ICS()], [comp_norm_test()], [comp_simu_test()],
#' [dist_simu_test()] and 
#' [print()][print.ICS_Out()], [plot()][plot.ICS_Out()], [summary()][summary.ICS_Out()] methods
#'
#' @examples
#' # ReliabilityData example: the observations 414 and 512 are suspected to be outliers  
#' library(REPPlab)
#' data(ReliabilityData)
#' # For demo purpose only small mDist value, but as extreme quantiles
#' # are of interest mDist should be much larger. Also number of cores used
#' # should be larger if available
#' icsOutlierDA <- ICS_outlier(ReliabilityData, S1 = ICS_tM, S2 = ICS_cov, 
#' level_dist = 0.01, n_dist = 50, n_cores = 1)
#' icsOutlierDA
#' summary(icsOutlierDA)
#' plot(icsOutlierDA)
#' 
#' \dontrun{
#'   # For using several cores and for using a scatter function from a different package
#'   # Using the parallel package to detect automatically the number of cores
#'   library(parallel)
#'   # ICS with MCD estimates and the usual estimates
#'   # Need to create a wrapper for the CovMcd function to return first the location estimate
#'   # and the scatter estimate secondly.
#'   data(HTP)
#'  library(ICSClust)
#'   # For demo purpose only small m value, should select the first seven components
#'   icsOutlier <- ICS_outlier(HTP, S1 = ICS_mcd_rwt, S2 = ICS_cov,
#'                             S1_args = list(location = TRUE, alpha = 0.75),
#'                             n_eig = 50, level_test = 0.05, adjust = TRUE,
#'                             level_dist = 0.025, n_dist = 50,
#'                             n_cores =  detectCores()-1, iseed = 123,
#'                             pkg = c("ICSOutlier", "ICSClust"))
#'   icsOutlier
#' }
#' 
#' # Exemple of no direction and hence also no outlier
#' set.seed(123)
#' X = rmvnorm(500, rep(0, 2), diag(rep(0.1,2)))
#' icsOutlierJB <- ICS_outlier(X, test = "jarque.test", level_dist = 0.01,
#'                             level_test = 0.01, n_dist = 100, n_cores = 1)
#' summary(icsOutlierJB)
#' plot(icsOutlierJB)
#' rm(.Random.seed)
#' 
#' # Example of no outlier
#' set.seed(123)
#' X = matrix(rweibull(1000, 4, 4), 500, 2)
#' X = apply(X,2, function(x){ifelse(x<5 & x>2, x, runif(sum(!(x<5 & x>2)), 5, 5.5))})
#' icsOutlierAG <- ICS_outlier(X, test = "anscombe.test", level_dist = 0.01,
#'                             level_test = 0.05, n_dist = 100, n_cores = 1)
#' summary(icsOutlierAG)
#' plot(icsOutlierAG)
#' rm(.Random.seed)
ICS_outlier <- function(X,
                        S1 = ICS_cov,
                        S2 = ICS_cov4,
                        S1_args = list(),
                        S2_args = list(),
                        ICS_algorithm = c("whiten", "standard", "QR"),
                        method = "norm_test", test = "agostino.test", 
                        n_eig = 10000, level_test = 0.05, adjust = TRUE, 
                        level_dist = 0.025, n_dist = 10000, type = "smallprop",
                        n_cores = NULL, iseed = NULL, pkg = "ICSOutlier", q_type = 7, ...)
{
  # match algorithm
  algorithm <- match.arg(ICS_algorithm)
  # choose method - interpolation should be added when available
  method <- match.arg(method, c("norm_test", "simulation"))
  
  #if (!inherits(object, "ICS")) stop("'object' must be of class 'ICS'")
  if (!(inherits(X, "data.frame") | inherits(X, "matrix"))){
    stop("'X' must be of class 'data.frame' or 'matrix'")
  } 
  
  # S1 <- get(object$S1_label)
  # S2 <- get(object$S2_label)
  if (!is.function(S1)) stop(paste("S1 must be specified as a function"))
  if (!is.function(S2)) stop(paste("S2 must be specified as a function"))
  
  # ICS
  object <- tryCatch({
    ICS(X, S1 = S1, S2 = S2, S1_args = S1_args, S2_args = S2_args,
        algorithm = algorithm,
        center = TRUE, fix_signs = "scores")
  },
  warning = function(w) stop(w),
  error = function(e) stop(e))
  rownames(object$scores) <- rownames(X)
  
  row_names <- rownames(object$scores)
  
  #type <- match.arg(type, c("smallprop", "cluster", "all"))
  type <- match.arg(type, c("smallprop"))
  
  res_method <- switch(method,
                       norm_test = {
                         comp_norm_test(
                           object,
                           test = test,
                           level = level_test,
                           adjust = adjust,
                           type = type
                         )
                       },
                       simulation = {
                         comp_simu_test (
                           object,
                           S1 = S1,
                           S2 = S2,
                           S1_args = S1_args,
                           S2_args = S2_args,
                           m = n_eig,
                           level = level_test,
                           adjust = adjust,
                           type = type,
                           n_cores = n_cores,
                           iseed = iseed,
                           pkg = pkg,
                           q_type = q_type,
                           ...
                         )
                       })
  
  n <- nrow(object$scores)
  p <- ncol(object$scores)
  
  if(sum(res_method$index < 0.5)) {
    outliers <- rep(0L, n)
    names(outliers) <- row_names
    IC_distances <- rep(0L, n)
    names(IC_distances) <- row_names
    IC_distances_quantile <- rep(0, n)
  } else {
    IC_distances_quantile <-
      dist_simu_test(
        object,
        S1 = S1,
        S2 = S2,
        S1_args = S1_args,
        S2_args = S2_args,
        m = n_dist,
        index = res_method$index,
        level = level_dist,
        n_cores = n_cores,
        iseed = iseed,
        pkg = pkg,
        q_type = q_type,
        ...
      )
    IC_distances <- ics_distances(object, index = res_method$index)
    outliers <- as.integer(IC_distances > IC_distances_quantile)
    names(outliers) <- row_names
    
  }
  
  res <- list(outliers = outliers, 
              ics_distances = IC_distances, 
              ics_dist_cutoff = IC_distances_quantile, 
              level_dist = level_dist, 
              level_test = level_test, 
              method = method,
              index = res_method$index,
              test = res_method$test,
              criterion = res_method$criterion,
              adjust = res_method$adjust,
              type =  res_method$type,
              n_dist = as.integer(n_dist),
              n_eig = as.integer(n_eig),
              S1_label = object$S1_label, 
              S2_label = object$S2_label)
  
  class(res) <- "ICS_Out"
  res
}


#'  Summary of an 'ICS_Out' Object
#'  
#'  Summarizes an 'ICS_Out' object in an informative way.
#'
#'
#' @param object object object of class `"ICS_Out"`.
#' @param ... additional arguments passed to [summary()]
#'
#'
#' @return An object of class `"ICS_Out_summary"` with the following components:
#' - `comps`: Vector giving the indices of the ICS components selected.
#' - `method`: Name of the method used to decide upon the number of ICS components.
#' - `test`: he name of the normality test as specified in the function call.
#' - `S1_label`:  Name of S1.
#' - `S2_label`:  Name of S2.
#' - `level_test`: The level for deciding upon the cut-off value for the ICS distances.
#' - `level_dist`: The initial level for selecting the invariant coordinates.
#' - `nb_outliers`: the number of observations identified as outliers.
#' 
#' @export
#' 
#' @name summary.ICS_Out
#' @method summary ICS_Out
#' @author Aurore Archimbaud and Klaus Nordhausen
#' 
summary.ICS_Out <- function(object, ...) {
  out <- list(comps = object$index,
              method = object$method,
              test = object$test,
              S1_label = object$S1_label, 
              S2_label = object$S2_label,
              level_test = object$level_test,
              level_dist = object$level_dist,
              nb_outliers = sum(object$outliers)
  )
  
  class(out) <- "ICS_Out_summary"
  out
}

#' Print of an `ICS_Out_summary` object
#'
#' Prints an `ICS_Out_summary` object in an informative way.
#' 
#' @param x object of class `"ICS_Out_summary"`.
#' @param ...  additional arguments, not used.
#'
#' @return The supplied object of class `"ICS_Out_summary"` is returned invisibly.
#' @noRd
#' @name print.ICS_Out_summary
#' @method print ICS_Out_summary
#' @author Aurore Archimbaud and Klaus Nordhausen
print.ICS_Out_summary <-  function(x,  ...) {
  if (sum(x$comps) < 0.5)
    ncomps <- 0
  else
    ncomps <- length(x$comps)
  if (x$method == "norm.test")
    method <- paste(x$method, " (", x$test, ")", sep = "")
  else
    method <- x$method
  cat("\nICS based on two scatter matrices and two location estimates\n")
  cat("S1: ", x$S1_label)
  cat("\nS2: ", x$S2_label)
  cat("\n")
  cat("\nSearching for a small proportion of outliers\n")
  cat("\n")
  #cat(paste("Components selected: ", ncomps, sep = "") )
  #cat("\n")
  cat(paste(
    "Components selected at nominal level ",
    x$level_test,
    ": ",
    ncomps,
    sep = ""
  ))
  cat("\n")
  cat(paste("Selection method: ", x$method, sep = ""))
  cat("\n")
  #cat(paste("Number of outliers: ", sum(object@outliers), sep = ""))
  cat(paste(
    "Number of outliers at nominal level ",
    x$level_dist,
    ": ",
    x$nb_outliers,
    sep = ""
  ))
  cat("\n")
  invisible(x)
}

#' Vector of Outlier Indicators
#' 
#'  Short statement about how many components are selected for the outlier detection and how many outliers are detected.
#'
#' @param x object object of class `"ICS_Out"`.
#' @param ... additional arguments, not used.
#'
#' @return The supplied object of class `"ICS_Out_summary"` is returned invisibly.
#' @export
#' @name print.ICS_Out
#' @method print ICS_Out
#' @author Aurore Archimbaud and Klaus Nordhausen
#'
print.ICS_Out <- function(x, ...)
{
  comps <- x$index
  if (sum(comps) < 0.5) {
    ncomps <- 0
    print(paste(ncomps, " components were selected and no outliers were detected.", sep = ""))
  } else {
    ncomps <- length(comps)
    print(paste(ncomps, " components were selected and ", sum(x$outliers), " outliers were detected.", sep = ""))
  }
  # return x invisibly
  invisible(x)
}


#' Distances Plot for an 'ICS_Out' Object
#' 
#'  Distances plot for an 'ICS_Out' object visualizing the separation of the outliers from the good data points.
#'
#' @param x object of class `"ICS_Out"`.
#' @param pch.out plotting symbol for the outliers.
#' @param pch.good plotting symbol for the 'good' data points.
#' @param col.out color for the outliers.
#' @param col.good color for the 'good' data points.
#' @param col.cut color for cut-off line.
#' @param lwd.cut lwd value for cut-off line.
#' @param lty.cut lty value for cut-off line.
#' @param xlab default x-axis label.
#' @param ylab default y-axis label.
#' @param ... other arguments for \code{plot}
#'
#' @return A plot is displayed.
#' 
#' @details For the figure the IC distances are plotted versus their index. The cut-off value for distances is given as a horizontal line and all observations above the line are considered as outliers.
#' @export
#' 
#' @importFrom grDevices grey
#' @importFrom graphics abline text
#' @importFrom methods hasArg new 
#' @author Aurore Archimbaud and Klaus Nordhausen
#' 
#' @name plot.ICS_Out
#' @method plot ICS_Out
plot.ICS_Out <- function(x, pch.out = 16, pch.good = 4, col.out = 1, col.good = grey(0.5), col.cut = 1, 
                        lwd.cut = 1, lty.cut = 1, xlab = "Observation Number", ylab = "ICS distances", ...)
{
  
  YESylim <- hasArg("ylim") 
  if (sum(x$outliers)>0.5){
    colPoints <- ifelse(x$outliers == 1L, col.out, col.good)
    pchPoints <- ifelse(x$outliers == 1L, pch.out, pch.good)
    plot(x$ics_distances, col = colPoints, pch = pchPoints, xlab = xlab, ylab = ylab, ...)
    abline(h = x$ics_dist_cutoff, col = col.cut, lwd = lwd.cut, lty = lty.cut)
  } else {
    if (sum(x$ics_distances) > 0){
      if (!YESylim) {
        plot(x$ics_distances, col = col.good, pch = pch.good, xlab = xlab, ylab = ylab, ylim = c(0, x$ics_dist_cutoff), ...)
      }  else {
        plot(x$ics_distances, col = col.good, pch = pch.good, xlab = xlab, ylab = ylab, ...)
      }
      abline(h = x$ics_dist_cutoff, col = col.cut, lwd = lwd.cut, lty = lty.cut)
    } else {
      plot(x$ics_distances, xlab = xlab, ylab = ylab, type = "n", ...)
      text(length(x$ics_distances)/2, 0, labels = c("No components have been selected for outlier detection.\n There is nothing to plot."))
    }
  }
}
