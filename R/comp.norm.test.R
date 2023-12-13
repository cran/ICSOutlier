
#' @export
comp.norm.test <-
function(object, test = "agostino.test", type = "smallprop", level = 0.05, adjust = TRUE)
    {
    if (!inherits(object, "ics2")) stop("'object' must be of class ics2")
    test <- match.arg(test, c("jarque.test", "anscombe.test", "bonett.test", "agostino.test", "shapiro.test"))
    #type <- match.arg(type, c("smallprop", "cluster", "all"))
    type <- match.arg(type, c("smallprop"))
    
    test.res <- apply(object@Scores, 2, test)
    test.pvals <- sapply(test.res, function(x) x$p.value)
    #test.pvals
    p <- ncol(object@Scores)
    if (adjust == TRUE) levels <- level/1:p else levels <- rep(level, p)
    decisions <- (test.pvals <= levels)
    k <- match(FALSE, decisions)-1
    if (is.na(k)) { index <- 1:p } else { if (k == 0) index <- 0 else index <- 1:k}
    RES <- list(index = index, test = test, criterion = test.pvals, levels = levels, adjust = adjust, type = type)
    RES
    }


#' Selection of Nonnormal Invariant Components Using Marginal Normality Tests
#'
#'
#' Identifies invariant coordinates that are non normal using univariate normality tests.
#' 
#' @param object object of class \code{"ICS"} where both \code{S1} and \code{S2} are specified as functions. 
#' The sample size and the dimension of interest are also obtained from the object.
#' @param test name of the normality test to be used. Possibilites are \code{"jarque.test"}, 
#' \code{"anscombe.test"}, \code{"bonett.test"}, \code{"agostino.test"}, \code{"shapiro.test"}.
#' Default is \code{"agostino.test"}.
#' @param type currently the only option is \code{"smallprop"}. See details.
#' @param level the initial level used to make a decision based on the test p-values. See details.
#' @param adjust logical. If \code{TRUE}, the quantiles levels are adjusted. Default is \code{TRUE}. See details.
#' 
#' 
#' @details Currently the only available \code{type} is \code{"smallprop"} which detects which of the components follow a univariately normal distribution. It starts from the first component and stops when a component is detected as gaussian. Five tests for univariate normality are available. See [normal_crit()][ICSClust::normal_crit()] function for more general cases. 
#' 
#' If \code{adjust = FALSE} all tests are performed at the same \code{level}. This leads however often to too many components. Therefore some multiple testing adjustments might be useful. The current default adjusts the level for the jth component as \code{level}/j. 
#' 
#' Note that the function is seldomly called directly by the user but internally by [ICS_outlier()]. 
#'
#' @return A list containing:
#' - `index`:  integer vector indicating the indices of the selected components.
#' - `test`: string with the name of the normality test used.
#' - `criterion`: vector of the p-values from the marginal normality tests for each component.
#' - `levels`: vector of the levels used for the decision for each component.
#' - `adjust`: logical. \code{TRUE} if adjusted.
#' - `type`: \code{type} used
#' 
#' @export
#' 
#' @references
#' Archimbaud, A., Nordhausen, K. and Ruiz-Gazen, A. (2018), ICS for multivariate outlier detection with application to quality control. Computational Statistics & Data Analysis, 128:184-199. ISSN 0167-9473.  \doi{10.1016/j.csda.2018.06.011}. 

#' @author Aurore Archimbaud and Klaus Nordhausen
#' @import moments
#' @importFrom stats shapiro.test
#' 
#' @seealso [ICS()][ICS::ICS()], [comp_simu_test()], [jarque.test()][moments::jarque.test()], 
#' [anscombe.test()][moments::anscombe.test()], [bonett.test()][moments::bonett.test()], [bonett.test()][moments::agostino.test()], 
#'  [shapiro.test()][stats::shapiro.test()]
#'
#' @examples
#' 
#' Z <- rmvnorm(1000, rep(0, 6))
#' # Add 20 outliers on the first component
#' Z[1:20, 1] <- Z[1:20, 1] + 10
#' pairs(Z)
#' icsZ <- ICS(Z)
#' # The shift located outliers can be displayed in one dimension
#' comp_norm_test(icsZ)
#' # Only one invariant component is non normal and selected.
#' comp_norm_test(icsZ, test = "bonett.test")
#' 
#' # Example with no outlier
#' Z0 <- rmvnorm(1000, rep(0, 6))
#' pairs(Z0)
#' icsZ0 <-ICS(Z0)
#' # Should select no component
#' comp_norm_test(icsZ0, level = 0.01)$index
#' 
comp_norm_test <-
    function(object, test = "agostino.test", type = "smallprop", level = 0.05, adjust = TRUE)
    {
        if (!inherits(object, "ICS")) stop("'object' must be of class 'ICS'")
        test <- match.arg(test, c("jarque.test", "anscombe.test", "bonett.test", "agostino.test", "shapiro.test"))
        #type <- match.arg(type, c("smallprop", "cluster", "all"))
        type <- match.arg(type, c("smallprop"))
        
        test_res <- apply(object$scores, 2, test)
        test_pvals <- sapply(test_res, function(x) x$p.value)
        p <- ncol(object$scores)
        if (adjust == TRUE) levels <- level/1:p else levels <- rep(level, p)
        decisions <- (test_pvals <= levels)
        k <- match(FALSE, decisions)-1
        if (is.na(k)) { index <- 1:p } else { if (k == 0) index <- 0 else index <- 1:k}
        res <- list(index = index, test = test, criterion = test_pvals, 
                    levels = levels, adjust = adjust, type = type)
        res
    }
