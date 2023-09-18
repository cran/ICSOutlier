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
