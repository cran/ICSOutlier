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
