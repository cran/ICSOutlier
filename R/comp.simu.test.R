comp.simu.test <- function(object, m = 10000, type = "smallprop",
level = 0.05, adjust = TRUE, ncores = NULL, iseed = NULL, pkg = "ICSOutlier", qtype = 7, ...)
    {
    if (class(object) != "ics2")
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
