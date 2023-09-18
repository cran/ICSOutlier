dist.simu.test <- function(object, index, m = 10000, level = 0.025, ncores = NULL, iseed=NULL, pkg = "ICSOutlier", qtype = 7,...)
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
    
    if(!is.null(ncores) && ncores > 1){
    
    DOTS <- list(...)
    if (is.null(DOTS$na.rm)) {Qna.rm<-FALSE} else {Qna.rm<-DOTS$na.rm}
    if (is.null(DOTS$names)) {Qnames<-7} else {Qnames<-DOTS$names}

       if(is.null(iseed)){
                if (exists(".Random.seed", envir=globalenv())){
                       oldseed <- get(".Random.seed", envir=globalenv())
                       rm(.Random.seed, envir=globalenv())
                       on.exit(assign(".Random.seed", oldseed, envir=globalenv()))
                        }
                    }
    
    type <-  "PSOCK"
    cl <- makeCluster(ncores, type=type)
    clusterExport(cl, c("n","m","MEAN","S1","S2","object","index","LEVEL","qtype","Qna.rm","Qnames","pkg","iseed"),envir=environment())
    clusterEvalQ(cl, lapply(pkg, require, character.only = TRUE))
    clusterSetRNGStream(cl = cl, iseed = iseed)
         DISTS <-parSapply(cl, 1:m, function(i,...) {
                  quantile(dist.simu.test.internal(rmvnorm(n,
                  mean = MEAN), S1 = S1, S2 = S2, S1args = object@S1args,
                  S2args = object@S2args, index = index), probs = LEVEL, type=qtype, na.rm=Qna.rm, names=Qnames) } )
    stopCluster(cl)
    
    } else {
    
     DISTS <- replicate(m, quantile(dist.simu.test.internal(rmvnorm(n, 
        mean = MEAN), S1 = S1, S2 = S2, S1args = object@S1args, 
        S2args = object@S2args, index = index), probs = LEVEL, 
        ...))
    }

    if (length(LEVEL) > 1) {
        RES <- rowMeans(DISTS)
    } else {
        RES = mean(DISTS)
        names(RES) = paste0(LEVEL * 100, "%")
    }
    return(RES)
}
