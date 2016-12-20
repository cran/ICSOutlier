dist.simu.test <- function(object, index, m = 10000, level = 0.025, ...)
        {
        LEVEL <- 1-level
        nameObject <- deparse(substitute(object))
        if (class(object) != "ics2") stop("'object' must be of class ics2")

        S1 <- get(object@S1name)
        S2 <- get(object@S2name)
        if (!is.function(S1)) stop(paste("S1 in '", nameObject, ", must be a specified as a function"))
        if (!is.function(S2)) stop(paste("S2 in '", nameObject, ", must be a specified as a function"))

        n <- nrow(object@Scores)
        p <- ncol(object@Scores)

        MEAN <- rep(0, p)
        DISTS <- replicate(m, quantile(dist.simu.test.internal(rmvnorm(n, mean = MEAN), S1 = S1, S2 = S2, 
                S1args = object@S1args, S2args = object@S2args, index = index), probs = LEVEL, ...))
        if (length(LEVEL)>1){
          RES <- rowMeans(DISTS)
        }else{
          RES = mean(DISTS)
          names(RES) = paste0(LEVEL*100, "%")
        }
        return(RES)
        }
