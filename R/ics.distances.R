ics.distances <-
function(object, index = NULL)
     {
    if (!inherits(object, "ics2")) stop("'object' must be of class 'ics2'" )
    if(is.null(index)) index <- 1:ncol(object@Scores)
    DIST <- rowSums(object@Scores[, index, drop = FALSE]^2)
    return(DIST)
    }
