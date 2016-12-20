setClass("icsOut", representation(outliers = "integer",  
                                  ics.distances = "numeric", 
                                  ics.dist.cutoff = "numeric", 
                                  level.dist = "numeric", 
                                  level.test = "numeric",
                                  method = "character",
                                  index = "numeric",
                                  test = "character",
                                  criterion = "numeric",
                                  adjust = "logical",
                                  type = "character",
                                  mDist = "integer",
                                  mEig = "integer",
                                  S1name = "character",
                                  S2name = "character"))




setMethod("show", signature(object = "icsOut"),
function(object)
    {
    comps <- object@index
    if (sum(comps)<0.5) {
        ncomps <- 0
        print(paste(ncomps, " components were selected and no outliers were detected.", sep = ""))
        } else {
        ncomps <- length(comps)
        print(paste(ncomps, " components were selected and ", sum(object@outliers), " outliers were detected.", sep = ""))
        }
    #invisible(object@outliers)
    }
)


setMethod("summary", signature(object = "icsOut"),
function(object, digits = 4)
    {
    comps <- object@index
    if (sum(comps)<0.5) ncomps <- 0 else ncomps <- length(comps) 
    if (object@method == "norm.test") METHOD <- paste(object@method, " (", object@test, ")", sep = "") else METHOD <- object@method
    cat("\nICS based on two scatter matrices and two location estimates\n")
    cat("S1: ", object@S1name)
    cat("\nS2: ", object@S2name)
    cat("\n")
    cat("\nSearching for a small proportion of outliers\n")
    cat("\n")
    #cat(paste("Components selected: ", ncomps, sep = "") )
    #cat("\n")
    cat(paste("Components selected at nominal level ", object@level.test, ": ", ncomps, sep = "") )
    cat("\n")
    cat(paste("Selection method: ", METHOD, sep = ""))
    cat("\n")
    #cat(paste("Number of outliers: ", sum(object@outliers), sep = "")) 
    cat(paste("Number of outliers at nominal level ", object@level.dist, ": ", sum(object@outliers), sep = "")) 
    cat("\n")  
    invisible(object)
    }
)

#setMethod("plot", signature(x = "icsOut", y = "missing"),
#function(x, pch.out = 16, pch.good = 4, col.out = 1, col.good = grey(0.5), col.cut = 1, 
#         lwd.cut = 1, lty.cut = 1, xlab = "Observation Number", ylab = "ICS distances", ...)
#    {
#    colPoints <- ifelse(x@outliers == 1L, col.out, col.good)
#    pchPoints <- ifelse(x@outliers == 1L, pch.out, pch.good)
#    
#    plot(x@ics.distances, col = colPoints, pch = pchPoints, xlab = xlab, ylab = ylab, ...)
#    abline(h = x@ics.dist.cutoff, col = col.cut, lwd = lwd.cut, lty = lty.cut)
#    }
#)


setMethod("plot", signature(x = "icsOut", y = "missing"),
function(x, pch.out = 16, pch.good = 4, col.out = 1, col.good = grey(0.5), col.cut = 1, 
         lwd.cut = 1, lty.cut = 1, xlab = "Observation Number", ylab = "ICS distances", ...)
    {
    
    YESylim <- hasArg("ylim") 
    if (sum(x@outliers)>0.5){
            colPoints <- ifelse(x@outliers == 1L, col.out, col.good)
            pchPoints <- ifelse(x@outliers == 1L, pch.out, pch.good)
            plot(x@ics.distances, col = colPoints, pch = pchPoints, xlab = xlab, ylab = ylab, ...)
            abline(h = x@ics.dist.cutoff, col = col.cut, lwd = lwd.cut, lty = lty.cut)
            } else {
                if (sum(x@ics.distances) > 0){
                    if (!YESylim) {
                        plot(x@ics.distances, col = col.good, pch = pch.good, xlab = xlab, ylab = ylab, ylim = c(0, x@ics.dist.cutoff), ...)
                        }  else {
                        plot(x@ics.distances, col = col.good, pch = pch.good, xlab = xlab, ylab = ylab, ...)
                        }
                    abline(h = x@ics.dist.cutoff, col = col.cut, lwd = lwd.cut, lty = lty.cut)
                    } else {
                    plot(x@ics.distances, xlab = xlab, ylab = ylab, type = "n", ...)
                    text(length(x@ics.distances)/2, 0, labels = c("No components have been selected for outlier detection.\n There is nothing to plot."))
                    }
            }
    }
)
