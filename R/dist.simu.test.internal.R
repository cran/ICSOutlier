dist.simu.test.internal <-
  function(X, S1, S2, S1args, S2args, index)
  {
    ICS <- ics2(
      X,
      S1 = S1,
      S2 = S2,
      S1args = S1args,
      S2args = S2args
    )
    ics.distances(ICS, index = index)
  }


#'
#' @param X a numeric matrix or data frame containing the data to be
#' transformed.
#' @param S1 an object of class \code{"ICS_scatter"} or a function that 
#' contains the location vector and scatter matrix as \code{location} and \code{scatter} components.
#' @param S2 an object of class \code{"ICS_scatter"} or a function that 
#' contains the location vector and scatter matrix as \code{location} and \code{scatter} components.
#' @param S1_args a list containing additional arguments for \code{S1}. 
#' @param S2_args a list containing additional arguments for \code{S2}. 
#' @param index integer vector specifying which components are used to compute the  \code{\link{ics.distances}}.
#'
#' 
#' @seealso [ics_distances()], [ICS::ICS()]
#' 
#' @import ICS
#' 
#' @noRd
#'
dist_simu_test_internal <-
  function(X, S1, S2, S1_args, S2_args, index)
  {
    ICS_res <-
      ICS(
        X,
        S1 = S1,
        S2 = S2,
        S1_args = S1_args,
        S2_args = S2_args,
        center = TRUE,
        fix_signs = "scores"
      )
    ics_distances(ICS_res, index = index)
  }
