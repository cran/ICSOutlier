test_that("dist_simu_test - error", {
  X <- iris[,-5]
  icsX <- ICS(X)

  expect_error(dist_simu_test(icsX, S1 = ICS_cov, S2= ICS_cov4,
                          index = 1, m = 500),
                 "'ICS' object must has been computed with `center` = TRUE.")

})
