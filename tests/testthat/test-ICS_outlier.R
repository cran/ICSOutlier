test_that("ics.outlier and ICS_outlier - D'Agostino test", {
  X <- iris[,-5]
  
  # ics.outlier
  out_ics.outlier <- ics.outlier(ics2(X, S1 = MeanCov, S2 = Mean3Cov4), 
                                 test = "agostino.test",
                                 level.test = 0.1, level.dist = 0.05,
                                 mDist = 100, iseed = 05102023,
                                 ncores = 2)
  
  # ICS_outlier
  out_ICS_outlier <- ICS_outlier(X, S1 = ICS_cov, S2 = ICS_cov4, 
                                 S2_args = list(location = "mean3"),
                                 test = "agostino.test",
                                 level_test = 0.1, level_dist = 0.05,
                                 n_dist = 100, iseed = 05102023,
                                 n_cores = 2)
  
  # Choice of components 
  expect_equal(out_ics.outlier@index, out_ICS_outlier$index)
  expect_equal(out_ics.outlier@criterion, out_ICS_outlier$criterion)
  
  
  # Outliers
  expect_equal(out_ics.outlier@outliers,
               out_ICS_outlier$outliers)
  
  
  # Distances
  expect_equal(as.vector(out_ics.outlier@ics.distances),
               as.vector(out_ICS_outlier$ics_distances))
})

test_that("ics.outlier and ICS_outlier - Jarque Bera test", {
  X <- iris[,-5]
  
  # ics.outlier
  out_ics.outlier <- ics.outlier(ics2(X, S1 = MeanCov, S2 = Mean3Cov4), 
                                 test = "jarque.test",
                                 level.test = 0.1, level.dist = 0.05,
                                 mDist = 100, iseed = 05102023,
                                 ncores = 2)
  
  # ICS_outlier
  out_ICS_outlier <- ICS_outlier(X, S1 = ICS_cov, S2 = ICS_cov4, 
                                 S2_args = list(location = "mean3"),
                                 test = "jarque.test",
                                 level_test = 0.1, level_dist = 0.05,
                                 n_dist = 100, iseed = 05102023,
                                 n_cores = 2)
  
  # Choice of components 
  expect_equal(out_ics.outlier@index, out_ICS_outlier$index)
  expect_equal(out_ics.outlier@criterion, out_ICS_outlier$criterion)
  
  
  # Outliers
  expect_equal(out_ics.outlier@outliers,
               out_ICS_outlier$outliers)
  
  
  # Distances
  expect_equal(as.vector(out_ics.outlier@ics.distances),
               as.vector(out_ICS_outlier$ics_distances))
})


test_that("ics.outlier and ICS_outlier - Anscombe test", {
  X <- iris[,-5]
  
  # ics.outlier
  out_ics.outlier <- ics.outlier(ics2(X, S1 = MeanCov, S2 = Mean3Cov4), 
                                 test = "anscombe.test",
                                 level.test = 0.1, level.dist = 0.05,
                                 mDist = 100, iseed = 05102023,
                                 ncores = 2)
  
  # ICS_outlier
  out_ICS_outlier <- ICS_outlier(X, S1 = ICS_cov, S2 = ICS_cov4, 
                                 S2_args = list(location = "mean3"),
                                 test = "anscombe.test",
                                 level_test = 0.1, level_dist = 0.05,
                                 n_dist = 100, iseed = 05102023,
                                 n_cores = 2)
  
  # Choice of components 
  expect_equal(out_ics.outlier@index, out_ICS_outlier$index)
  expect_equal(out_ics.outlier@criterion, out_ICS_outlier$criterion)
  
  
  # Outliers
  expect_equal(out_ics.outlier@outliers,
               out_ICS_outlier$outliers)
  
  
  # Distances
  expect_equal(as.vector(out_ics.outlier@ics.distances),
               as.vector(out_ICS_outlier$ics_distances))
})

test_that("ics.outlier and ICS_outlier - Bonett test", {
  X <- iris[,-5]
  
  # ics.outlier
  out_ics.outlier <- ics.outlier(ics2(X, S1 = MeanCov, S2 = Mean3Cov4), 
                                 test = "bonett.test",
                                 level.test = 0.1, level.dist = 0.05,
                                 mDist = 100, iseed = 05102023,
                                 ncores = 2)
  
  # ICS_outlier
  out_ICS_outlier <- ICS_outlier(X, S1 = ICS_cov, S2 = ICS_cov4, 
                                 S2_args = list(location = "mean3"),
                                 test = "bonett.test",
                                 level_test = 0.1, level_dist = 0.05,
                                 n_dist = 100, iseed = 05102023,
                                 n_cores = 2)
  
  # Choice of components 
  expect_equal(out_ics.outlier@index, out_ICS_outlier$index)
  expect_equal(out_ics.outlier@criterion, out_ICS_outlier$criterion)
  
  
  # Outliers
  expect_equal(out_ics.outlier@outliers,
               out_ICS_outlier$outliers)
  
  
  # Distances
  expect_equal(as.vector(out_ics.outlier@ics.distances),
               as.vector(out_ICS_outlier$ics_distances))
})



test_that("ics.outlier and ICS_outlier - Shapiro test", {
  X <- iris[,-5]
  
  # ics.outlier
  out_ics.outlier <- ics.outlier(ics2(X, S1 = MeanCov, S2 = Mean3Cov4), 
                                 test = "shapiro.test",
                                 level.test = 0.1, level.dist = 0.05,
                                 mDist = 100, iseed = 05102023,
                                 ncores = 2)
  
  # ICS_outlier
  out_ICS_outlier <- ICS_outlier(X, S1 = ICS_cov, S2 = ICS_cov4, 
                                 S2_args = list(location = "mean3"),
                                 test = "shapiro.test",
                                 level_test = 0.1, level_dist = 0.05,
                                 n_dist = 100, iseed = 05102023,
                                 n_cores = 2)
  
  # Choice of components 
  expect_equal(out_ics.outlier@index, out_ICS_outlier$index)
  expect_equal(out_ics.outlier@criterion, out_ICS_outlier$criterion)
  
  
  # Outliers
  expect_equal(out_ics.outlier@outliers,
               out_ICS_outlier$outliers)
  
  
  # Distances
  expect_equal(as.vector(out_ics.outlier@ics.distances),
               as.vector(out_ICS_outlier$ics_distances))
})


test_that("ics.outlier and ICS_outlier - simulation test", {
  X <- iris[,-5]
  
  # ics.outlier
  out_ics.outlier <- ics.outlier(ics2(X, S1 = MeanCov, S2 = Mean3Cov4), 
                                 method = "simulation",
                                 mEig = 100,  level.test = 0.1,
                                 level.dist = 0.05,
                                 mDist = 100, iseed = 05102023,
                                 ncores = 2)
  
  # ICS_outlier
  out_ICS_outlier <- ICS_outlier(X, S1 = ICS_cov, S2 = ICS_cov4, 
                                 S2_args = list(location = "mean3"),
                                 method = "simulation",
                                 n_eig = 100,
                                 level_test = 0.1, level_dist = 0.05,
                                 n_dist = 100, iseed = 05102023,
                                 n_cores = 2)
  
  # Choice of components 
  expect_equal(out_ics.outlier@index, out_ICS_outlier$index)
  expect_equal(out_ics.outlier@criterion, out_ICS_outlier$criterion)
  
  
  # Outliers
  expect_equal(out_ics.outlier@outliers,
               out_ICS_outlier$outliers)
  
  
  # Distances
  expect_equal(as.vector(out_ics.outlier@ics.distances),
               as.vector(out_ICS_outlier$ics_distances))
})



test_that("ICS_outlier - error if not functions", {
  X <- iris[,-5]
  expect_error(ICS_outlier(X, S1 = ICS_cov, S2 = cov4(X)), 
               "S2 must be specified as a function")
  
  expect_error(ICS_outlier(X, S1 = cov(X), S2 = ICS_cov4), 
               "S1 must be specified as a function")
  
})


test_that("ICS_outlier - different scatters", {
  X <- iris[,-5]
  out_ICS_outlier <- ICS_outlier(X, S1 = ICS_cov, S2 = ICS_cov4, 
                                 test = "agostino.test",
                                 level_test = 0.1, level_dist = 0.05,
                                 n_dist = 100, iseed = 05102023,
                                 n_cores = 2)
  
  out_ICS_outlier2 <- ICS_outlier(X, S1 = MeanCov, S2 = ICS_cov4,
                           test = "agostino.test",
                           level_test = 0.1, level_dist = 0.05,
                           n_dist = 100, iseed = 05102023,
                           n_cores = 2)
  
  # Choice of components 
  expect_equal(out_ICS_outlier$index, out_ICS_outlier2$index)
  expect_equal(out_ICS_outlier$criterion, out_ICS_outlier2$criterion)
  
  
  # Outliers
  expect_equal(out_ICS_outlier$outliers,
               out_ICS_outlier2$outliers)
  
  
  # Distances
  expect_equal(as.vector(out_ICS_outlier2$ics_distances),
               as.vector(out_ICS_outlier$ics_distances))
})
