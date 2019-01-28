context("fast2PL")

test_that("fast2PL", {
  #set.seed(1)
  #M <- 1000
  #N <- 10
  #y <- matrix(rbinom(M * N, 1, 0.5), M, N)
  #impact <- as.factor(rep(c("0", "1"), each = M / 2))
  #mod <- fast2PL(y)
  #gmod <- fast2PL(y, impact = impact)
  #expect_equal(as.vector(mod$ll), -6917.756, tolerance = 1e-4)
  #expect_equal(as.vector(mod$ipars),
  #  c(0.21269, -0.05631, 0.17384, 0.13845, -0.14945, 0.26901, 0.43349, 0.15692, 0.53225, -0.46662, -0.08093, -0.07610, -0.07660, -0.03215, 0.00804, 0.08963, 0.10881, -0.08859, -0.00847, -0.02949),
  #  tolerance = 1e-4)
  #expect_equal(as.vector(gmod$ll), -6916.637, tolerance = 1e-4)
  #expect_equal(as.vector(gmod$ipars),
  #  c(0.29489, 0.11037, 0.05508, 0.10423, 0.02586, 0.22470, 0.28907, -0.00667, 0.14693, -0.57797, -0.13980, -0.09798, -0.08691, -0.05253, 0.00294, 0.04596, 0.05122, -0.08675, -0.03677, 0.07845),
  #  tolerance = 1e-4)
  #expect_equal(as.vector(gmod$mu),
  #  c(0.00000, 0.39220),
  #  tolerance = 1e-4)
  #expect_equal(as.vector(gmod$sg),
  #  c(1.00000, 1.40878),
  #  tolerance = 1e-4)
})

