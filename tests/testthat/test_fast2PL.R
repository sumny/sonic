context("fast2PL")

test_that("fast2PL", {
  data(LSAT7, package = "sonic")
  impact <- as.factor(rep(c(1, 2), dim(LSAT7)[1] / 2))
  mod <- fast2PL(LSAT7)
  gmod <- fast2PL(LSAT7, impact = impact)
  expect_equal(as.vector(mod$ll), -2658.805, tolerance = 1e-2)
  expect_equal(as.vector(mod$ipars),
    c(0.9890370, 1.0796213, 1.7038597, 0.7662276, 0.7360880,
    1.8568123, 0.8078543, 1.8035323, 0.4862891, 1.8547752), tolerance = 1e-2)
  expect_equal(as.vector(gmod$ll), -2658.799, tolerance = 1e-2)
  expect_equal(as.vector(gmod$ipars),
    c(0.9909031, 1.0844801, 1.7124429, 0.7675232, 0.7382976,
    1.8539057, 0.8058388, 1.8013606, 0.4844957, 1.8530637), tolerance = 1e-2)
  expect_equal(as.vector(gmod$mu),
    c(0.00000, 0.003723521), tolerance = 1e-2)
  expect_equal(as.vector(gmod$sg),
    c(1.00000, 0.9927122), tolerance = 1e-2)
})

