library(Rcpp)
library(lrm)
sourceCpp("~/lrm/Cplusplus/lrm_Cpp.cpp")

test_that("lrm works", {
  m1 = lrm(mtcars$mpg, mtcars$wt)
  expect_equivalent(round(lrm.estimate(m1)[,1],3),
               round(lm(mpg~wt, data=mtcars)$coefficients,3))

  m2 = lrm(CO2$uptake, CO2$conc)
  expect_equivalent(round(lrm.estimate(m2)[,1],3),
                    round(lm(uptake~conc, data=CO2)$coefficients,3))
})
