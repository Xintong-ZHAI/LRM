% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LRM.R
\name{lrm.GLH}
\alias{lrm.GLH}
\title{LRM_GLH}
\usage{
lrm.GLH(lr.model, test.matrix, interpretation = TRUE)
}
\arguments{
\item{lr.model}{the fitted linear regression model}

\item{test.matrix}{contrast matrix}

\item{interpretation}{print the result of hypothesis testing}
}
\value{
hypothesis testing result
}
\description{
Global hypothesis test of a subset of parameters in a fitted linear regression model
}
\examples{
library(Rcpp)
sourceCpp("~/lrm/Cplusplus/lrm_Cpp.cpp")
m = lrm(mtcars$mpg, mtcars$wt)
T = matrix(c(0,1),nrow=1, ncol=2)
lrm.GLH(m, T)

}
