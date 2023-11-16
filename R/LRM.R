#'LRM
#'
#'Fit a linear regression model to inference and predict
#'
#'@param y the dependent variable
#'
#'@param x the independent invariable
#'
#'@param intercept model include the intercept term or not
#'
#'@return fitted model
#'
#'@examples
#'lrm(mtcars$mpg, mtcars$wt)
#'
#'@export
#'
#'
#'
library(Rcpp)
sourceCpp("./Cplusplus/lrm_Cpp.cpp")
lrm <- function(y,x,intercept=TRUE){
  # estimation
  # Rcpp to get X inverse
  X <<- as.matrix(x)
  Y <<- as.matrix(y)
  if(intercept == TRUE){
    X <<- cbind(rep(1,n),X)
  }
  n <<- nrow(X)
  p <<- ncol(X)
  beta.hat <<- lrm_coef(X, Y)
  Y.hat <<- X%*%beta.hat
  residual <<- Y-Y.hat
  sigmasquare.hat <<- crossprod(residual)/(n-p)
  var.beta.hat <<- diag( solve(crossprod(residual)) )*c(sigmasquare.hat)
  se.beta.hat <<- sqrt(var.beta.hat)
  t.stat <<- c(beta.hat/se.beta.hat)
  pvalue.t <<- c(2*( 1-pt(q=abs(t.stat),df=n-p) ))

  Y.bar <<- mean(Y)
  SSY <<- sum((Y-Y.bar)^2)
  SSR <<- sum((Y.hat-Y.bar)^2)
  SSE <<- sum((Y-Y.hat)^2)
  MSR <<- SSR/(p-1)
  MSE <<- SSE/(n-p)
}



#'LRM_estimate
#'
#'Estimate the parameters of a linear regression model
#'
#'@param lr.model fitted linear regression using lrm
#'
#'@return estimated parameters and standard error
#'
#'@examples
#'m = lrm(mtcars$mpg, mtcars$wt)
#'lrm.estimate(m)
#'
#'@export
#'
#'
#'
lrm.estimate <- function(lr.model){
  cbind(Estimate = c(beta.hat),
        Std.Err = se.beta.hat)
}



#'LRM_Ftest
#'
#'Test the sigificance of parameters of a fitted linear regression model
#'
#'@param lr.model the fitted linear regression model
#'
#'@param alpha the significance level
#'
#'@param interpretation print the result of hypothesis testing
#'
#'@param anova print anova table
#'
#'@return hypothesis testing result
#'
#'@examples
#'m = lrm(mtcars$mpg, mtcars$wt)
#'lrm.Ftest(m)
#'
#'@export
#'
#'
#'
lrm.Ftest <- function(lr.model, alpha=0.05, interpretation=TRUE, anova=TRUE){
  # hypothesis testing
  F.stat <<- MSR/MSE
  pvalue.F <<- 1-pf(q=F.stat, df1=(p-1), df2=(n-p))
  R.square <<- SSR/SSY

  if(anova==TRUE){
    #Source_vec <- c(model,error,total)
    SS_vec <- c(SSR, SSE, SSY)
    df_vec <- c(p-1, n-p, n-1)
    MS_vec <- c(MSR, MSE, NA)
    F.stat_vec <- c(F.stat, NA, NA)
    p.value_vec <- c(pvalue.F, NA, NA)

    anova.mat <- cbind(SS = SS_vec, df = df_vec, MS = MS_vec, F_statistic = F.stat_vec, p_value = p.value_vec)
    rownames(anova.mat) <- c("model", "error", "total")
    print(anova.mat)
  }

  if(interpretation==TRUE){
    if(pvalue.F < alpha){
      print("Under 95% significance level, We have significance evidence to reject the null hypothesis, which means that coefficients are significantly different from 0.")
    }
  }
}



#'LRM_partialtest
#'
#'Test the significance of one parameter of a fitted linear regression model
#'
#'@param lr.model the fitted linear regression model
#'
#'@param test.variable the chosen parameter to be tested
#'
#'@param alpha the significance level
#'
#'@param interpretation print the result of hypothesis testing
#'
#'@param anova print anova table
#'
#'@return hypothesis testing result
#'
#'@examples
#'m = lrm(mtcars$mpg, mtcars$wt)
#'lrm.partialtest(m, 2)
#'
#'@export
#'
#'
#'
lrm.partialtest <- function(lr.model, test.variable, alpha=0.05, confidence.interval=TRUE, interpretation=TRUE){
  # test.variable is the index of beta
  i <- test.variable
  t.mat <- cbind(Estimate = c(beta.hat[i]),Std.Err = c(se.beta.hat[i]),t_statistic = c(t.stat[i]),p_value = c(pvalue.t[i]))
  print(t.mat)

  if(interpretation==TRUE){
    if(pvalue.t[i] < alpha){
      print("Under 95% significance level, We have significance evidence to reject the null hypothesis, which means that beta is significantly different from 0.")
    }
  }

  if(confidence.interval==TRUE){
    lower.bound <- beta.hat[i]-t.stat[i]*se.beta.hat[i]
    upper.bound <- beta.hat[i]+t.stat[i]*se.beta.hat[i]

    cbind("2.5%" = c(lower.bound),
          "97.5%" = c(upper.bound))
  }
}



#'LRM_GLH
#'
#'Global hypothesis test of a subset of parameters in a fitted linear regression model
#'
#'@param lr.model the fitted linear regression model
#'
#'@param test.matrix contrast matrix
#'
#'@param interpretation print the result of hypothesis testing
#'
#'@return hypothesis testing result
#'
#'@examples
#'m = lrm(mtcars$mpg, mtcars$wt)
#'T = as.matrix(c(1,-1))
#'lrm.GLH(m)
#'
#'@export
#'
#'
#'
lrm.GLH <- function(lr.model, test.matrix, interpretation=TRUE){
  c <- rep(0, nrow(test.matrix))
  rank.T <- qr(test.matrix)$rank
  vec.GLH <- test.matrix%*%beta.hat-c
  mat.GLH <- test.matrix%*%solve(crossprod(X))%*%t(test.matrix)
  F.stat.GLH <- t(vec.GLH)%*%solve(mat.GLH)%*%vec.GLH/rank.T/MSE
  pvalue.F.GLH <- 1-pf(q=F.stat.GLH, df1=rank.T, df2=(n-p))
  GLH.mat <- cbind(F_statistic = c(F.stat.GLH), pvalue=c(pvalue.F.GLH))
  print(GLH.mat)
}
