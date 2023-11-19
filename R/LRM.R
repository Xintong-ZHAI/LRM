#'LRM
#'
#'@description
#'This package is used to fit a linear regression model to inference and predict
#'C++ is used to accelerate calculation
#'
#'@importFrom stats pf pt qt
#'
#'@param y the dependent variable
#'
#'@param x the independent invariable
#'
#'@param test.variable the chosen parameter to be tested, default is the first term
#'
#'@param test.matrix contrast matrix for GLH, default is an empty matrix
#'
#'@param intercept "TRUE" means model include the intercept term, "FALSE" means not
#'
#'@param estimate "TRUE" means calculate the point estimation of parameters, "FALSE" means not
#'
#'@param alpha the significance level, default is 0.05
#'
#'@param interpretation "TRUE" means print the result of hypothesis testing, "FALSE" means not
#'
#'@param anova "TRUE" means print anova table, "FALSE" means not
#'
#'@param Ftest "TRUE" means conduct F test, "FALSE" means not
#'
#'@param partialtest "TRUE" means conduct partial test, "FALSE" means not
#'
#'@param confidence.interval "TRUE" means print the confidence interval, "FALSE" means not
#'
#'@param GLH "TRUE" means conduct GLH, "FLASE" means not
#'
#'@return result of one specific hypothesis testing
#'
#'@examples
#'library(Rcpp)
#'sourceCpp("~/lrm/src/lrm_Cpp.cpp")
#'lrm(mtcars$mpg, mtcars$wt)
#'
#'@export
#'
#'
#'
lrm <- function(y, x, test.variable=1, test.matrix=matrix(c(0)), intercept=TRUE, estimate=TRUE,
                alpha=0.05, interpretation=FALSE, anova=FALSE, Ftest=FALSE, partialtest=FALSE,
                confidence.interval=FALSE, GLH=FALSE){
  # estimation
  # Rcpp to get X inverse
  X <- as.matrix(x)
  Y <- as.matrix(y)
  n <- nrow(X)
  if(intercept == TRUE){
    X <- cbind(rep(1,n),X)
  }
  p <- ncol(X)
  beta.hat <- lrm_coef(X, Y)
  Y.hat <- X%*%beta.hat
  residual <- Y-Y.hat
  sigmasquare.hat <- crossprod(residual)/(n-p)
  var.beta.hat <- diag( solve(crossprod(X)) )*c(sigmasquare.hat)
  se.beta.hat <- sqrt(var.beta.hat)
  t.stat <- c(beta.hat/se.beta.hat)
  pvalue.t <- c(2*( 1-pt(q=abs(t.stat),df=n-p) ))

  Y.bar <- mean(Y)
  SSY <- sum((Y-Y.bar)^2)
  SSR <- sum((Y.hat-Y.bar)^2)
  SSE <- sum((Y-Y.hat)^2)
  MSR <- SSR/(p-1)
  MSE <- SSE/(n-p)

  if(estimate == TRUE){
    result <- cbind(Estimate = c(beta.hat),
                    Std.Err = se.beta.hat)
  }

  if(Ftest == TRUE){
    # hypothesis testing
    F.stat <- MSR/MSE
    pvalue.F <- 1-pf(q=F.stat, df1=(p-1), df2=(n-p))
    R.square <- SSR/SSY

    if(anova==TRUE){
      #Source_vec <- c(model,error,total)
      SS_vec <- c(SSR, SSE, SSY)
      df_vec <- c(p-1, n-p, n-1)
      MS_vec <- c(MSR, MSE, NA)
      F.stat_vec <- c(F.stat, NA, NA)
      p.value_vec <- c(pvalue.F, NA, NA)

      anova.mat <- cbind(SS = SS_vec, df = df_vec, MS = MS_vec, F_statistic = F.stat_vec, p_value = p.value_vec)
      rownames(anova.mat) <- c("model", "error", "total")
      result <- anova.mat
    }

    if(interpretation==TRUE){
      if(pvalue.F < alpha){
        print("Under 95% significance level, We have significance evidence to reject the null hypothesis, which means that coefficients are significantly different from 0.")
      }
    }
  }

  if(partialtest==TRUE){
    # test.variable is the index of beta
    i <- test.variable
    t.mat <- cbind(Estimate = c(beta.hat[i]),Std.Err = c(se.beta.hat[i]),t_statistic = c(t.stat[i]),p_value = c(pvalue.t[i]))
    result <- t.mat

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

  if(GLH==TRUE){
    c <- rep(0, nrow(test.matrix))
    rank.T <- qr(test.matrix)$rank
    vec.GLH <- test.matrix%*%beta.hat-c
    mat.GLH <- test.matrix%*%solve(crossprod(X))%*%t(test.matrix)
    F.stat.GLH <- t(vec.GLH)%*%solve(mat.GLH)%*%vec.GLH/rank.T/MSE
    pvalue.F.GLH <- 1-pf(q=F.stat.GLH, df1=rank.T, df2=(n-p))
    GLH.mat <- cbind(F_statistic = c(F.stat.GLH), pvalue=c(pvalue.F.GLH))
    result <- GLH.mat
  }

  return(result)
}


