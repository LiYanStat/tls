################################################################################
##
##   R package tls by Yan Li, Kun Chen and Jun Yan
##   Copyright (C) 2018-2018
##
##   This file is part of the R package tls.
##
##   The R package tls is free software: You can redistribute it and/or
##   modify it under the terms of the GNU General Public License as published
##   by the Free Software Foundation, either version 3 of the License, or
##   any later version (at your option). See the GNU General Public License
##   at <http://www.gnu.org/licenses/> for details.
##
##   The R package tls is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##
################################################################################

##' Fitting error-in-variables models via total least squares.
##'
##' It can be used to carry out regression models that account for
##' measurement errors in the independent variables.
##'
##' This function should be used with care. Confidence interval estimation
##' is given by normal approximation or bootstrap. The normal approximation and
##' bootstrap are proper when all the error terms are independent from normal
##' distribution with zero mean and equal variance (see the references for
##' more details).
##'
##' @param formula an object of class "formula" (or one that can be
##'     coerced to that class): a symbolic description of the model to
##'     be fitted.
##' @param data an optional data frame, list or environment (or object
##'     coercible by as.data.frame to a data frame) containing the
##'     variables in the model.
##' @param method method for computing confidence interval
##' @param conf.level the confidence level for confidence interval.
##' @param ... Optional arguments for future usage.
##' @return \code{tls} returns parameters of the fitted model including
##'     estimations of coefficient, corresponding estimated standard errors
##'     and confidence intervals.
##' @author Yan Li
##' @references \itemize{
##' \item  Gleser, Estimation in a Multivariate "Errors in Variables"
##' Regression Model: Large Sample Results, 1981, Ann. Stat.
##' \item Golub and Laon, An Analysis of the Total Least Squares Problem,
##' 1980, SIAM J. Numer. Anal.
##' \item Pesta, Total least squares and bootstrapping with
##' applications in calibration, 2012, Statistics.}
##' @examples
##' library(tls)
##' set.seed(100)
##' X.1 <- sqrt(1:100)
##' X.tilde.1 <- rnorm(100) + X.1
##' X.2 <- sample(X.1, size = length(X.1), replace = FALSE)
##' X.tilde.2 <- rnorm(100) + X.2
##' Y <- rnorm(100) + X.1 + X.2
##' data <- data.frame(Y = Y, X.1 = X.tilde.1, X.2 = X.tilde.2)
##' tls(Y ~ X.1 + X.2 - 1, data = data)
##' @import stats
##' @importFrom utils tail
##' @export tls
tls <- function(formula, data, method = c("normal", "bootstrap"),
                conf.level = 0.95, ...) {
  ## check whether formula is "formula"
  if (!("formula" %in% class(formula))) {
    formula <- as.formula(formula)
  }
  ## check the method for confidence interval
  method <- match.arg(method)
  if (! method %in% c("normal", "bootstrap")) {
    stop("Unknown method for computing confidence interval.")
  }
  ## check the confidence level
  if (conf.level >= 1 | conf.level <= 0) {
    stop("Confidence level should between 0 and 1.")
  }
  ## extract model
  model <- terms(formula)
  ## check intercept
  intercept <- attr(model, "intercept")
  if (intercept) {
    stop("Using total least squares requires model without intercept.")
  }
  ## check intercept
  if (all(attr(model, "order") != 1)) {
    stop("Using total least squares requires model without interactions.")
  }
  ## extract covariates
  covars <- as.character(attr(model, "term.labels"))
  ## extract repsonse
  if (attr(model, "response") == 1) {
    response <- as.character(attr(model, "variables"))[[2]]
  } else {
    stop("Cannot extract response from formula, please provide response.")
  }
  if(all(c(response, covars) %in% colnames(data))) {
    X <- as.matrix(data[, covars])
    Y <- as.matrix(data[, response])
  } else {
    stop("Cannot extract variables from data, please check data colnames.")
  }

  output <- tlsLm(X, Y, method, conf.level)
  beta.hat <- output$beta.hat
  ci.estim <- output$ci
  sd.estim <- output$sd

  beta.hat <- as.vector(beta.hat)
  names(beta.hat) <- colnames(X)
  ## confidence interval
  rownames(ci.estim) <- colnames(X)

  colnames(ci.estim) <- c(paste0(50 - conf.level * 50 , "% lower bound"),
                          paste0(conf.level * 50 + 50, "% upper bound"))
  ## result
  result <- list(coefficient = beta.hat,
                 confidence.interval = ci.estim,
                 sd.est = sd.estim)
  result
}

## estimate via Total least square approach
tlsLm <- function(X, Y, method, conf.level) {
  ## input:
  ##   X: n*p matrix, including p predictors
  ##   Y: n*1 matrix, the observations
  ## output:
  ##   a list containing the point estimate and confidence interval of
  ##   coefficient for each predictor.
  if (! is.matrix(Y)) {
    stop("Y should be a n*1 matrix")
  }
  if (dim(X)[1] != dim(Y)[1])  {  ## check size of X and Y
    stop("sizes of inputs X, Y are not consistent")
  }
  n <- dim(Y)[1]  ## number of observations
  m <- dim(X)[2]  ## number of predictors
  ## internal function for estimate
  Estls <- function(X, Y) {
    M <- cbind(X, Y)
    lambda <- tail(eigen(t(M) %*% M)$values, 1)
    sigma2.hat <- lambda / n
    Delta.hat <- (t(X) %*% X - lambda * diag(m)) / n
    ## beta.hat for the tls regression for adjusted X and Y with equal variance
    beta.hat <- as.vector(solve(t(X) %*% X - lambda * diag(m)) %*% t(X) %*% Y)
    ## [I|beta.hat]
    I.b <- cbind(diag(m), beta.hat)
    ## var.hat for beta.hat
    Var.hat <- sigma2.hat * (1 + sum(beta.hat^2)) *
      (solve(Delta.hat) + sigma2.hat * solve(Delta.hat) %*%
         solve(I.b %*% t(I.b)) %*% solve(Delta.hat)) / n
    list(beta.hat = beta.hat, Var.hat = Var.hat)
  }
  ## normal approximation
  tmp.res <- Estls(X, Y)
  beta.hat <- tmp.res$beta.hat
  var.hat <- tmp.res$Var.hat

  if (method == "normal") {
    ## standard deviation from normal approximation
    sd <- sqrt(diag(var.hat))
    ## compute z value
    z <- qnorm(0.5 + conf.level / 2)
    ## confidence interval from asymptotic normal distribution
    ci <- cbind(beta.hat - z * sd,
                     beta.hat + z * sd)
  } else {
    ## nonparamatric bootstrap
    B <- 5000
    for(i in 1:1000) {
      beta.s <- tryCatch({
        resample <- sapply(1:B,
                           function(x) {
                             sample(1:n, size = n, replace = TRUE)
                           })
        beta.s <- apply(resample, 2,
                        function(x) {
                          Xs <- X[x, ]
                          Ys <- Y[x, ]
                          Estls(Xs, Ys)$beta.hat
                        })
        if(ncol(X) == 1) {
          matrix(beta.s, ncol = B)
        } else {
          beta.s
        }
      }, error = function(e) {
        as.matrix(c(0, 0))
      })
      if (ncol(beta.s) == B) {
        break
      }
    }

    alpha <- 1 - conf.level
    sd <- t(apply(beta.s, 1, sd))
    ci <- t(apply(beta.s, 1,
                  function(x) {
                    quantile(x, c(alpha / 2, alpha / 2 + conf.level))
                  }))
  }
  list(beta.hat = beta.hat, ci = ci, sd = sd)
}
