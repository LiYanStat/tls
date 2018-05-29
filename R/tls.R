##' Total least squares for error-in-variables models.
##'
##' This function estimates the coefficients and confidence intervals for
##' error-in-variables via total least squares.
##'
##' @param X predictors, n*p matrix.
##' @param Y response, n*1 vector.
##' @param conf.level confidence level for confidence interval.
##' @return a list of the fitted model including point estimate and
##' interval estimate of coefficients and corresponding estimate of
##' standard error.
##' @author Yan Li
##' @references \itemize{
##' \item  Gleser, Estimation in a Multivariate "Errors in Variables"
##' Regression Model: Large Sample Results, 1981, Ann. Stat.
##' \item Golub and Laon, An Analysis of the Total Least Squares Problem,
##' 1980, SIAM J. Numer. Anal.
##' \item Pesta, Total least squares and bootstrapping with
##' applications in calibration, 2012, Statistics.}
##' @details This function should be used with care. Confidence interval estimation
##' is given by normal approximation or bootstrap. The normal approximation and
##' bootstrap are proper when all the error terms are independent from normal
##' distribution with zero mean and equal variance (see the references for
##' more details).
##' @examples
##' library(tls)
##' library(MASS)
##' set.seed(100)
##' X.1 <- sqrt(1:100)
##' X.tilde.1 <- mvrnorm(mu = X.1, Sigma = diag(length(X.1)) * 2)
##' X.2 <- sample(X.1, size = length(X.1), replace = FALSE)
##' X.tilde.2 <- mvrnorm(mu = X.2, Sigma = diag(length(X.2)) * 2)
##' Y <- mvrnorm(mu = X.1 + X.2, Sigma = diag(length(X.1)) * 2)
##' X.tilde <- cbind(X.tilde.1, X.tilde.2)
##' tls(X.tilde, Y)
##' @import MASS stats utils methods
##' @export tls
tls <- function(X, Y, conf.level = 0.95) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  output <- tlsLm(X, Y, conf.level)
  beta.hat <- output$beta.hat
  ci.estim <- output$ci
  sd.estim <- output$sd

  beta.hat <- as.vector(beta.hat)
  names(beta.hat) <- colnames(X)
  ## confidence interval
  rownames(ci.estim) <- colnames(X)
  colnames(ci.estim) <- c("Boot lower bound", "Boot upper bound",
                          "Norm lower bound", "Norm upper bound")
  ## result
  result <- list(coefficient = beta.hat,
                 confidence.interval = ci.estim,
                 sd.est = sd.estim)
  result
}

## estimate via Total least square approach
tlsLm <- function(X, Y, conf.level) {
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
  ## standard deviation from normal approximation
  sd.norm <- sqrt(diag(var.hat))
  ## compute z value
  z <- qnorm(0.5 + conf.level / 2)
  ## confidence interval from asymptotic normal distribution
  ci.norm <- cbind(beta.hat - z * sd.norm,
                   beta.hat + z * sd.norm)
  ## nonparamatric bootstrap
  B <- 1000
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
        matrix(beta.s, ncol = 1000)
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

  ci.ordboot <- t(apply(beta.s, 1,
                        function(x) {
                          quantile(x, c(alpha / 2, alpha / 2 + conf.level))
                        }))

  list(beta.hat = beta.hat, ci = cbind(ci.ordboot, ci.norm), sd = sd.norm)
}
