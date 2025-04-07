#!/usr/local/bin/Rscript
library(CVXR)
library(MASS)
library(Matrix)
library(mvtnorm)
library(expm)
library(kernlab)
library(glmnet)
library(hal9001)

qcqp <- function(b, V, pen, roughness, tol = 1e-03) {
  
  lam.old <- 1/roughness
  lam.converge <- FALSE
  step <- 1e-06
  max.iter <- 50
  iter <- 0

  while(!lam.converge) {
    iter <- iter + 1
    a <- ginv(V + lam.old * pen) %*% b
    a.1 <- ginv(V + (lam.old + step) * pen) %*% b
    a.2 <- ginv(V + (lam.old - step) * pen) %*% b

    deriv.approx <-  (((t(a.1) %*% pen %*% a.1)/(t(a.1) %*% V %*% a.1)) -
                      ((t(a.2) %*% pen %*% a.2)/(t(a.2) %*% V %*% a.2)))/
                     (2 * step)
    lam <- lam.old - c((t(a) %*% pen %*% a)/(t(a) %*% V %*% a) - roughness)/c(deriv.approx)
    if(lam <= 0) {lam <- lam.old * .5}

    lam.converge <- abs(log(lam) - log(lam.old)) < tol
    lam.old <- lam
    if(iter > max.iter) break
  }
  
  a <- (ginv(V + lam.old * pen) %*% b)
  test.stat <-  c(abs(sum(b * a))/sqrt(t(a) %*% V %*% a))
  out <- list(test.stat = test.stat,
              lam = lam,
              coef = a/c(sqrt(t(a) %*% V %*% a)))
  # out <- test.stat
  return(out)
}

#' Conditional density estimation
#'
#' Estimate the conditional density using kernel smoothing
#' 
#' @param x Predictor of interest, an n-dimensional vector.
#' @param W Covariates for adjustment, an n x p matrix.
#' @param no.folds Number of folds to be used for selection of bandwidth via cross-validation
#' @param no.eval Number of evaluation points (levels of predictor of interest) at which the conditional density will be estimated.
#'                Estimates at intermediate points are estimated using linear interpolation.
#' @param bw Value for bandwidth, if selected a priori. Set to NULL by default.
#'
#'  @return A list containing the following: 
#'         \describe{
#'            \item{cond.dens} Estimate of conditional density of predictor of interest, given covariates.
#'            \item{marg.dens} Estimate of marginal density of predictor of interest.
#'            \item{prop.score} Ratio of marginal density to conditional density.
#'          }             
#' @return An n-dimensional vector containing the conditional density, marginal density, and their ratior.
#'         
#' @export
EstConditionalDensity <- function(x, W, no.folds = 10, no.eval = 10, bw = NULL, continuous = TRUE) {
  n <- length(x)
  x.eval <- seq(min(x), max(x), length.out = no.eval)
  cond.dens <- matrix(NA, n, n) # set cond.dens[i,j] = cond.dens(a[i], w[j])
  
  if(continuous) {
    # use bandwidth that is optimal for kernel estimation of margianl density
    # is sensible choice if con'l density is about as smooth in a as marginal density
    # could try other choices via cross-validation, but this is comp'l cheaper
    # con'l density estimation is generally hard, and needs further study
    
    if(is.null(bw)) {
      # select bandwidth via cross-validation
      folds <- sample(cut(1:n, breaks = no.folds, labels = FALSE))
      bw.seq <- exp(seq(log(.01 * sd(x)), log(25 * sd(x)), length.out = 100))
    
      pred.error <- matrix(NA, nrow = 100, ncol = no.folds) # prediction error
      for(i in 1:no.folds) {
        for(j in 1:100) {
          x.train <- x[folds != i]
          x.test <- x[folds == i]
          
          n.test <- length(x.test)
          fit.test <- numeric(n.test)
          
          for(k in 1:n.test) {
            kern <- dnorm(x[folds != i], mean = x.test[k], sd = bw.seq[j])
            fit.test[k] <- mean(kern)
          }
          
          pred.error[j,i] <- mean(-log(fit.test))
        }
      }
      avg.pred.error <- apply(pred.error, 1, mean)
      se.pred.error <- apply(pred.error, 1, sd)/sqrt(no.folds)
      bw.opt <- max(bw.seq[avg.pred.error < min(avg.pred.error) + se.pred.error[which.min(avg.pred.error)]])
      # bw.opt <- bw.seq[which.min(apply(pred.error, 1, mean))]
     } else {
     # set bandwidth to a pre-specified value, if one is provided
      bw.opt <- bw
     }
    
    # marg.dens <- numeric(n)
    # for(i in 1:n) {
    #   marg.dens[i] <- mean(dnorm(x, mean = x[i], sd = 2 * bw.opt))
    # }
    
    f.eval <- matrix(NA, nrow = n, ncol = no.eval)
    for(j in 1:no.eval) {
      kern.xj <- dnorm(x, mean = x.eval[j], sd = bw.opt)
      # hal.fit <- fit_hal(X = w, Y = kern.aj)
      hal.fit <- fit_hal(X = W, Y = kern.xj,
                         family = "poisson", return_x_basis = TRUE)
      f.eval[,j] <- predict(hal.fit, new_data = W)
      # f.eval[,j] <- exp(hal.fit$X.basis %*% hal.fit$coefs)
    }
    
    # set cond.dens[i,j] = cond.dens(a[i], w[j])
    for(i in 1:n) {
      k <- min(which(x.eval - x[i] >= 0))
      for(j in 1:n) {
        if(k == 1) {
          cond.dens[i,j] <- f.eval[j,1]
        } else if(k == no.eval) {
          cond.dens[i,j] <- f.eval[j,no.eval]
        } else {
          t <- (x.eval[k] - x[i])/(x.eval[k] - x.eval[k-1])
          cond.dens[i,j] <- f.eval[j,k] * (1-t) + f.eval[j,k-1] * t
        }
      }
    }
  } else { # x is binary
    hal.fit <- fit_hal(X = W, Y = x, family = "binomial")
    predict.hal.fit <- predict(hal.fit, new_data = W)
    for(i in 1:n) {
      cond.dens[i,] <- predict.hal.fit * x[i] + (1 - predict.hal.fit) * (1 - x[i])
    }
  }
  
  marg.dens <- apply(cond.dens, 1, mean)
  out <- list(cond.dens = cond.dens, marg.dens = marg.dens)
  
  return(out)
}

KernelEigenBasis <- function(x, W, tau = sqrt(.2/ncol(W)), eigenval.min = 1e-03) {
  
  weights <- c(sd(x)^(-1), tau * (apply(W, 2, sd))^(-1))
  
  K <- kernelMatrix(x = cbind(x, W) %*% diag(weights), kernel = rbfdot(sigma = .5))
  
  # Obtain the eigen-decomposition of the kernel matrix
  K.eigenvals <- eigen(K, only.values = TRUE)
  basis.dim <- min(c(max(which(K.eigenvals$values >= eigenval.min)), 100))
  K.svd <- svd(K, nu = basis.dim, nv = basis.dim)
  
  
  eigenbasis <- K.svd$u # eigen-vectors
  eigenvals <- K.svd$d[1:basis.dim] # eigen-values
  
  out <- list(K = K,
              eigenbasis = eigenbasis,
              eigenvals = eigenvals,
              weights = weights)
  return(out)
}

CenterEigenBasis <- function(x, W, weights, eigenbasis, eigenvals, continuous = TRUE) {
  
  n <- nrow(W)
  
  dens.est <- EstConditionalDensity(x = x, W = W, continuous = continuous)
  cond.dens <- dens.est$cond.dens
  marg.dens <- dens.est$marg.dens
  
  eigenbasis.pred <- matrix(NA, nrow = n, ncol = ncol(eigenbasis))
  
  for(l in 1:ncol(eigenbasis)) {
    # evaluate l-th eigen function at each (X[j], W[,i])
    eigenfunction.eval <- matrix(NA, nrow = n, ncol = n)
    for(i in 1:n) {
      for(j in 1:n) {
        diff.mat <- matrix(c(x[j], W[i,]), nrow = nrow(W), ncol = ncol(W) + 1, byrow = T) -
                    cbind(x, W)
        k.vec <- exp(-diff.mat^2 %*% weights^2/2)
        eigenfunction.eval[i, j] <- sum(k.vec * eigenbasis[,l])/eigenvals[l]
      }
      # estimate conditional mean
      prop.score <- cond.dens[,i]/marg.dens
      eigenbasis.pred[i,l] <- mean(eigenfunction.eval[i,] * prop.score)
    }
  }
  
  out <- eigenbasis.pred
  return(out)
}

VimpEstNuisance <- function(y, x, W, tau = sqrt(.2/ncol(W)),
                            eigenval.min = 1e-03, continuous = TRUE) {
  
  n <- length(y)
  
  # Estimate conditional mean of Y given W
  hal.fit.y <- fit_hal(X = W, Y = y,
                       family = "gaussian", return_x_basis = TRUE)
  y.pred <- predict(hal.fit.y, new_data = W)
  
  # Construct basis expansion for f(X, W)
  eigenquantities <- KernelEigenBasis(x = x, W = W, 
                                      tau = tau, eigenval.min = eigenval.min)
  eigenbasis <- eigenquantities$eigenbasis
  eigenvals <- eigenquantities$eigenvals
  weights <- eigenquantities$weights
  
  # Estimate the conditional mean of all basis functions
  eigenbasis.pred <- CenterEigenBasis(x = x, W = W, continuous = continuous, weights = weights,
                                      eigenbasis = eigenbasis, eigenvals = eigenvals)
  
  V <- t(eigenbasis) %*% eigenbasis/n
  
  out <- list(y.pred = y.pred,
              eigenbasis = eigenbasis,
              eigenbasis.pred = eigenbasis.pred,
              eigenvals = eigenvals,
              V = V)
  
  return(out)
}


# Select roughness using cross-validation
VimpSelectRoughness <- function(y, x, W,
                                y.pred, eigenbasis, eigenbasis.pred,
                                V, eigenvals) {
  
  y.cent <- y - y.pred
  eigenbasis.cent <- eigenbasis - eigenbasis.pred
  
  pen.mat.sqrtinv <- diag(eigenvals^(1/2))
  pen.mat.sqrt <- diag(eigenvals^(-1/2))
  pen.mat <- diag(eigenvals^(-1))
  
  cv.glmnet.fit <- cv.glmnet(x = eigenbasis.cent %*% pen.mat.sqrtinv, y = y.cent, 
                             family = "gaussian", alpha = 0, standardize = FALSE)
  lam.1se <- cv.glmnet.fit$lambda.1se
  lam.cv <- cv.glmnet.fit$lambda.min
  glmnet.fit.1se <- glmnet(x = eigenbasis.cent %*% pen.mat.sqrtinv, y = y.cent,
                           family = "gaussian", lambda = lam.1se, alpha = 0, standardize = FALSE)
  glmnet.fit.cv <- glmnet(x = eigenbasis.cent %*% pen.mat.sqrtinv, y = y.cent,
                          family = "gaussian", lambda = lam.cv, alpha = 0, standardize = FALSE)
    
  coef.1se <- as.vector(pen.mat.sqrt %*% as.vector(glmnet.fit.1se$beta))
  coef.cv <- as.vector(pen.mat.sqrt %*% as.vector(glmnet.fit.cv$beta))
  
  # this is a very stupid rule and will need to change in the future. but that's ok
  roughness.1se <- c(t(coef.1se) %*% pen.mat %*% coef.1se)/c(t(coef.1se) %*% V %*% coef.1se)
  roughness.1se <- max(c(roughness.1se, 100))
  roughness.cv <- c(t(coef.cv) %*% pen.mat %*% coef.cv)/c(t(coef.cv) %*% V %*% coef.cv)
  roughness.cv <- max(c(roughness.cv, 100))

  # out <- list(roughness.1se, roughness.cv)
  out <- roughness.1se
  return(out)
}

# estimate nuisance parameters
# return nuisance, plug-in, and EIF
VimpOneStep <- function(y, x, W,
                        y.pred, eigenbasis, eigenbasis.pred,
                        V, eigenvals, roughness) {
  n <- length(y)
  
  y.cent <- y - y.pred
  eigenbasis.cent <- eigenbasis - eigenbasis.pred

  pen.mat <- diag(eigenvals^(-1))
  gram <- t(eigenbasis) %*% eigenbasis/n
  inner.prod <- t(y.cent) %*% (eigenbasis.cent)/n
  
  cvx.result <- qcqp(c(inner.prod), V = V, pen = pen.mat, roughness = roughness, tol = 1e-03)
  cvx.val <- cvx.result$test.stat
  cvx.sol <- cvx.result$coef
  coef.opt <- sum(c(inner.prod) * c(cvx.sol)) * cvx.sol
  
  null.gof <- mean(y.cent^2)
  alt.gof <- mean((y.cent - eigenbasis %*% coef.opt)^2 + 
                   2 * y.cent * eigenbasis.pred %*% coef.opt)
  iif <- null.gof - alt.gof
  
  eif <- (y.cent^2) - ((y.cent - eigenbasis %*% coef.opt)^2 + 
                        2 * y.cent * eigenbasis.pred %*% coef.opt)
  eif <- c(eif - mean(eif))
  
  out <- list(eif = eif, 
              iif = iif)
  
  return(out)
}

# Get bias-corrected and "naive" estimators for the improvement in fit for the vimp problem
VimpTest <- function(y, x, W, 
                     y.pred, eigenbasis, eigenbasis.pred,
                     V, eigenvals, roughness, no.bs = 1000, alpha = .05, 
                     alt.eif, iif) {
  
  y.cent <- y - y.pred
  eigenbasis.cent <- eigenbasis - eigenbasis.pred
  
  n <- length(y)
  null.bs.samples <- numeric(no.bs)
  alt.bs.samples <- numeric(no.bs)
  
  pen.mat <- diag(eigenvals^(-1))
  null.eif <- -2 * diag(c(y.cent)) %*% eigenbasis.cent +
               2 * matrix(apply(diag(c(y.cent)) %*% eigenbasis, 2, mean),
                          nrow = n, ncol = ncol(eigenbasis.cent),
                          byrow = TRUE)
  
  null.eif <- (diag(n) - matrix(1/n, nrow  = n, ncol = n)) %*% null.eif
  
  # Generate multiplier bootstrap samples under zero importance and pos importance settings
  for(i in 1:no.bs) {
    mult <- 2 * rbinom(n, 1, .5) - 1
    # mult <- rnorm(n)
    mult <- mult - mean(mult)
    mult.null.eif <- diag(mult) %*% null.eif
    mult.null.eif.avg <- c(apply(mult.null.eif, 2, mean))
    
    cvx.result <- qcqp(mult.null.eif.avg, V = V, pen = pen.mat, roughness = roughness, tol = 1e-03)
    null.bs.samples[i] <- cvx.result$test.stat^2/4
    
    alt.bs.samples[i] <- mean(alt.eif * mult)
    
    # if(i %% 100 == 0) print(paste0(i, " out of ", no.bs))
    
  }
  
  p.val <- mean(iif < null.bs.samples, na.rm = TRUE)
  
  # Get BS interval
  # weight <- ifelse(iif <= log(n) * mean(null.bs.samples), 
  #                  1, mean(iif <= null.bs.samples, na.rm = TRUE))
  weight <- ifelse(iif <= log(n) * mean(null.bs.samples), 1, 0)
  bs.cent <- weight * sqrt(null.bs.samples) + (1 - weight) * abs((2 * sqrt(iif))^(-1) * alt.bs.samples)
  Sn <- quantile(bs.cent, 1 - alpha, na.rm = TRUE)
  ci.upper <- (Sn + sqrt(iif))^2
  ci.lower <- max(c(sqrt(iif) - Sn, 0))^2
  
  
  # Get BS interval and estimate with bias correction
  bias <- mean(null.bs.samples)
  iif.debias <- ifelse(iif <= log(n) * mean(null.bs.samples),
                       max(c(iif - bias, 0)),
                       iif)
  bs.cent.debias <- NULL
  if(iif <= log(n) * mean(null.bs.samples)) {
    bs.cent.debias <- sqrt(apply(cbind(null.bs.samples - bias, 0), 1, max))
  } else {
    bs.cent.debias <- abs((2 * sqrt(iif.debias))^(-1) * alt.bs.samples)
  }
                           
  # bs.cent.debias <- weight * sqrt(apply(cbind(null.bs.samples, 0), 1, max)) + 
  #                   (1 - weight) * abs((2 * sqrt(iif.debias))^(-1) * alt.bs.samples)
  Sn.debias <- quantile(bs.cent.debias, 1 - alpha, na.rm = TRUE)
  ci.upper.debias <- (sqrt(iif.debias) + Sn.debias)^2
  ci.lower.debias <- max(c(sqrt(iif.debias) - Sn.debias, 0))^2
  
  out <- list(iif = iif,
              p.val = p.val,
              ci.lower = ci.lower,
              ci.upper = ci.upper,
              iif.debias = iif.debias,
              ci.lower.debias = ci.lower.debias,
              ci.upper.debias = ci.upper.debias,
              null.bs.samples = null.bs.samples,
              alt.bs.samples = alt.bs.samples,
              bs.cent = bs.cent,
              null.eif = null.eif,
              alt.eif = alt.eif)
  return(out)
}


Vimp <- function(y, x, W, continuous = TRUE,
                 tau = sqrt(.2/ncol(W)),
                 eigenval.min = 1e-03, roughness = NULL, no.bs = 1000,
                 alpha = .05) {
  
  nuisance <- VimpEstNuisance(y = y, x = x, W = W, tau = tau, continuous = continuous,
                              eigenval.min = eigenval.min)
  
  if(is.null(roughness)) {
    roughness <- VimpSelectRoughness(y = y, x = x, W = W,
                                     y.pred = nuisance$y.pred,
                                     eigenbasis = nuisance$eigenbasis, 
                                     eigenbasis.pred = nuisance$eigenbasis.pred,
                                     V = nuisance$V, eigenvals = nuisance$eigenvals)
  }
  
  onestep <- VimpOneStep(y = y, x = x, W = W,
                         y.pred = nuisance$y.pred, 
                         eigenbasis = nuisance$eigenbasis,
                         eigenbasis.pred = nuisance$eigenbasis.pred,
                         V = nuisance$V, eigenvals = nuisance$eigenvals,
                         roughness = roughness)
  
  test <- VimpTest(y = y, x = x, W = W, 
                   y.pred = nuisance$y.pred,
                   eigenbasis = nuisance$eigenbasis,
                   eigenbasis.pred = nuisance$eigenbasis.pred,
                   V = nuisance$V, eigenvals = nuisance$eigenvals,
                   roughness = roughness, no.bs = no.bs, alpha = alpha,
                   iif = onestep$iif, alt.eif = onestep$eif)
  
  out <- list(estimate = onestep$iif,
              p.val = test$p.val,
              ci.lower = test$ci.lower,
              ci.upper = test$ci.upper,
              estimate.debias = test$iif.debias,
              ci.lower.debias = test$ci.lower.debias,
              ci.upper.debias = test$ci.upper.debias,
              null.bs.samples = test$null.bs.samples,
              alt.bs.samples = test$alt.bs.samples,
              roughness = roughness)
  
  return(out)
  
}

VimpSplit <- function(y, x, W, alpha = .05) {
  n <- length(y)
  
  split <- sample(1:2, n, replace = TRUE)
  
  n.red <- sum(split == 1)
  y.red <- y[split == 1]
  W.red <- W[split == 1,]
  
  n.full <- sum(split == 2)
  y.full <- y[split == 2]
  W.full <- W[split == 2,]
  x.full <- x[split == 2]
  
  # Reduced fit
  hal.fit.red <- fit_hal(X = W.red, Y = y.red, yolo = FALSE,
                         family = "gaussian")
  est.red <- predict(hal.fit.red, new_data= W.red)
  
  # Full fit
  hal.fit.full <- fit_hal(X = cbind(W.full, x.full), Y = y.full, yolo = FALSE,
                         family = "gaussian")
  est.full <- predict(hal.fit.full, new_data = cbind(W.full, x.full))
  
  # Perform test
  iif <- (mean((y.red - est.red)^2) - mean((y.full - est.full)^2))
  # iif <- max(c(iif, 0))
  se <-  sqrt(var((y.red - est.red)^2)/n.red + var((y.full - est.full)^2)/n.full)
  p.val <- pchisq(iif^2/se^2, df = 1, ncp = 0, lower.tail = FALSE)
  ci.lower <- iif - qnorm(1-alpha/2, mean = 0, sd = 1) * se
  ci.upper <- iif + qnorm(1-alpha/2, mean = 0, sd = 1) * se
  # ci.lower <- max(c(ci.lower, 0))
  
  # Return output
  out <- list(estimate = iif, 
              se = se,
              p.val = p.val,
              ci.upper = ci.upper,
              ci.lower = ci.lower)
}


### Examples
# n <- 100
# rho <- .5
# p <- 3
# Sigma <- matrix(rho, p, p) + (1 - rho) * diag(p)
# Predictors <- 2 * (pnorm(rmvnorm(n, mean = rep(0, p), sigma = Sigma)) - .5)
# y <- 2 * sin(pi * (Predictors[,1])) + 2 * exp(4 * (Predictors[,1] + Predictors[,2]))/
#                                          (exp(4 * (Predictors[,1] + Predictors[,2])) + 1) + rnorm(n)
# 
# out <- Vimp(y = y, x = as.matrix(Predictors[,3]), W = Predictors[,1:2],
#             roughness = 1000)
# out.adap <- Vimp(y = y, x = as.matrix(Predictors[,3]), W = Predictors[,1:2])
# 
# 
# out.split <- VimpSplit(y = y, x = as.matrix(Predictors[,3]), W = Predictors[,1:2])
