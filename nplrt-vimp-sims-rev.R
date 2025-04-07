#!/usr/local/bin/Rscript
source("nplrt-vimp-rev.R")

job.id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
seed <- 206 * job.id
set.seed(seed)

no.reps <- 1
ns <- c(100, 200, 400, 800, 1600)
noise.sd <- 1

#########################
### Simulation Settings
#########################

cond.dens <- function(x, rho = .5) {
  Sigma <- matrix(rho, p, p) + (1 - rho) * diag(p)
  Sigma.red <- Sigma[1:2, 1:2]
  cond.dens <- dmvnorm(qnorm((x + 1)/2), sigma = Sigma) * 
               (2 * dmvnorm(qnorm((x[1:2] + 1)/2), sigma = Sigma.red) * dnorm(qnorm((x[3] + 1)/2)))^(-1)
  return(cond.dens)
}

# theta0 takes as input a 3-dim vector x
# const <- mean((seq(-1, 1, length.out = 10000) - .5)^2)
const <- 0

theta0.list <- list()
null.list <- list()

# Setting 1: Null Holds
theta0.list[[1]] <- function(x, noise.sd = 2) {
  out <- 2 * sin(pi * (x[1])) - 2 * ((x[2] - .5)^2 - const) +
  rnorm(1) * noise.sd
  return(out)
}

null.list[[1]] <- function(x) {
  out <- 2 * sin(pi * (x[1])) - 2 * ((x[2] - .5)^2 - const)
  return(out)
}

# Setting 2: Alt Holds (additive)
theta0.list[[2]] <- function(x, noise.sd = 2) {
  # out <- 2 * sin(pi * (x[1])) - 2 * ((x[2] - .5)^2 - const) + 2 * exp(-x[3]^2 * 4) +
  # runif(1, -2, 2) * noise.sd
  out <- 2 * sin(pi * (x[1])) - 2 * ((x[2] - .5)^2 - const) + 2 * exp(-x[3]^2 * 4) +
  rnorm(1) * noise.sd
  return(out)
}

null.list[[2]] <- function(x, rho = .5) {
  x3 <- seq(-.9995, .9995, length.out = 1000)
  out <- 2 * sin(pi * (x[1])) - 2 * ((x[2] - .5)^2 - const) + 
         mean(apply(cbind(x[1], x[2], x3), 1, cond.dens) * 2 * exp(-x3^2 * 4)) * 1.999
  return(out)
}

# Setting 3: Alt Holds (interaction)
theta0.list[[3]] <- function(x, noise.sd = 2) {
  # out <- 2 * sin(pi * (x[1])) - 2 * ((x[2] - .5)^2 - const) + 3 * exp(-(x[1] + x[3])/2) +
  # runif(1, -2, 2) * noise.sd
  out <- 2 * sin(pi * (x[1])) - 2 * ((x[2] - .5)^2 - const) + 3 * exp(-(x[1] + x[3])/2) +
  rnorm(1) * noise.sd
  return(out)
}

null.list[[3]] <- function(x, rho = .5) {
  x3 <- seq(-.9995, .9995, length.out = 1000)
  out <- 2 * sin(pi * (x[1])) - 2 * ((x[2] - .5)^2 - const) +
         mean(apply(cbind(x[1], x[2], x3), 1, cond.dens) * 3 * exp(-(x[1] + x3)/2)) * 1.999
  return(out)
}

# g <- function(x) {exp(-x^2 * 4)}
# jkl <- f2.null - f1
# asdf <- U.red %*% t(U.red) %*% g(X[,3])
# 
# jkl2 <- U.red %*% t(U.red) %*% jkl
# asdf2.fit <- fit_hal(Y = g(X[,3]), X = X[,1:2], yolo = FALSE)
# asdf2 <- predict(asdf2.fit, new_data = X[,1:2])

#################

# Monte Carlo approx of the improvement-in-fit
n <- 5000
rho <- .5
p <- 3
Sigma <- matrix(rho, p, p) + (1 - rho) * diag(p)
X <- 2 * (pnorm(rmvnorm(n, mean = rep(0, p), sigma = Sigma)) - .5)

# Construct a kernel matrix, measuring similarity between rows of X0
tau <- sqrt(.2/2)
weights <- c(tau * (apply(X[,-3], 2, sd))^(-1), sd(X[,3])^(-1))
K <- matrix(NA, nrow = n, ncol = n)
K <- kernelMatrix(rbfdot(sigma = .5), X %*% diag(weights))

K.red <- matrix(NA, nrow = n, ncol = n)
K.red <- kernelMatrix(rbfdot(sigma = .5), X[,-3])

# Obtain the eigen-decomposition of the kernel matrices
eigen.min <- 1e-03

# K.eigenvals <- eigen(K, only.values = TRUE)
# basis.dim <- max(which(K.eigenvals$values >= eigen.min))
K.svd <- svd(K, nu = 100, nv = 100)
U <- K.svd$u # eigen-vectors
Lam <- K.svd$d[1:100] # eigen-values

# Calculate the projections
f1 <- apply(X, 1, function(x) theta0.list[[1]](x, noise.sd = 0))
# f1 <- f1 - mean(f1)
f2 <- apply(X, 1, function(x) theta0.list[[2]](x, noise.sd = 0))
# f2 <- f2 - mean(f2)
f3 <- apply(X, 1, function(x) theta0.list[[3]](x, noise.sd = 0))
# f3 <- f3 - mean(f3)

f1.null <- apply(X, 1, function(x) {null.list[[1]](x)})
# f1.null <- f1.null - mean(f1.null)
f2.null <- apply(X, 1, function(x) {null.list[[2]](x)})
# f2.null <- f2.null - mean(f2.null)
f3.null <- apply(X, 1, function(x) {null.list[[3]](x)})
# f3.null <- f3.null - mean(f3.null)

coef.1 <- t(U) %*% (f1 - f1.null)
coef.2 <- t(U) %*% (f2 - f2.null)
coef.3 <- t(U) %*% (f3 - f3.null)

# Calculate the improvement-in-fit
true.iif <- c(var(f1 - f1.null), var(f2 - f2.null), var(f3 - f3.null))

# Calculate oracle smoothness
oracle.roughness <- c(100, 
                      sum(coef.2^2/Lam)/var(f2 - f2.null) * 1.05,
                      sum(coef.3^2/Lam)/var(f3 - f3.null) * 1.05)

# oracle.roughness <- c(200, 200, 200)

#####################
### Storing Output
#####################

p.vals.oracle <- array(NA, dim = c(no.reps, length(ns), 3))
p.vals.adap <- array(NA, dim = c(no.reps, length(ns), 3))
p.vals.split <- array(NA, dim = c(no.reps, length(ns), 3))

ests.oracle <- array(NA, dim = c(no.reps, length(ns), 3))
ests.adap <- array(NA, dim = c(no.reps, length(ns), 3))
ests.oracle.debias <- array(NA, dim = c(no.reps, length(ns), 3))
ests.adap.debias <- array(NA, dim = c(no.reps, length(ns), 3))
ests.split <- array(NA, dim = c(no.reps, length(ns), 3))

ci.lower.oracle <- array(NA, dim = c(no.reps, length(ns), 3))
ci.lower.adap <- array(NA, dim = c(no.reps, length(ns), 3))
ci.lower.oracle.debias <- array(NA, dim = c(no.reps, length(ns), 3))
ci.lower.adap.debias <- array(NA, dim = c(no.reps, length(ns), 3))
ci.lower.split <- array(NA, dim = c(no.reps, length(ns), 3))

ci.upper.oracle <- array(NA, dim = c(no.reps, length(ns), 3))
ci.upper.adap <- array(NA, dim = c(no.reps, length(ns), 3))
ci.upper.oracle.debias <- array(NA, dim = c(no.reps, length(ns), 3))
ci.upper.adap.debias <- array(NA, dim = c(no.reps, length(ns), 3))
ci.upper.split <- array(NA, dim = c(no.reps, length(ns), 3))

roughness.adap <- array(NA, dim = c(no.reps, length(ns), 3))

#####################################
### Running Simulation
#####################################

for(k in 1:3) {
  for(j in 1:length(ns)) {
    n <- ns[j] #sample size
    for(i in 1:no.reps) {
      # generate predictors
      p <- 3
      Sigma <- matrix(rho, 3, 3) + (1 - rho) * diag(3)
      X <- 2 * (pnorm(rmvnorm(n, mean = rep(0, 3), sigma = Sigma)) - .5)
  
      # generate outcome
      y <- apply(X, 1, function(x) theta0.list[[k]](x, noise.sd = noise.sd))
      
      # Vimp for X3
      vimp.oracle <- NULL
      vimp.adapt <- NULL
      vimp.split <- NULL
      try(vimp.oracle <- Vimp(y = y, x = X[,3], W = X[,c(1,2)],
                              roughness = oracle.roughness[3]))
      try(vimp.adapt <- Vimp(y = y, x = X[,3], W = X[,c(1,2)],
                             roughness = NULL))
      try(vimp.split <- VimpSplit(y = y, x =  X[,3], W = X[,c(1,2)]))
      
      # Store results
      if(!is.null(vimp.oracle)) {
        p.vals.oracle[i, j, k] <- vimp.oracle$p.val
        ests.oracle[i, j, k] <- vimp.oracle$estimate
        ci.lower.oracle[i, j, k] <- vimp.oracle$ci.lower
        ci.upper.oracle[i, j, k] <- vimp.oracle$ci.upper
        ests.oracle.debias[i, j, k] <- vimp.oracle$estimate.debias
        ci.lower.oracle.debias[i, j, k] <- vimp.oracle$ci.lower.debias
        ci.upper.oracle.debias[i, j, k] <- vimp.oracle$ci.upper.debias
      }
      
      if(!is.null(vimp.adapt)) {
        p.vals.adap[i, j, k] <- vimp.adapt$p.val
        ests.adap[i, j, k] <- vimp.adapt$estimate
        ci.lower.adap[i, j, k] <- vimp.adapt$ci.lower
        ci.upper.adap[i, j, k] <- vimp.adapt$ci.upper
        ests.adap.debias[i, j, k] <- vimp.adapt$estimate.debias
        ci.lower.adap.debias[i, j, k] <- vimp.adapt$ci.lower.debias
        ci.upper.adap.debias[i, j, k] <- vimp.adapt$ci.upper.debias
        roughness.adap[i, j, k] <- vimp.adapt$roughness
      }
      if(!is.null(vimp.split)) {
        p.vals.split[i, j, k] <- vimp.split$p.val
        ests.split[i, j, k] <- vimp.split$estimate
        ci.lower.split[i, j, k] <- vimp.split$ci.lower
        ci.upper.split[i, j, k] <- vimp.split$ci.upper
      }
      
    }
  }
}

keep <- list(p.vals.oracle = p.vals.oracle,
             p.vals.adap = p.vals.adap,
             p.vals.split = p.vals.split,
             ests.oracle = ests.oracle,
             ests.adap = ests.adap,
             ests.oracle.debias = ests.oracle.debias,
             ests.adap.debias = ests.adap.debias,
             ests.split = ests.split,
             ci.lower.oracle = ci.lower.oracle,
             ci.upper.oracle = ci.upper.oracle,
             ci.lower.adap = ci.lower.adap,
             ci.upper.adap = ci.upper.adap,
             ci.lower.split = ci.lower.split,
             ci.upper.split = ci.upper.split,
             ci.lower.oracle.debias = ci.lower.oracle.debias,
             ci.upper.oracle.debias = ci.upper.oracle.debias,
             ci.lower.adap.debias = ci.lower.adap.debias,
             ci.upper.adap.debias = ci.upper.adap.debias,
             roughness.adap = roughness.adap,
             true.iif = true.iif,
             oracle.roughness = oracle.roughness,
             seed = seed)

fname <- paste0("/home/ahudson/nplrt/results-rev/out-11-28-24/nplrt-results-rev3-", job.id, ".rda")
save(keep, file = fname)
