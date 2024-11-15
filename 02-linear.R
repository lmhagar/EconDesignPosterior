## BEGIN SETUP ##

## load necessary packages
require(qrng)
require(foreach)
require(doParallel)
require(doSNOW)

## function to calculate the logit of the posterior probability
## u is point from pseudorandom sequence, params are the draw from Psi_0
## or Psi_1, deltas is the interval c(delta_L, delta_U), n_val is the 
## sample size presently explored, hyper is a list of the hyperparameters, 
## and q is the constant for imbalanced sample size determination
## Note: we can use mappings between sufficient statistics and the unit
## hypercube for this example
logitP <- function(u, params, delta, n_val, q, 
                       hyper = list(rep(0,3), 0.01*diag(3), 1, 1)){
  
  ## get coefficient values for data generation
  beta0 <- params[1]
  beta1 <- params[2]
  beta2 <- params[3]
  
  ## get values to generate covariates and error terms
  mu2 <- params[4]
  sigma2 <- params[5]
  sigmae <- params[6]
  
  ## save hyperparameters
  mu0 <- hyper[[1]]
  lambda0 <- hyper[[2]]
  a0 <- hyper[[3]]
  b0 <- hyper[[4]]
  
  ## obtain sample sizes for groups A (n1) and B (n2)
  n1 <- q*n_val
  n2 <- n_val

  ## generate group-specific sample means for x2
  n <- n1 + n2
  sum_x1 <- n1
  x21_bar <- qnorm(u[1], mu2, sigma2/sqrt(n1))
  x22_bar <- qnorm(u[2], mu2, sigma2/sqrt(n2))
  
  ## generate sample variance for x2 (both groups)
  ss_x2 <- sigma2^2*qchisq(u[3], n1 + n2 - 1)
  
  ## use algebra to get sufficient statistics
  sum_x2 <- n1*x21_bar + n2*x22_bar
  sum_x1x2 <- n1*x21_bar
  sum_x22 <- ss_x2 + sum_x2^2/n
  
  ## generate group-specific sample means for error terms
  eps21_bar <- qnorm(u[4], 0, sigmae/sqrt(n1))
  eps22_bar <- qnorm(u[5], 0, sigmae/sqrt(n2))
  
  ## use algebra to get sufficient statistics
  sum_eps <- n1*eps21_bar + n2*eps22_bar
  sum_x1eps <- n1*eps21_bar
  
  ## use Barlett decomposition to get remaining components of the 
  ## covariance matrix for x2 and varepsilon
  ss_epsx2 <- sigma2*sigmae*sqrt(qchisq(u[3], n1 + n2 - 1))*qnorm(u[6])
  ss_eps <- sigmae^2*(qchisq(u[7], n - 2) + qnorm(u[6])^2)
  sum_epsx2 <- ss_epsx2 + sum_x2*sum_eps/(n1 + n2)
  sum_eps2 <- ss_eps + sum_eps^2/n
  
  ## use algebra to get sufficient statistics
  sum_y <- beta0*(n1 + n2) + beta1*n1 + beta2*sum_x2 + sum_eps
  sum_yx1 <- beta0*n1 + beta1*n1 + beta2*n1*x21_bar + n1*eps21_bar
  sum_yx2 <- beta0*sum_x2 + beta1*sum_x1x2 + beta2*sum_x22 + sum_epsx2
  sum_y2 <- n*beta0^2 + beta1^2*n1 + beta2^2*sum_x22 + sum_eps2 +
    2*beta0*beta1*n1 + 2*beta0*beta2*sum_x2 + 2*beta0*sum_eps +
    2*beta1*beta2*sum_x1x2 + 2*beta1*sum_x1eps + 2*beta2*sum_epsx2
  
  ## compute the following statistics to get the location and
  ## scale parameter of the t-distribution posterior
  XtX <- rbind(c(n1 + n2, n1, sum_x2),
               c(n1, n1, sum_x1x2),
               c(sum_x2, sum_x1x2, sum_x22))
  
  Xty <- c(sum_y, sum_yx1, sum_yx2)
  lambdaN <- XtX + lambda0
  muN <- solve(lambdaN)%*%t(t(mu0)%*%lambda0 + Xty)
  
  aN <- a0 + 0.5*(n1 + n2)
  bN <- b0 + 0.5*(sum_y2 + t(mu0)%*%lambda0%*%mu0 - t(muN)%*%lambdaN%*%muN)
  
  ## compute posterior probability
  realP <- pt((delta[2]-muN[2])/sqrt(solve(lambdaN)[2,2]*bN/aN), 2*aN) - pt((delta[1]-muN[2])/sqrt(solve(lambdaN)[2,2]*bN/aN), 2*aN)
  
  ## slight perturbation of the posterior probability if it is too close to 0 or 1
  ## this ensures that the logit of the posterior probability is finite.
  if (realP > 1 - 10^(-7)){
    realP <- pt((delta[2]-muN[2])/sqrt(solve(lambdaN)[2,2]*bN/aN), 2*aN, lower.tail = FALSE) + pt((delta[1]-muN[2])/sqrt(solve(lambdaN)[2,2]*bN/aN), 2*aN)
    logitP <- log(realP) - log(1 - realP)
    logitP <- -1*logitP
  }
  else if (realP < .Machine$double.eps){
    realP <- .Machine$double.eps
    logitP <- log(realP) - log(1 - realP)
  }
  else{
    logitP <- log(realP) - log(1 - realP)
  }
  
  ## return the logit of the posterior probability
  return(logitP)
  
}

## Code to implement Algorithm 2 for the illustrative example. We explore sample sizes
## and return the optimal design (i.e., the (n, gamma) combination)
## Sobol' sequences can be used for more efficiency, but we use pseudorandom ones 
## in this paper, corresponding to independently generated samples
findGamma <- function(eta_plus1, eta_plus0, pwr, typeI, deltas, q, 
                      hypers, m = 8192, segments = 10,
                      seed_H1 = 1, seed_H0 = 2, prng = TRUE, boot = FALSE,
                      contour = FALSE, b1.low = 9, b1.high = 12,
                      n_integer = TRUE){
  
  ## eta_plus1: parameter draws from Psi_1
  ## eta_plus0: parameter draws from Psi_0
  ## pwr: target power (1 - beta)
  ## typeI: desired type I error rate (alpha)
  ## deltas: interval to define H1
  ## q: constant for imbalanced sample size determination
  ## hypers: hyperparamters for conjugate prior
  ## m: number of simulation repetitions
  ## segments: number of segments used to calculate linear approximations for assurance
  ## seed_H1: seed to generate Sobol' sequence for Psi_1
  ## seed_H0: seed to generate Sobol' sequence for Psi_0
  ## prng: TRUE is pseudorandom sequence is used (otherwise, use Sobol' sequence)
  ## boot: when TRUE, returns sampling distribution estimates for bootstrap
  ## contour: when TRUE, returns slopes and intercepts that construct contour plot
  ## b1.low: lower endpoint for uniform distribution for beta1 (specific to this Psi1)
  ## b1.high: upper endpoint for uniform distribution for beta1 (specific to this Psi1)
  ## n_integer: when FALSE, noninteger sample size recommendations are allowed (to check CI coverage)
  
  ## generate Sobol' sequence for each sampling distribution
  ## pseudorandom option is used by default
  if (prng == FALSE){
    sob_H1 <- sobol(m, d = 8, randomize = "digital.shift", seed = seed_H1)
    sob_H0 <- sobol(m, d = 7, randomize = "digital.shift", seed = seed_H0)
  }
  else{
    set.seed(seed_H1); sob_H1 <- matrix(runif(8*m), ncol = 8)
    set.seed(seed_H0); sob_H0 <- matrix(runif(7*m), ncol = 7)
  }
  
  ## generate effect sizes from the design prior
  eta_plus1[,2] <- b1.low + (b1.high - b1.low)*sob_H1[,8]
  
  ## uniroot without error checking
  uu <- function (fun, lower, upper, f_lower = NULL, f_upper = NULL, maxiter = 1000, tol = 1e-4, tol2 = 0.01, ...)
  {
    f <- function(x) fun(x, ...)
    x1 <- lower
    if (!is.null(f_lower)){
      f1 <- f_lower
    }
    else{
      f1 <- f(x1)
    }
    if (f1 > 0){return(x1)}
    x2 <- upper
    if (!is.null(f_upper)){
      f2 <- f_upper
    }
    else{
      f2 <- f(x2)
    }
    f2 <- f(x2)
    if (f2 < 0){return(x2)}
    x3 <- 0.5 * (lower + upper)
    niter <- 1
    while (niter <= maxiter) {
      f3 <- f(x3)
      if (abs(f3) < tol) {
        x0 <- x3
        return(x0)
      }
      if (f1 * f3 < 0) {
        upper <- x3}
      else {lower <- x3}
      if ((upper - lower) < tol2 * max(abs(upper), 1)) {
        x0 <- 0.5 * (lower + upper)
        return(x0)
      }
      denom <- (f2 - f1) * (f3 - f1) * (f2 - f3)
      numer <- x3 * (f1 - f2) * (f2 - f3 + f1) + f2 * x1 *
        (f2 - f3) + f1 * x2 * (f3 - f1)
      if (denom == 0) {
        dx <- upper - lower
      }
      else {
        dx <- f3 * numer/denom
      }
      x <- x3 + dx
      if ((upper - x) * (x - lower) < 0) {
        dx <- 0.5 * (upper - lower)
        x <- lower + dx
      }
      if (x1 < x3) {
        x2 <- x3
        f2 <- f3
      }
      else {
        x1 <- x3
        f1 <- f3
      }
      niter <- niter + 1
      if (abs(x - x3) < tol2) {
        x0 <- x
        return(x0)
      }
      x3 <- x
    }
    return(x0)
  }
  
  ## get initial sample size recommendation using asymptotic results
  EXtX <- rbind(c(q + 1, q, (q+1)*eta_plus1[1,4]),
                c(q, q, q*eta_plus1[1,4]),
                c((q+1)*eta_plus1[1,4], q*eta_plus1[1,4], (q+1)*(eta_plus1[1,5]^2 + eta_plus1[1,4]^2)))
  
  if (!is.finite(deltas[2])){
    n0 <- ceiling((qnorm(pwr) + qnorm(1 -typeI))^2*eta_plus1[1,6]^2*solve(EXtX)[2,2]/(0.5*(b1.high + b1.low)-deltas[1])^2)
  }
  else if (!is.finite(deltas[1])){
    n0 <- ceiling((qnorm(pwr) + qnorm(1 -typeI))^2*eta_plus1[1,6]^2*solve(EXtX)[2,2]/(0.5*(b1.high + b1.low)-deltas[2])^2)
  }
  else{
    ## for equivalence test, find more conservative upper bound using criteria for credible intervals
    theta0 <- 0.5*(b1.high + b1.low)
    a_cons <- (deltas[2] - theta0)/sqrt(eta_plus1[1,6]^2*solve(EXtX)[2,2])
    b_cons <- (deltas[1] - theta0)/sqrt(eta_plus1[1,6]^2*solve(EXtX)[2,2])
    c_cons <- qnorm(typeI/2)
    ## lower bound for root-finding algorithm
    lower_cons <- -2*c_cons*sqrt(eta_plus1[1,6]^2*solve(EXtX)[2,2])/(deltas[2] - deltas[1])
    upper_cons <- lower_cons
    
    fn_ci = function(n_sq, a, b, c, pwr){
      return(pnorm(a*n_sq + c) - pnorm(b*n_sq - c) - pwr)}
    
    upper_large <- FALSE
    while(upper_large == FALSE){
      upper_cons <- 10*upper_cons
      upper_check <- fn_ci(n_sq = sqrt(upper_cons), a = a_cons, b = b_cons, c = c_cons, pwr = pwr)
      if (upper_check > 0){
        upper_large <- TRUE
      }
    }
    
    upper_n <- (uu(fn_ci, a = a_cons, b = b_cons, c = c_cons, pwr = pwr, 
                   lower = lower_cons, upper = upper_cons))^2
    
    fn_start = function(n_trans, a, b, u, gam)
    {
      return(gam + pnorm(b, n_trans*qnorm(u), n_trans) - pnorm(a, n_trans*qnorm(u), n_trans))
    }
    uuu <- as.numeric(sobol(128, d = 1, randomize = "digital.shift", seed = seed_H1))
    start_rep <- NULL
    for (i in 1:length(uuu)){
      check_upper_n <- fn_start(n_trans = sqrt(eta_plus1[1,6]^2*solve(EXtX)[2,2])/sqrt(upper_n), 
                                a = deltas[2] - theta0, u = uuu[i],
                                b = deltas[1] - theta0, gam = 1 - typeI)
      if (check_upper_n > 0){
        start_rep[i] <- upper_n
      }
      else{
        start_rep[i] <- (uu(fn_start, a = deltas[2] - theta0, u = uuu[i],
                            b = deltas[1] - theta0, gam = 1 - typeI, 
                            lower = (sqrt(eta_plus1[1,6]^2*solve(EXtX)[2,2])/sqrt(upper_n)), upper = 1000)/sqrt(eta_plus1[1,6]^2*solve(EXtX)[2,2]))^(-2)
      }
    }
    n0 <- as.numeric(quantile(start_rep, pwr))
  }
  
  n <- n0
  
  H0_vec <- NULL
  H1_vec <- NULL
  for (i in 1:m){
    H1_vec[i] <- logitP(n_val = n, q = q,
                        params = as.numeric(eta_plus1[i,]),
                        delta = deltas, u = sob_H1[i,1:7],
                        hyper = hypers)
    H0_vec[i] <- logitP(n_val = n, q = q,
                        params = as.numeric(eta_plus0[i,]),
                        delta = deltas, u = sob_H0[i,],
                        hyper = hypers)
  }
  
  ## save the logits of the probabilities
  H1_vec_n0 <- H1_vec
  H0_vec_n0 <- H0_vec
  n0 <- n
  
  ## obtain the indices of the order statistics
  mid1 <- ceiling((1 - typeI)*m)
  mid2 <- floor((1-pwr)*m)
  
  ## if we do not satisfy the criterion on the order statistics, we create the linear
  ## approximation using a larger second sample size 
  if (sort(H0_vec)[mid1] > sort(H1_vec)[mid2]){
    H1_low <- H1_vec; H0_low <- H0_vec
    lower <- n
    n_low <- n
    ## choose a larger sample size n(1) using linear approximations
    lower_temp <- n0
    upper_temp <- 2*n0
    
    ## get the limiting slopes from Theorem 1
    slopes_H1_an <- 0.5*pmin((eta_plus1[,2] - deltas[1])^2, (eta_plus1[,2] - deltas[2])^2)/(eta_plus1[,6]^2*solve(EXtX)[2,2])
    
    while (upper_temp - lower_temp > 1){
      n_temp <- ceiling(0.5*(lower_temp + upper_temp))
      
      H1_vec_temp <- H1_vec + slopes_H1_an*(n_temp - lower)
      if (sort(H0_vec)[mid1] <= sort(H1_vec_temp)[mid2]){
        upper_temp <- n_temp
      }
      else{
        lower_temp <- n_temp
      }
    }
    
    ## estimate sampling distributions at second sample size
    if (prng == FALSE){
      sob_H1 <- sobol(m, d = 8, randomize = "digital.shift", seed = seed_H1 + 1000)
      sob_H0 <- sobol(m, d = 7, randomize = "digital.shift", seed = seed_H0 + 1000)
    }
    else{
      set.seed(seed_H1 + 1000); sob_H1 <- matrix(runif(8*m), ncol = 8)
      set.seed(seed_H0 + 1000); sob_H0 <- matrix(runif(7*m), ncol = 7)
    }
    
    eta_plus_n0 <- eta_plus1[,2]
    eta_plus1[,2] <- b1.low + (b1.high - b1.low)*sob_H1[,8]
    eta_plus_n1 <- eta_plus1[,2]
    
    n <- max(upper_temp, n0 + 5)
    H1_vec <- NULL
    H0_vec <- NULL
    ## use Algorithm 1 to approximate the posterior probabilities at all m points
    for (i in 1:m){
      H1_vec[i] <- logitP(n_val = n, q = q,
                          params = as.numeric(eta_plus1[i,]),
                          delta = deltas, u = sob_H1[i,1:7],
                          hyper = hypers)
      H0_vec[i] <- logitP(n_val = n, q = q,
                          params = as.numeric(eta_plus0[i,]),
                          delta = deltas, u = sob_H0[i,],
                          hyper = hypers)
    }
    
    ## save the logits of the posterior probabilities
    H1_vec_n1 <- H1_vec
    H0_vec_n1 <- H0_vec
    n1 <- n
    
    n_high <- n
    H1_high <- H1_vec; H0_high <- H0_vec
    
    ## allow for assurance under H1; break the posterior
    ## probabilities under H1 into segments based on their theta values
    ## before sorting
    m_seg <- rep(floor(m/segments), segments)
    if (m%%segments != 0){
      m_seg[seq(1, m%%segments, 1)] <- rep(ceiling(m/segments), m%%segments)
    }
    seg_start <- c(0,cumsum(m_seg))
    
    H1_slope <- NULL
    H1_int <- NULL
    for (k in 1:segments){
      H1_high_temp <- H1_high[which(rank(eta_plus_n1, ties.method = "first") %in% 
                                      seq(seg_start[k]+1, seg_start[k+1], 1))]
      H1_low_temp <- H1_low[which(rank(eta_plus_n0, ties.method = "first") %in% 
                                    seq(seg_start[k]+1, seg_start[k+1], 1))]
      H1_slope_temp <- (sort(H1_high_temp) - sort(H1_low_temp))/(n_high-n_low)
      H1_int_temp <- sort(H1_low_temp) - H1_slope_temp*n_low
      H1_slope <- c(H1_slope, H1_slope_temp)
      H1_int <- c(H1_int, H1_int_temp)
    }
    
    ## obtain the linear approximations for all m points
    H0_slope <- (sort(H0_high) - sort(H0_low))/(n_high-n_low)
    H0_int <- sort(H0_low) - H0_slope*n_low
    
    ## if the criterion on the order statistics is still not satisfied for this
    ## larger sample size, we need to find a larger upper bound for the binary search
    ## we do not do this using Algorithm 1 -- just the linear approximations to the
    ## posterior probabilities on the logit scale.
    if (sort(H0_vec)[mid1] <= sort(H1_vec)[mid2]){
      upper <- n
    }
    else {
      n <- ceiling(3*n_low)
      H1_vec <- H1_int + H1_slope*n
      H0_vec <- H0_int + H0_slope*n
      upper <- n
      
      while(sort(H0_vec)[mid1] > sort(H1_vec)[mid2]){
        n <- ceiling(1.5*n)
        H1_vec <- H1_int + H1_slope*n
        H0_vec <- H0_int + H0_slope*n
        if (sort(H0_vec)[mid1] <= sort(H1_vec)[mid2]){
          upper <- n
        }
      }
    }
  }
  ## if we satisfy the criterion on the order statistics, we create the linear
  ## approximation using a smaller second sample size 
  else{
    H1_high <- H1_vec; H0_high <- H0_vec
    upper <- n
    n_high <- n
    ## choose a smaller sample size n(1) using linear approximations
    upper_temp <- n0
    lower_temp <- floor(n0/2)
    
    ## get the limiting slopes from Theorem 1
    slopes_H1_an <- 0.5*pmin((eta_plus1[,2] - deltas[1])^2, (eta_plus1[,2] - deltas[2])^2)/(eta_plus1[,6]^2*solve(EXtX)[2,2])
    
    while (upper_temp - lower_temp > 1){
      n_temp <- ceiling(0.5*(lower_temp + upper_temp))
      
      H1_vec_temp <- H1_vec + slopes_H1_an*(n_temp - upper)
      if (sort(H0_vec)[mid1] <= sort(H1_vec_temp)[mid2]){
        upper_temp <- n_temp
      }
      else{
        lower_temp <- n_temp
      }
    }
    
    ## generate estimates of sampling distribution at second sample size
    if (prng == FALSE){
      sob_H1 <- sobol(m, d = 8, randomize = "digital.shift", seed = seed_H1 + 1000)
      sob_H0 <- sobol(m, d = 7, randomize = "digital.shift", seed = seed_H0 + 1000)
    }
    else{
      set.seed(seed_H1 + 1000); sob_H1 <- matrix(runif(8*m), ncol = 8)
      set.seed(seed_H0 + 1000); sob_H0 <- matrix(runif(7*m), ncol = 7)
    }
    
    eta_plus_n0 <- eta_plus1[,2]
    eta_plus1[,2] <- b1.low + (b1.high - b1.low)*sob_H1[,8]
    eta_plus_n1 <- eta_plus1[,2]
    
    n <- min(lower_temp, n0 - 5)
    H1_vec <- NULL
    H0_vec <- NULL
    ## use Algorithm 1 to approximate the posterior probabilities at all m points
    for (i in 1:m){
      H1_vec[i] <- logitP(n_val = n, q = q,
                          params = as.numeric(eta_plus1[i,]),
                          delta = deltas, u = sob_H1[i,1:7],
                          hyper = hypers)
      H0_vec[i] <- logitP(n_val = n, q = q,
                          params = as.numeric(eta_plus0[i,]),
                          delta = deltas, u = sob_H0[i,],
                          hyper = hypers)
    }
    
    ## save the logits of the posterior probabilities
    H1_vec_n1 <- H1_vec
    H0_vec_n1 <- H0_vec
    n1 <- n
    
    n_low <- n
    H1_low <- H1_vec; H0_low <- H0_vec
    
    ## allow for assurance under H1; break the posterior
    ## probabilities under H1 into segments based on their theta values
    ## before sorting
    m_seg <- rep(floor(m/segments), segments)
    if (m%%segments != 0){
      m_seg[seq(1, m%%segments, 1)] <- rep(ceiling(m/segments), m%%segments)
    }
    seg_start <- c(0,cumsum(m_seg))
    
    H1_slope <- NULL
    H1_int <- NULL
    for (k in 1:segments){
      H1_high_temp <- H1_high[which(rank(eta_plus_n0, ties.method = "first") %in% 
                                      seq(seg_start[k]+1, seg_start[k+1], 1))]
      H1_low_temp <- H1_low[which(rank(eta_plus_n1, ties.method = "first") %in% 
                                    seq(seg_start[k]+1, seg_start[k+1], 1))]
      H1_slope_temp <- (sort(H1_high_temp) - sort(H1_low_temp))/(n_high-n_low)
      H1_int_temp <- sort(H1_low_temp) - H1_slope_temp*n_low
      H1_slope <- c(H1_slope, H1_slope_temp)
      H1_int <- c(H1_int, H1_int_temp)
    }
    
    ## obtain the linear approximations for all m points
    H0_slope <- (sort(H0_high) - sort(H0_low))/(n_high-n_low)
    H0_int <- sort(H0_low) - H0_slope*n_low
    
    ## if the criterion on the order statistics is still satisfied for this
    ## larger sample size, we need to find a smaller lower bound for the binary search
    if (sort(H0_vec)[mid1] > sort(H1_vec)[mid2]){
      lower <- n
    }
    else{
      n <- floor(0.25*n_high)
      H1_vec <- H1_int + H1_slope*n
      H0_vec <- H0_int + H0_slope*n
      lower <- n
      
      while(sort(H0_vec)[mid1] <= sort(H1_vec)[mid2]){
        n <- floor(0.5*n)
        H1_vec <- H1_int + H1_slope*n
        H0_vec <- H0_int + H0_slope*n
        if (sort(H0_vec)[mid1] > sort(H1_vec)[mid2]){
          lower <- n
        }
      }
    }
  }
  
  if (n_integer == TRUE){
    tolerance <- 1
  }
  else {
    tolerance <- 0.01
  }
  ## implement binary search given these bounds for the sample size
  while ((upper - lower) > tolerance){
    if (n_integer == TRUE){
      n <- ceiling(0.5*(upper + lower))
    }
    else {
      n <- 0.5*(upper + lower)
    }
  
    ## use the linear approximations to estimate sampling 
    ## distributions for new sample sizes
    H1_vec <- H1_int + H1_slope*n
    H0_vec <- H0_int + H0_slope*n
    
    if (sort(H0_vec)[mid1] <= sort(H1_vec)[mid2]){
      upper <- n
      H1_final <- H1_vec; H0_final <- H0_vec
    }
    else{
      lower <- n
    }
  }
  
  n <- upper
  ## get the final estimates for the sampling distribution
  if (n == n_high){
    H1_final <- H1_high; H0_final <- H0_high
  }
  
  ## first component of list is (n, gamma, relevant order statistic of H1 probabilities, type I error at n2
  ## estimate, power estimate at n2, n0, n1)
  ## components 2 to 5 of the list are only needed for the contour plot (they are the slopes and intercepts of 
  ## the linear approximations under H0 and H1)
  ## components 6 to 11 are used to implement the bootstrap (they are the estimates of the sampling distributions
  ## and the effect sizes used to generate the data for assurance purposes)
  results <- list(c(n, 1/(1 + exp(-as.numeric(sort(H0_final)[mid1]))),
                    1/(1 + exp(-as.numeric(sort(H1_final)[mid2]))),
                    as.numeric(mean(H0_final > sort(H0_final)[mid1])), as.numeric(mean(H1_final > sort(H0_final)[mid1])),
                    sort(c(n0, n1))), H0_slope, H0_int, H1_slope, H1_int,
                    H0_low, H0_high, H1_low, H1_high, eta_plus_n0, eta_plus_n1)
  
  ## return only the necessary results
  if (contour == TRUE){
    if (boot == TRUE){
      return(results)
    }
    else {
      return(lapply(list(1, 2, 3, 4, 5), function(i){results[[i]]}))
    }
  }
  else if (boot == TRUE){
    return(lapply(list(1, 6, 7, 8, 9, 10, 11), function(i){results[[i]]}))
  }
  else{
    return(results[[1]][1:5])
  }
}

## set up hyperparameters and draws from Psi1 and Psi0 for numerical studies
## we overwrite the beta1 value in FindGamma()
eta_plus1s <- c(3 - 115*0.25, 10, 0.25, 115, 14.5, sqrt(10.7^2 - 0.25^2*(14.5^2)))
eta_plus1s <- matrix(rep(eta_plus1s, 10000), nrow = 10000, byrow = TRUE)

eta_plus0s <- c(3 - 115*0.25, 5, 0.25, 115, 14.5, sqrt(10.7^2 - 0.25^2*(14.5^2)))
eta_plus0s <- matrix(rep(eta_plus0s, 10000), nrow = 10000, byrow = TRUE)
hyper_ex <- list(rep(0,3), 0.01*diag(3), 1, 1)

## repeat the sample size calculation for the gamma example 1000 times
## with the different Sobol' sequences
## Use boot = TRUE to implement the bootstrap procedure that includes rounding
temp <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng"),
                .options.snow=opts, .errorhandling = "remove") %dopar% {
                  temp_res <- findGamma(eta_plus1 = eta_plus1s, eta_plus0 = eta_plus0s, 
                                        pwr = 0.8, typeI = 0.05, deltas = c(5, Inf), q = 2, 
                                        hypers = hyper_ex, m = 10000, boot = TRUE,
                                        seed_H1 = i, seed_H0 = i + 2000, contour = FALSE)
                  
                  H0_low_back <- temp_res[[2]]; H0_high_back <- temp_res[[3]]
                  H1_low_back <- temp_res[[4]]; H1_high_back <- temp_res[[5]]
                  
                  mm <- c(10000)
                  typeI <- 0.05; pwr <- 0.8
                  segments <- 10
                  set.seed(i + 4000)
                  
                  ## create matrix to save the results and update the indices of
                  ## the order statistics to reflect new number of resamples
                  res_matrix <- NULL
                  mid1 <- ceiling((1 - typeI)*mm)
                  mid2 <- floor((1-pwr)*mm)
                  
                  for (ii in 1:1000){
                    ## resample from each estimated sampling distribution
                    H1_high <- H1_high_back[sample(1:mm, mm, replace = TRUE)]
                    H1_low <- H1_low_back[sample(1:mm, mm, replace = TRUE)]
                    H0_high <- H0_high_back[sample(1:mm, mm, replace = TRUE)]
                    H0_low <- H0_low_back[sample(1:mm, mm, replace = TRUE)]
                    
                    ## create new linear approximations
                    m_seg <- rep(floor(mm/segments), segments)
                    if (m%%segments != 0){
                      m_seg[seq(1, mm%%segments, 1)] <- rep(ceiling(mm/segments), mm%%segments)
                    }
                    seg_start <- c(0,cumsum(m_seg))
                    
                    n_low <- temp_res[[1]][6]; n_high <- temp_res[[1]][7]
                    
                    H1_slope <- NULL
                    H1_int <- NULL
                    H1_high_temp <- NULL
                    H1_low_temp <- NULL
                    for (k in 1:segments){
                      H1_high_temp <- H1_high[which(rank(temp_res[[6]], ties.method = "first") %in%
                                                      seq(seg_start[k]+1, seg_start[k+1], 1))]
                      H1_low_temp <- H1_low[which(rank(temp_res[[7]], ties.method = "first") %in%
                                                    seq(seg_start[k]+1, seg_start[k+1], 1))]
                      H1_slope_temp <- (sort(H1_high_temp) - sort(H1_low_temp))/(n_high-n_low)
                      H1_int_temp <- sort(H1_low_temp) - H1_slope_temp*n_low
                      H1_slope <- c(H1_slope, H1_slope_temp)
                      H1_int <- c(H1_int, H1_int_temp)
                    }
                    
                    ## obtain the linear approximations for all m points
                    H0_slope <- (sort(H0_high) - sort(H0_low))/(n_high-n_low)
                    H0_int <- sort(H0_low) - H0_slope*n_low
                    
                    ## implement binary search with new linear approximations
                    lower <- floor(0.5*temp_res[[1]][6]); upper <- 2*temp_res[[1]][7]
                    while ((upper - lower) > 1){
                      n <- ceiling(0.5*(upper + lower))
                      
                      H1_vec <- H1_int + H1_slope*n
                      H0_vec <- H0_int + H0_slope*n
                      
                      if (sort(H0_vec)[mid1] <= sort(H1_vec)[mid2]){
                        upper <- n
                        H1_final <- H1_vec; H0_final <- H0_vec
                      } else{
                        lower <- n
                      }
                    }
                    
                    n <- upper
                    
                    ## save bootstrap results for each simulation repetition
                    res_matrix <- rbind(res_matrix, c(n, 1/(1 + exp(-as.numeric(sort(H0_final)[mid1]))),
                                                      1/(1 + exp(-as.numeric(sort(H1_final)[mid2]))),
                                                      as.numeric(mean(H0_final > sort(H0_final)[mid1])), as.numeric(mean(H1_final > sort(H0_final)[mid1])),
                                                      sort(c(n0, n1))))
                  }
                  
                  ## return results from Algorithm 2 with bootstrap CIs
                  c(i, temp_res[[1]],
                    quantile(res_matrix[,1], c(0.025, 0.975)), quantile(res_matrix[,2], c(0.025, 0.975)))
}

## output the optimal design and summary for each simulation repetition
write.csv(temp, "lin_boot_10k_summary.csv", row.names = FALSE)

## get 95% bootstrap CI for sample size n_B
quantile(temp[,2], c(0.025, 0.975))
## repeat for gamma
round(quantile(temp[,3], c(0.025, 0.975)),4)

## get estimated power and the type I error rate and the median 
## recommendations for n_B and gamma

## change the number of simulation repetitions to 4096
eta_plus1s <- c(3 - 115*0.25, 10, 0.25, 115, 14.5, sqrt(10.7^2 - 0.25^2*(14.5^2)))
eta_plus1s <- matrix(rep(eta_plus1s, 4096), nrow = 4096, byrow = TRUE)

eta_plus0s <- c(3 - 115*0.25, 5, 0.25, 115, 14.5, sqrt(10.7^2 - 0.25^2*(14.5^2)))
eta_plus0s <- matrix(rep(eta_plus0s, 4096), nrow = 4096, byrow = TRUE)
hyper_ex <- list(rep(0,3), 0.01*diag(3), 1, 1)

n <- 35
lA <- log(0.9564) - log(1 - 0.9564)
m <- 4096
deltas <- c(5, Inf)
temp35p95 <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng"),
                     .options.snow=opts, .errorhandling = "remove") %dopar% {
                       
                       sob_H1 <- sobol(m, d = 8, randomize = "digital.shift", seed = i)
                       sob_H0 <- sobol(m, d = 7, randomize = "digital.shift", seed = i + 1000)
                       
                       eta_plus1s[,2] <- 9 + (3)*sob_H1[,8]
                       
                       H1_vec <- NULL
                       H0_vec <- NULL
                       ## implement Algorithm 2 for all points
                       for (i in 1:nrow(sob_H1)){
                         H1_vec[i] <- logitP(n_val = n, q = 2,
                                             params = as.numeric(eta_plus1s[i,]),
                                             delta = deltas, u = sob_H1[i,1:7],
                                             hyper = hyper_ex)
                         H0_vec[i] <- logitP(n_val = n, q = 2,
                                             params = as.numeric(eta_plus0s[i,]),
                                             delta = deltas, u = sob_H0[i,],
                                             hyper = hyper_ex)
                       }
                       c(i, mean(H0_vec >= lA), mean(H1_vec >= lA))
                     }

write.csv(temp35p95, "lin_4096_35p95_summary.csv", row.names = FALSE)

## calculate type I error rate for this design
mean(temp35p95[,2])
## calculate power for this design
mean(temp35p95[,3])

## calculate coverage of the bootstrap confidence intervals
## when accounting for the rounding up

## coverage for sample size n_B
mean(temp[,9] <= median(temp[,2]) & temp[,10] >= median(temp[,2]))

## coverage for critical value gamma
mean(temp[,11] <= median(temp[,3]) & temp[,12] >= median(temp[,3]))

## change the number of simulation repetitions back to 10000
eta_plus1s <- c(3 - 115*0.25, 10, 0.25, 115, 14.5, sqrt(10.7^2 - 0.25^2*(14.5^2)))
eta_plus1s <- matrix(rep(eta_plus1s, 10000), nrow = 10000, byrow = TRUE)

eta_plus0s <- c(3 - 115*0.25, 5, 0.25, 115, 14.5, sqrt(10.7^2 - 0.25^2*(14.5^2)))
eta_plus0s <- matrix(rep(eta_plus0s, 10000), nrow = 10000, byrow = TRUE)
hyper_ex <- list(rep(0,3), 0.01*diag(3), 1, 1)

## repeat the sample size calculation for the gamma example 1000 times
## with the different Sobol' sequences
## Use boot = TRUE to implement the bootstrap procedure that allows noninteger n
tempFrac <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng"),
                .options.snow=opts, .errorhandling = "remove") %dopar% {
                  temp_res <- findGamma(eta_plus1 = eta_plus1s, eta_plus0 = eta_plus0s, 
                                        pwr = 0.8, typeI = 0.05, deltas = c(5, Inf), q = 2, 
                                        hypers = hyper_ex, m = 10000, boot = TRUE,
                                        seed_H1 = i, seed_H0 = i + 2000, contour = FALSE,
                                        n_integer = FALSE)
                  
                  H0_low_back <- temp_res[[2]]; H0_high_back <- temp_res[[3]]
                  H1_low_back <- temp_res[[4]]; H1_high_back <- temp_res[[5]]
                  
                  mm <- c(10000)
                  typeI <- 0.05; pwr <- 0.8
                  segments <- 10
                  set.seed(i + 4000)
                  
                  ## create matrix to save the results and update the indices of
                  ## the order statistics to reflect new number of resamples
                  res_matrix <- NULL
                  mid1 <- ceiling((1 - typeI)*mm)
                  mid2 <- floor((1-pwr)*mm)
                  
                  for (ii in 1:1000){
                    ## resample from each estimated sampling distribution
                    H1_high <- H1_high_back[sample(1:mm, mm, replace = TRUE)]
                    H1_low <- H1_low_back[sample(1:mm, mm, replace = TRUE)]
                    H0_high <- H0_high_back[sample(1:mm, mm, replace = TRUE)]
                    H0_low <- H0_low_back[sample(1:mm, mm, replace = TRUE)]
                    
                    ## create new linear approximations
                    m_seg <- rep(floor(mm/segments), segments)
                    if (m%%segments != 0){
                      m_seg[seq(1, mm%%segments, 1)] <- rep(ceiling(mm/segments), mm%%segments)
                    }
                    seg_start <- c(0,cumsum(m_seg))
                    
                    n_low <- temp_res[[1]][6]; n_high <- temp_res[[1]][7]
                    
                    H1_slope <- NULL
                    H1_int <- NULL
                    H1_high_temp <- NULL
                    H1_low_temp <- NULL
                    for (k in 1:segments){
                      H1_high_temp <- H1_high[which(rank(temp_res[[6]], ties.method = "first") %in%
                                                      seq(seg_start[k]+1, seg_start[k+1], 1))]
                      H1_low_temp <- H1_low[which(rank(temp_res[[7]], ties.method = "first") %in%
                                                    seq(seg_start[k]+1, seg_start[k+1], 1))]
                      H1_slope_temp <- (sort(H1_high_temp) - sort(H1_low_temp))/(n_high-n_low)
                      H1_int_temp <- sort(H1_low_temp) - H1_slope_temp*n_low
                      H1_slope <- c(H1_slope, H1_slope_temp)
                      H1_int <- c(H1_int, H1_int_temp)
                    }
                    
                    ## obtain the linear approximations for all m points
                    H0_slope <- (sort(H0_high) - sort(H0_low))/(n_high-n_low)
                    H0_int <- sort(H0_low) - H0_slope*n_low
                    
                    ## implement binary search with new linear approximations
                    lower <- floor(0.5*temp_res[[1]][6]); upper <- 2*temp_res[[1]][7]
                    while ((upper - lower) > 0.01){
                      n <- 0.5*(upper + lower)
                      
                      H1_vec <- H1_int + H1_slope*n
                      H0_vec <- H0_int + H0_slope*n
                      
                      if (sort(H0_vec)[mid1] <= sort(H1_vec)[mid2]){
                        upper <- n
                        H1_final <- H1_vec; H0_final <- H0_vec
                      } else{
                        lower <- n
                      }
                    }
                    
                    n <- upper
                    
                    ## save bootstrap results for each simulation repetition
                    res_matrix <- rbind(res_matrix, c(n, 1/(1 + exp(-as.numeric(sort(H0_final)[mid1]))),
                                                      1/(1 + exp(-as.numeric(sort(H1_final)[mid2]))),
                                                      as.numeric(mean(H0_final > sort(H0_final)[mid1])), as.numeric(mean(H1_final > sort(H0_final)[mid1])),
                                                      sort(c(n0, n1))))
                  }
                  
                  ## return results from Algorithm 2 with bootstrap CIs
                  c(i, temp_res[[1]],
                    quantile(res_matrix[,1], c(0.025, 0.975)), quantile(res_matrix[,2], c(0.025, 0.975)))
                }

## output the optimal design and summary for each simulation repetition
write.csv(tempFrac, "lin_boot_10k_frac.csv", row.names = FALSE)

## find median recommendations for sample size and critical value
median(tempFrac[,2]); median(tempFrac[,3])

## calculate coverage for new confidence intervals
## coverage for sample size n_B
mean(tempFrac[,9] <= median(tempFrac[,2]) & tempFrac[,10] >= median(tempFrac[,2]))

## coverage for critical value gamma
mean(tempFrac[,11] <= median(tempFrac[,3]) & tempFrac[,12] >= median(tempFrac[,3]))

## get estimated power and the type I error rate for the analytical 
## sample size recommendation with fixed gamma

## change the number of simulation repetitions to 4096
eta_plus1s <- c(3 - 115*0.25, 10, 0.25, 115, 14.5, sqrt(10.7^2 - 0.25^2*(14.5^2)))
eta_plus1s <- matrix(rep(eta_plus1s, 4096), nrow = 4096, byrow = TRUE)

eta_plus0s <- c(3 - 115*0.25, 5, 0.25, 115, 14.5, sqrt(10.7^2 - 0.25^2*(14.5^2)))
eta_plus0s <- matrix(rep(eta_plus0s, 4096), nrow = 4096, byrow = TRUE)
hyper_ex <- list(rep(0,3), 0.01*diag(3), 1, 1)

## define parameters for this setting
n <- 32
lA <- log(0.95) - log(0.05)
m <- 4096
deltas <- c(5, Inf)
temp32p95 <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng"),
                     .options.snow=opts, .errorhandling = "remove") %dopar% {
                       
                       sob_H1 <- sobol(m, d = 8, randomize = "digital.shift", seed = i)
                       sob_H0 <- sobol(m, d = 7, randomize = "digital.shift", seed = i + 1000)
                       
                       eta_plus1s[,2] <- 9 + (3)*sob_H1[,8]
                       
                       H1_vec <- NULL
                       H0_vec <- NULL
                       ## implement Algorithm 2 for all points
                       for (i in 1:nrow(sob_H1)){
                         H1_vec[i] <- logitP(n_val = n, q = 2,
                                             params = as.numeric(eta_plus1s[i,]),
                                             delta = deltas, u = sob_H1[i,1:7],
                                             hyper = hyper_ex)
                         H0_vec[i] <- logitP(n_val = n, q = 2,
                                             params = as.numeric(eta_plus0s[i,]),
                                             delta = deltas, u = sob_H0[i,],
                                             hyper = hyper_ex)
                       }
                       c(i, mean(H0_vec >= lA), mean(H1_vec >= lA))
                     }

write.csv(temp32p95, "lin_4096_32p95_summary.csv", row.names = FALSE)

## calculate type I error rate for this design
mean(temp32p95[,2])
## calculate power for this design
mean(temp32p95[,3])

## Code to implement the modified version of Algorithm 2 for the illustrative example.
## This version that finds n given fixed gamma.
findn <- function(eta_plus1, eta_plus0, pwr, typeI, deltas, q, 
                      hypers, m = 8192, segments = 10,
                      seed_H1 = 1, seed_H0 = 2, prng = TRUE, boot = FALSE,
                      contour = FALSE, b1.low = 9, b1.high = 12,
                      n_integer = TRUE){
  
  ## the inputs are the same as findGamma()
  
  ## however, the sampling distribution of posterior probabilities under H1 is always compared to
  ## the logit of typeI (for all sample sizes considered). Thus, the criterion for the type I error
  ## rate may not always be satisfied
  lA <- log(1 - typeI) - log(typeI)
  
  ## generate Sobol' sequence for each sampling distribution
  ## pseudorandom option is used by default
  if (prng == FALSE){
    sob_H1 <- sobol(m, d = 8, randomize = "digital.shift", seed = seed_H1)
    sob_H0 <- sobol(m, d = 7, randomize = "digital.shift", seed = seed_H0)
  }
  else{
    set.seed(seed_H1); sob_H1 <- matrix(runif(8*m), ncol = 8)
    set.seed(seed_H0); sob_H0 <- matrix(runif(7*m), ncol = 7)
  }
  
  ## generate effect sizes from the design prior
  eta_plus1[,2] <- b1.low + (b1.high - b1.low)*sob_H1[,8]
  
  ## uniroot without error checking
  uu <- function (fun, lower, upper, f_lower = NULL, f_upper = NULL, maxiter = 1000, tol = 1e-4, tol2 = 0.01, ...)
  {
    f <- function(x) fun(x, ...)
    x1 <- lower
    if (!is.null(f_lower)){
      f1 <- f_lower
    }
    else{
      f1 <- f(x1)
    }
    if (f1 > 0){return(x1)}
    x2 <- upper
    if (!is.null(f_upper)){
      f2 <- f_upper
    }
    else{
      f2 <- f(x2)
    }
    f2 <- f(x2)
    if (f2 < 0){return(x2)}
    x3 <- 0.5 * (lower + upper)
    niter <- 1
    while (niter <= maxiter) {
      f3 <- f(x3)
      if (abs(f3) < tol) {
        x0 <- x3
        return(x0)
      }
      if (f1 * f3 < 0) {
        upper <- x3}
      else {lower <- x3}
      if ((upper - lower) < tol2 * max(abs(upper), 1)) {
        x0 <- 0.5 * (lower + upper)
        return(x0)
      }
      denom <- (f2 - f1) * (f3 - f1) * (f2 - f3)
      numer <- x3 * (f1 - f2) * (f2 - f3 + f1) + f2 * x1 *
        (f2 - f3) + f1 * x2 * (f3 - f1)
      if (denom == 0) {
        dx <- upper - lower
      }
      else {
        dx <- f3 * numer/denom
      }
      x <- x3 + dx
      if ((upper - x) * (x - lower) < 0) {
        dx <- 0.5 * (upper - lower)
        x <- lower + dx
      }
      if (x1 < x3) {
        x2 <- x3
        f2 <- f3
      }
      else {
        x1 <- x3
        f1 <- f3
      }
      niter <- niter + 1
      if (abs(x - x3) < tol2) {
        x0 <- x
        return(x0)
      }
      x3 <- x
    }
    return(x0)
  }
  
  ## get initial sample size recommendation using asymptotic results
  EXtX <- rbind(c(q + 1, q, (q+1)*eta_plus1[1,4]),
                c(q, q, q*eta_plus1[1,4]),
                c((q+1)*eta_plus1[1,4], q*eta_plus1[1,4], (q+1)*(eta_plus1[1,5]^2 + eta_plus1[1,4]^2)))
  
  if (!is.finite(deltas[2])){
    n0 <- ceiling((qnorm(pwr) + qnorm(1 -typeI))^2*eta_plus1[1,6]^2*solve(EXtX)[2,2]/(0.5*(b1.high + b1.low)-deltas[1])^2)
  }
  else if (!is.finite(deltas[1])){
    n0 <- ceiling((qnorm(pwr) + qnorm(1 -typeI))^2*eta_plus1[1,6]^2*solve(EXtX)[2,2]/(0.5*(b1.high + b1.low)-deltas[2])^2)
  }
  else{
    ## for equivalence test, find more conservative upper bound using criteria for credible intervals
    theta0 <- 0.5*(b1.high + b1.low)
    a_cons <- (deltas[2] - theta0)/sqrt(eta_plus1[1,6]^2*solve(EXtX)[2,2])
    b_cons <- (deltas[1] - theta0)/sqrt(eta_plus1[1,6]^2*solve(EXtX)[2,2])
    c_cons <- qnorm(typeI/2)
    ## lower bound for root-finding algorithm
    lower_cons <- -2*c_cons*sqrt(eta_plus1[1,6]^2*solve(EXtX)[2,2])/(deltas[2] - deltas[1])
    upper_cons <- lower_cons
    
    fn_ci = function(n_sq, a, b, c, pwr){
      return(pnorm(a*n_sq + c) - pnorm(b*n_sq - c) - pwr)}
    
    upper_large <- FALSE
    while(upper_large == FALSE){
      upper_cons <- 10*upper_cons
      upper_check <- fn_ci(n_sq = sqrt(upper_cons), a = a_cons, b = b_cons, c = c_cons, pwr = pwr)
      if (upper_check > 0){
        upper_large <- TRUE
      }
    }
    
    upper_n <- (uu(fn_ci, a = a_cons, b = b_cons, c = c_cons, pwr = pwr, 
                   lower = lower_cons, upper = upper_cons))^2
    
    fn_start = function(n_trans, a, b, u, gam)
    {
      return(gam + pnorm(b, n_trans*qnorm(u), n_trans) - pnorm(a, n_trans*qnorm(u), n_trans))
    }
    uuu <- as.numeric(sobol(128, d = 1, randomize = "digital.shift", seed = seed_H1))
    start_rep <- NULL
    for (i in 1:length(uuu)){
      check_upper_n <- fn_start(n_trans = sqrt(eta_plus1[1,6]^2*solve(EXtX)[2,2])/sqrt(upper_n), 
                                a = deltas[2] - theta0, u = uuu[i],
                                b = deltas[1] - theta0, gam = 1 - typeI)
      if (check_upper_n > 0){
        start_rep[i] <- upper_n
      }
      else{
        start_rep[i] <- (uu(fn_start, a = deltas[2] - theta0, u = uuu[i],
                            b = deltas[1] - theta0, gam = 1 - typeI, 
                            lower = (sqrt(eta_plus1[1,6]^2*solve(EXtX)[2,2])/sqrt(upper_n)), upper = 1000)/sqrt(eta_plus1[1,6]^2*solve(EXtX)[2,2]))^(-2)
      }
    }
    n0 <- as.numeric(quantile(start_rep, pwr))
  }
  
  n <- n0
  
  H0_vec <- NULL
  H1_vec <- NULL
  for (i in 1:m){
    H1_vec[i] <- logitP(n_val = n, q = q,
                        params = as.numeric(eta_plus1[i,]),
                        delta = deltas, u = sob_H1[i,1:7],
                        hyper = hypers)
    H0_vec[i] <- logitP(n_val = n, q = q,
                        params = as.numeric(eta_plus0[i,]),
                        delta = deltas, u = sob_H0[i,],
                        hyper = hypers)
  }
  
  ## save the logits of the probabilities
  H1_vec_n0 <- H1_vec
  H0_vec_n0 <- H0_vec
  n0 <- n
  
  ## obtain the indices of the order statistics
  mid1 <- ceiling((1 - typeI)*m)
  mid2 <- floor((1-pwr)*m)
  
  ## if we do not satisfy the criterion on the order statistics, we create the linear
  ## approximation using a larger second sample size 
  if (lA > sort(H1_vec)[mid2]){
    H1_low <- H1_vec; H0_low <- H0_vec
    lower <- n
    n_low <- n
    ## choose a larger sample size n(1) using linear approximations
    lower_temp <- n0
    upper_temp <- 2*n0
    
    ## get the limiting slopes from Theorem 1
    slopes_H1_an <- 0.5*pmin((eta_plus1[,2] - deltas[1])^2, (eta_plus1[,2] - deltas[2])^2)/(eta_plus1[,6]^2*solve(EXtX)[2,2])
    
    while (upper_temp - lower_temp > 1){
      n_temp <- ceiling(0.5*(lower_temp + upper_temp))
      
      H1_vec_temp <- H1_vec + slopes_H1_an*(n_temp - lower)
      if (lA <= sort(H1_vec_temp)[mid2]){
        upper_temp <- n_temp
      }
      else{
        lower_temp <- n_temp
      }
    }
    
    ## estimate sampling distributions at second sample size
    if (prng == FALSE){
      sob_H1 <- sobol(m, d = 8, randomize = "digital.shift", seed = seed_H1 + 1000)
      sob_H0 <- sobol(m, d = 7, randomize = "digital.shift", seed = seed_H0 + 1000)
    }
    else{
      set.seed(seed_H1 + 1000); sob_H1 <- matrix(runif(8*m), ncol = 8)
      set.seed(seed_H0 + 1000); sob_H0 <- matrix(runif(7*m), ncol = 7)
    }
    
    eta_plus_n0 <- eta_plus1[,2]
    eta_plus1[,2] <- b1.low + (b1.high - b1.low)*sob_H1[,8]
    eta_plus_n1 <- eta_plus1[,2]
    
    n <- max(upper_temp, n0 + 5)
    H1_vec <- NULL
    H0_vec <- NULL
    ## use Algorithm 1 to approximate the posterior probabilities at all m points
    for (i in 1:m){
      H1_vec[i] <- logitP(n_val = n, q = q,
                          params = as.numeric(eta_plus1[i,]),
                          delta = deltas, u = sob_H1[i,1:7],
                          hyper = hypers)
      H0_vec[i] <- logitP(n_val = n, q = q,
                          params = as.numeric(eta_plus0[i,]),
                          delta = deltas, u = sob_H0[i,],
                          hyper = hypers)
    }
    
    ## save the logits of the posterior probabilities
    H1_vec_n1 <- H1_vec
    H0_vec_n1 <- H0_vec
    n1 <- n
    
    n_high <- n
    H1_high <- H1_vec; H0_high <- H0_vec
    
    ## allow for assurance under H1; break the posterior
    ## probabilities under H1 into segments based on their theta values
    ## before sorting
    m_seg <- rep(floor(m/segments), segments)
    if (m%%segments != 0){
      m_seg[seq(1, m%%segments, 1)] <- rep(ceiling(m/segments), m%%segments)
    }
    seg_start <- c(0,cumsum(m_seg))
    
    H1_slope <- NULL
    H1_int <- NULL
    for (k in 1:segments){
      H1_high_temp <- H1_high[which(rank(eta_plus_n1, ties.method = "first") %in% 
                                      seq(seg_start[k]+1, seg_start[k+1], 1))]
      H1_low_temp <- H1_low[which(rank(eta_plus_n0, ties.method = "first") %in% 
                                    seq(seg_start[k]+1, seg_start[k+1], 1))]
      H1_slope_temp <- (sort(H1_high_temp) - sort(H1_low_temp))/(n_high-n_low)
      H1_int_temp <- sort(H1_low_temp) - H1_slope_temp*n_low
      H1_slope <- c(H1_slope, H1_slope_temp)
      H1_int <- c(H1_int, H1_int_temp)
    }
    
    ## obtain the linear approximations for all m points
    H0_slope <- (sort(H0_high) - sort(H0_low))/(n_high-n_low)
    H0_int <- sort(H0_low) - H0_slope*n_low
    
    ## if the criterion on the order statistics is still not satisfied for this
    ## larger sample size, we need to find a larger upper bound for the binary search
    ## we do not do this using Algorithm 1 -- just the linear approximations to the
    ## posterior probabilities on the logit scale.
    if (lA <= sort(H1_vec)[mid2]){
      upper <- n
    }
    else {
      n <- ceiling(3*n_low)
      H1_vec <- H1_int + H1_slope*n
      H0_vec <- H0_int + H0_slope*n
      upper <- n
      
      while(lA > sort(H1_vec)[mid2]){
        n <- ceiling(1.5*n)
        H1_vec <- H1_int + H1_slope*n
        H0_vec <- H0_int + H0_slope*n
        if (lA <= sort(H1_vec)[mid2]){
          upper <- n
        }
      }
    }
  }
  ## if we satisfy the criterion on the order statistics, we create the linear
  ## approximation using a smaller second sample size 
  else{
    H1_high <- H1_vec; H0_high <- H0_vec
    upper <- n
    n_high <- n
    ## choose a smaller sample size n(1) using linear approximations
    upper_temp <- n0
    lower_temp <- floor(n0/2)
    
    ## get the limiting slopes from Theorem 1
    slopes_H1_an <- 0.5*pmin((eta_plus1[,2] - deltas[1])^2, (eta_plus1[,2] - deltas[2])^2)/(eta_plus1[,6]^2*solve(EXtX)[2,2])
    
    while (upper_temp - lower_temp > 1){
      n_temp <- ceiling(0.5*(lower_temp + upper_temp))
      
      H1_vec_temp <- H1_vec + slopes_H1_an*(n_temp - upper)
      if (lA <= sort(H1_vec_temp)[mid2]){
        upper_temp <- n_temp
      }
      else{
        lower_temp <- n_temp
      }
    }
    
    ## generate estimates of sampling distribution at second sample size
    if (prng == FALSE){
      sob_H1 <- sobol(m, d = 8, randomize = "digital.shift", seed = seed_H1 + 1000)
      sob_H0 <- sobol(m, d = 7, randomize = "digital.shift", seed = seed_H0 + 1000)
    }
    else{
      set.seed(seed_H1 + 1000); sob_H1 <- matrix(runif(8*m), ncol = 8)
      set.seed(seed_H0 + 1000); sob_H0 <- matrix(runif(7*m), ncol = 7)
    }
    
    eta_plus_n0 <- eta_plus1[,2]
    eta_plus1[,2] <- b1.low + (b1.high - b1.low)*sob_H1[,8]
    eta_plus_n1 <- eta_plus1[,2]
    
    n <- min(lower_temp, n0 - 5)
    H1_vec <- NULL
    H0_vec <- NULL
    ## use Algorithm 1 to approximate the posterior probabilities at all m points
    for (i in 1:m){
      H1_vec[i] <- logitP(n_val = n, q = q,
                          params = as.numeric(eta_plus1[i,]),
                          delta = deltas, u = sob_H1[i,1:7],
                          hyper = hypers)
      H0_vec[i] <- logitP(n_val = n, q = q,
                          params = as.numeric(eta_plus0[i,]),
                          delta = deltas, u = sob_H0[i,],
                          hyper = hypers)
    }
    
    ## save the logits of the posterior probabilities
    H1_vec_n1 <- H1_vec
    H0_vec_n1 <- H0_vec
    n1 <- n
    
    n_low <- n
    H1_low <- H1_vec; H0_low <- H0_vec
    
    ## allow for assurance under H1; break the posterior
    ## probabilities under H1 into segments based on their theta values
    ## before sorting
    m_seg <- rep(floor(m/segments), segments)
    if (m%%segments != 0){
      m_seg[seq(1, m%%segments, 1)] <- rep(ceiling(m/segments), m%%segments)
    }
    seg_start <- c(0,cumsum(m_seg))
    
    H1_slope <- NULL
    H1_int <- NULL
    for (k in 1:segments){
      H1_high_temp <- H1_high[which(rank(eta_plus_n0, ties.method = "first") %in% 
                                      seq(seg_start[k]+1, seg_start[k+1], 1))]
      H1_low_temp <- H1_low[which(rank(eta_plus_n1, ties.method = "first") %in% 
                                    seq(seg_start[k]+1, seg_start[k+1], 1))]
      H1_slope_temp <- (sort(H1_high_temp) - sort(H1_low_temp))/(n_high-n_low)
      H1_int_temp <- sort(H1_low_temp) - H1_slope_temp*n_low
      H1_slope <- c(H1_slope, H1_slope_temp)
      H1_int <- c(H1_int, H1_int_temp)
    }
    
    ## obtain the linear approximations for all m points
    H0_slope <- (sort(H0_high) - sort(H0_low))/(n_high-n_low)
    H0_int <- sort(H0_low) - H0_slope*n_low
    
    ## if the criterion on the order statistics is still satisfied for this
    ## larger sample size, we need to find a smaller lower bound for the binary search
    if (lA > sort(H1_vec)[mid2]){
      lower <- n
    }
    else{
      n <- floor(0.25*n_high)
      H1_vec <- H1_int + H1_slope*n
      H0_vec <- H0_int + H0_slope*n
      lower <- n
      
      while(lA <= sort(H1_vec)[mid2]){
        n <- floor(0.5*n)
        H1_vec <- H1_int + H1_slope*n
        H0_vec <- H0_int + H0_slope*n
        if (lA > sort(H1_vec)[mid2]){
          lower <- n
        }
      }
    }
  }
  
  if (n_integer == TRUE){
    tolerance <- 1
  }
  else {
    tolerance <- 0.01
  }
  ## implement binary search given these bounds for the sample size
  while ((upper - lower) > tolerance){
    if (n_integer == TRUE){
      n <- ceiling(0.5*(upper + lower))
    }
    else {
      n <- 0.5*(upper + lower)
    }
    
    ## use the linear approximations to estimate sampling 
    ## distributions for new sample sizes
    H1_vec <- H1_int + H1_slope*n
    H0_vec <- H0_int + H0_slope*n
    
    if (lA <= sort(H1_vec)[mid2]){
      upper <- n
      H1_final <- H1_vec; H0_final <- H0_vec
    }
    else{
      lower <- n
    }
  }
  
  n <- upper
  ## get the final estimates for the sampling distribution
  if (n == n_high){
    H1_final <- H1_high; H0_final <- H0_high
  }
  
  ## first component of list is (n, gamma, relevant order statistic of H1 probabilities, type I error at n2
  ## estimate, power estimate at n2, n0, n1)
  ## components 2 to 5 of the list are only needed for the contour plot (they are the slopes and intercepts of 
  ## the linear approximations under H0 and H1)
  ## components 6 to 11 are used to implement the bootstrap (they are the estimates of the sampling distributions
  ## and the effect sizes used to generate the data for assurance purposes)
  results <- list(c(n, 1/(1 + exp(-as.numeric(lA))),
                    1/(1 + exp(-as.numeric(sort(H1_final)[mid2]))),
                    as.numeric(mean(H0_final > lA)), as.numeric(mean(H1_final > lA)),
                    sort(c(n0, n1))), H0_slope, H0_int, H1_slope, H1_int,
                  H0_low, H0_high, H1_low, H1_high, eta_plus_n0, eta_plus_n1)
  
  ## return only the necessary results
  if (contour == TRUE){
    if (boot == TRUE){
      return(results)
    }
    else {
      return(lapply(list(1, 2, 3, 4, 5), function(i){results[[i]]}))
    }
  }
  else if (boot == TRUE){
    return(lapply(list(1, 6, 7, 8, 9, 10, 11), function(i){results[[i]]}))
  }
  else{
    return(results[[1]][1:5])
  }
}

## change the number of simulation repetitions back to 10000
eta_plus1s <- c(3 - 115*0.25, 10, 0.25, 115, 14.5, sqrt(10.7^2 - 0.25^2*(14.5^2)))
eta_plus1s <- matrix(rep(eta_plus1s, 10000), nrow = 10000, byrow = TRUE)

eta_plus0s <- c(3 - 115*0.25, 5, 0.25, 115, 14.5, sqrt(10.7^2 - 0.25^2*(14.5^2)))
eta_plus0s <- matrix(rep(eta_plus0s, 10000), nrow = 10000, byrow = TRUE)
hyper_ex <- list(rep(0,3), 0.01*diag(3), 1, 1)

## sample size calculation for fixed critical value gamma
tempFixed <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng"),
                     .options.snow=opts, .errorhandling = "remove") %dopar% {
                       temp_res <- findn(eta_plus1 = eta_plus1s, eta_plus0 = eta_plus0s, 
                                         pwr = 0.8, typeI = 0.05, deltas = c(5, Inf), q = 2, 
                                         hypers = hyper_ex, m = 10000,
                                         seed_H1 = i, seed_H0 = i + 1000, contour = FALSE)
                       
                       c(i, as.numeric(unlist(temp_res)))
                     }

write.csv(tempFixed, "lin_10k_fixed_summary.csv", row.names = FALSE)

median(tempFixed[,2]) ## recommended sample size is 33

## now calculate the operating characteristics for n = 33 and gamma = 0.95
## change the number of simulation repetitions to 4096
eta_plus1s <- c(3 - 115*0.25, 10, 0.25, 115, 14.5, sqrt(10.7^2 - 0.25^2*(14.5^2)))
eta_plus1s <- matrix(rep(eta_plus1s, 4096), nrow = 4096, byrow = TRUE)

eta_plus0s <- c(3 - 115*0.25, 5, 0.25, 115, 14.5, sqrt(10.7^2 - 0.25^2*(14.5^2)))
eta_plus0s <- matrix(rep(eta_plus0s, 4096), nrow = 4096, byrow = TRUE)
hyper_ex <- list(rep(0,3), 0.01*diag(3), 1, 1)

n <- 33
lA <- log(0.95) - log(0.05)
m <- 4096
deltas <- c(5, Inf)
temp33p95 <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng"),
                     .options.snow=opts, .errorhandling = "remove") %dopar% {
                       
                       sob_H1 <- sobol(m, d = 8, randomize = "digital.shift", seed = i)
                       sob_H0 <- sobol(m, d = 7, randomize = "digital.shift", seed = i + 1000)
                       
                       eta_plus1s[,2] <- 9 + (3)*sob_H1[,8]
                       
                       H1_vec <- NULL
                       H0_vec <- NULL
                       ## implement Algorithm 2 for all points
                       for (i in 1:nrow(sob_H1)){
                         H1_vec[i] <- logitP(n_val = n, q = 2,
                                             params = as.numeric(eta_plus1s[i,]),
                                             delta = deltas, u = sob_H1[i,1:7],
                                             hyper = hyper_ex)
                         H0_vec[i] <- logitP(n_val = n, q = 2,
                                             params = as.numeric(eta_plus0s[i,]),
                                             delta = deltas, u = sob_H0[i,],
                                             hyper = hyper_ex)
                       }
                       c(i, mean(H0_vec >= lA), mean(H1_vec >= lA))
                     }

write.csv(temp33p95, "lin_4096_35p95_summary.csv", row.names = FALSE)

## calculate type I error rate for this design
mean(temp33p95[,2])
## calculate power for this design
mean(temp33p95[,3])

## change the number of simulation repetitions back to 10000
eta_plus1s <- c(3 - 115*0.25, 10, 0.25, 115, 14.5, sqrt(10.7^2 - 0.25^2*(14.5^2)))
eta_plus1s <- matrix(rep(eta_plus1s, 10000), nrow = 10000, byrow = TRUE)

eta_plus0s <- c(3 - 115*0.25, 5, 0.25, 115, 14.5, sqrt(10.7^2 - 0.25^2*(14.5^2)))
eta_plus0s <- matrix(rep(eta_plus0s, 10000), nrow = 10000, byrow = TRUE)
hyper_ex <- list(rep(0,3), 0.01*diag(3), 1, 1)

## repeat the sample size calculation for the gamma example 1000 times
## with the different Sobol' sequences (with contour = TRUE to return
## the slopes and intercepts we need to create the contour plots)
temp <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng"),
                .options.snow=opts, .errorhandling = "remove") %dopar% {
                  temp_res <- findGamma(eta_plus1 = eta_plus1s, eta_plus0 = eta_plus0s, 
                                        pwr = 0.8, typeI = 0.05, deltas = c(5, Inf), q = 2, 
                                        hypers = hyper_ex, m = 10000,
                                        seed_H1 = i, seed_H0 = i + 2000, contour = TRUE)
                  
                  c(i, as.numeric(unlist(temp_res)))
                }

## output the optimal design and summary for each simulation repetition
write.csv(temp[,1:8], "lin_contour_10k_summary.csv", row.names = FALSE)

## output the red and green posterior probabilities for each of the three sample
## sizes n^((0)) = small, n^((1)) = mid, and n^((2)) = large
write.csv(temp[,9:10008], "lin_10k_H0_slope.csv", row.names = FALSE)
write.csv(temp[,10009:20008], "lin_10k_H0_int.csv", row.names = FALSE)
write.csv(temp[,20009:30008], "lin_10k_H1_slope.csv", row.names = FALSE)
write.csv(temp[,30009:40008], "lin_10k_H1_int.csv", row.names = FALSE)

## contour plot confirmatory simulations with 81920 repetitions
resamps <- 81920
pb <- txtProgressBar(max = resamps, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## get parameter values from Psi1 and Psi0 to generate data from the linear model
samps <- seq(25, 45, 1)

delta_L <- 5
delta_U <- Inf
deltas <- c(delta_L, delta_U)

eta_plus1s <- c(3 - 115*0.25, 10, 0.25, 115, 14.5, sqrt(10.7^2 - 0.25^2*(14.5^2)))
eta_plus1s <- matrix(rep(eta_plus1s, resamps), nrow = resamps, byrow = TRUE)
## overwrite the beta1 value
eta_plus1s[,2] <- runif(resamps)*3 + 9

eta_plus0s <- c(3 - 115*0.25, 5, 0.25, 115, 14.5, sqrt(10.7^2 - 0.25^2*(14.5^2)))
eta_plus0s <- matrix(rep(eta_plus0s, resamps), nrow = resamps, byrow = TRUE)
hyper_ex <- list(rep(0,3), 0.01*diag(3), 1, 1)

for (k in 1:length(samps)){
  tic <- Sys.time()
  probs_in <- foreach(j=1:resamps, .combine='rbind', .packages = c("rjags", "coda"),
                      .options.snow=opts, .errorhandling = "remove") %dopar% {
                        
                        set.seed(k*resamps + j)
                        
                        temp <- logitP(n_val = samps[k], q = q,
                                       params = eta_plus1s[j,],
                                       delta = deltas, u = runif(7),
                                       hyper = hyper_ex)
                        
                        1/(1 + exp(-temp))
                        
                      }
  toc <- Sys.time()
  toc - tic
  
  ## output results to .csv file
  write.csv(probs_in, paste0("probs_H1_lin_", samps[k], ".csv"), row.names = FALSE)
}

## now consider the process under H0
params <- rbind(eta_plus0s, eta_plus0s, eta_plus0s, eta_plus0s, eta_plus0s)

for (k in 1:length(samps)){
  
  probs_out <- foreach(j=1:resamps, .combine='rbind', .packages = c("rjags", "coda"),
                       .options.snow=opts, .errorhandling = "remove") %dopar% {
                         
                         set.seed(k*resamps + j)
                         
                         temp <- logitP(n_val = samps[k], q = q,
                                        params = eta_plus0s[j,],
                                        delta = deltas, u = runif(7),
                                        hyper = hyper_ex)
                         
                         1/(1 + exp(-temp))
                       }
  
  ## output results to .csv file
  write.csv(probs_out, paste0("probs_H0_lin_", samps[k], ".csv"), row.names = FALSE)
}

## code to create contour plots
## create matrices for left contour plot in Figure 1 (based on single
## sample size calculation)
first_rep <- as.numeric(read.csv("lin_contour_10k_summary.csv")[960,])
n_low <- first_rep[7]; n_high <- first_rep[8]

## read in the posterior probabilities corresponding to n0 and n1
H0_slope <- unlist(read.csv("lin_10k_H0_slope.csv")[960,])
H0_int <- unlist(read.csv("lin_10k_H0_int.csv")[960,])
H1_slope <- unlist(read.csv("lin_10k_H1_slope.csv")[960,])
H1_int <- unlist(read.csv("lin_10k_H1_int.csv")[960,])

## approximate the sampling distributions of posterior probabilities
## on the logit scale using these approximations
for (i in seq(25, 45, 1)){
  assign(paste0("H1_vec_", i), H1_int + H1_slope*i)
  assign(paste0("H0_vec_", i), H0_int + H0_slope*i)
}

## create a vector of gamma values on the logit scale to compute power
## and type I error rate estimates
opt_gamma <- as.numeric(first_rep[3])
opt_gamma <- log(opt_gamma) - log(1 - opt_gamma)

gammas <- seq(log(0.93) - log(0.07), log(0.97) - log(0.03), length.out = 50)
gammas <- sort(c(opt_gamma, gammas))

x <- seq(25, 45, 1)
y <- 1/(1 + exp(-gammas))

## z matrix is for power
z_mat <- NULL
for (i in 1:length(x)){
  z_mat <- rbind(z_mat, get(paste0("H1_vec_", x[i])))
}
z <- NULL
for (j in 1:length(y)){
  z <- cbind(z, rowMeans(z_mat > gammas[j]))
}

## w matrix is for type I error
w_mat <- NULL
for (i in 1:length(x)){
  w_mat <- rbind(w_mat, get(paste0("H0_vec_", x[i])))
}
w <- NULL
for (j in 1:length(y)){
  w <- cbind(w, rowMeans(w_mat > gammas[j]))
}

write.csv(w, "w_mat1.csv", row.names = FALSE)
write.csv(z, "z_mat1.csv", row.names = FALSE)

## create matrices for center contour plot in Figure 1 (based on the average
## of 1000 sample size calculations)

## read in the sample sizes and posterior probabilities (logit scale) for 
## all 1000 repetitions
n_lows <- as.numeric(read.csv("lin_contour_10k_summary.csv")[,7])
n_highs <- as.numeric(read.csv("lin_contour_10k_summary.csv")[,8])

H0_slopes <- read.csv("lin_10k_H0_slope.csv")
H0_ints <- read.csv("lin_10k_H0_int.csv")
H1_slopes <- read.csv("lin_10k_H1_slope.csv")
H1_ints <- read.csv("lin_10k_H1_int.csv")

z_full <- matrix(0, nrow = 21, ncol = 50)
w_full <- matrix(0, nrow = 21, ncol = 50)
for (k in 1:1000){
  ## the process below repeats the process detailed for sample size calculation
  ## 1 above for all repetitions k = {1, 2, ..., 1000}
  print(k)
  H1_slope <- as.numeric(H1_slopes[k,])
  H0_slope <- as.numeric(H0_slopes[k,])
  H1_int <- as.numeric(H1_ints[k,])
  H0_int <- as.numeric(H0_ints[k,])
  
  for (i in seq(25, 45, 1)){
    assign(paste0("H1_vec_", i), H1_int + H1_slope*i)
    assign(paste0("H0_vec_", i), H0_int + H0_slope*i)
  }
  
  gammas <- seq(log(0.93) - log(0.07), log(0.97) - log(0.03), length.out = 50)
  
  x <- seq(25, 45, 1)
  y <- 1/(1 + exp(-gammas))
  
  z_mat <- NULL
  for (i in 1:length(x)){
    z_mat <- rbind(z_mat, get(paste0("H1_vec_", x[i])))
  }
  z <- NULL
  for (j in 1:length(y)){
    z <- cbind(z, rowMeans(z_mat > gammas[j]))
  }
  
  ## now we multiply the matrix for each sample size calculation by 1/1000
  ## to get the matrices corresponding to the average
  z_full <- 0.001*z + z_full
  
  w_mat <- NULL
  for (i in 1:length(x)){
    w_mat <- rbind(w_mat, get(paste0("H0_vec_", x[i])))
  }
  w <- NULL
  for (j in 1:length(y)){
    w <- cbind(w, rowMeans(w_mat > gammas[j]))
  }
  
  w_full <- 0.001*w + w_full
  
}
write.csv(z_full, "z_full_lin.csv", row.names = FALSE)
write.csv(w_full, "w_full_lin.csv", row.names = FALSE)

## create matrices for right contour plot in Figure 1 (based on
## simulating gamma data)
z_full2 <- matrix(0, nrow = 21, ncol = 50)
w_full2 <- matrix(0, nrow = 21, ncol = 50)
## convert the posterior probabilities for each approximated sampling
## distribution to the logit scale (with error checking to ensure
## to logits are finite)
for (i in seq(25, 45,2)){
  assign(paste0("H1_vec_", i), 
         pmin(pmax(as.numeric(unlist(read.csv(paste0("probs_H1_gamma_", i,".csv")))),
                   .Machine$double.eps), 1 - 10^(-7)))
  assign(paste0("H0_vec_", i), 
         pmin(pmax(as.numeric(unlist(read.csv(paste0("probs_H0_gamma_", i,".csv")))),
                   .Machine$double.eps), 1 - 10^(-7)))
}

## this process mirrors what was done to create the z and w matrices in 
## the previous two plots but with the estimates obtained by simulating data
gammas <- seq(log(0.93) - log(0.07), log(0.97) - log(0.03), length.out = 50)

x <- seq(275, 325, 2)
y <- 1/(1 + exp(-gammas))

z_mat <- NULL
for (i in 1:length(x)){
  z_mat <- rbind(z_mat, get(paste0("H1_vec_", x[i])))
}
z_mat <- log(z_mat) - log(1-z_mat)
z <- NULL
for (j in 1:length(y)){
  z <- cbind(z, rowMeans(z_mat > gammas[j]))
}

z_full2 <- z

w_mat <- NULL
for (i in 1:length(x)){
  w_mat <- rbind(w_mat, get(paste0("H0_vec_", x[i])))
}
w_mat <- log(w_mat) - log(1-w_mat)
w <- NULL
for (j in 1:length(y)){
  w <- cbind(w, rowMeans(w_mat > gammas[j]))
}

w_full2 <- w

## write output to a .csv file  
write.csv(z_full2, "z_full_mat2.csv", row.names = FALSE)
write.csv(w_full2, "w_full_mat2.csv", row.names = FALSE)

## create the three contour plots and output as .pdf file for the article
pdf(file = "Figure1.pdf",   # The directory you want to save the file in
    width = 7.5, 
    height = 5) 

par(mfrow=c(2,3), mar = c(3.75, 3.75, 2, 0.35) + 0.1, mgp=c(2.25,1,0))

## read in matrices for left plot
z <- matrix(unlist(read.csv("z_mat1.csv")), nrow = 21, ncol = 51)
w <- matrix(unlist(read.csv("w_mat1.csv")), nrow =21, ncol = 51)
gammas <- seq(log(0.93) - log(0.07), log(0.97) - log(0.03), length.out = 50)
gammas <- sort(c(opt_gamma, gammas))
y <- 1/(1 + exp(-gammas))
x <- seq(25, 45, 1)

contour(x, y, w, levels = c(seq(0.02, 0.030, 0.005), 0.045, seq(0.065, 0.08, 0.005)), 
        xlab = expression(italic("n")), xlim = c(25,45),
        ylab = expression(gamma),  main = "Type I Error Rate", labcex = 0.8, method = "edge", axes = FALSE, cex.lab = 1.25)
contour(x, y, w, levels = c(0.05), col = "firebrick", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, z, levels = c(0.8), col = "seagreen", add = TRUE, labcex = 0.8, labels = "", method = "edge")
contour(x, y, w, levels = c(0.055), col = "black", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, w, levels = c(0.06), col = "black", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, w, levels = c(0.04), col = "black", add = TRUE, labcex = 0.8, method = "edge")
points(x = first_rep[2], y = 1/(1 + exp(-opt_gamma)), pch = 19, col = adjustcolor("grey50", 0.75))
contour(x, y, w, levels = c(0.035), col = "black", add = TRUE, labcex = 0.8, method = "edge", labels = " ")
axis(side = 1, at = seq(25, 45, 5), cex.axis = 1.15)
axis(side = 2, at = seq(0.93, 0.97, 0.01), cex.axis = 1.15)
box()

gammas <- seq(log(0.93) - log(0.07), log(0.97) - log(0.03), length.out = 50)
y <- 1/(1 + exp(-gammas))

z_full <- matrix(unlist(read.csv("z_full_lin.csv")), nrow = 21, ncol = 50)
w_full <- matrix(unlist(read.csv("w_full_lin.csv")), nrow = 21, ncol = 50)

contour(x, y, w_full, levels = c(seq(0.02, 0.035, 0.005), 0.045, seq(0.065, 0.08, 0.005)), 
        xlab = expression(italic("n")), xlim = c(25,45),
        ylab = expression(gamma),  main = "Type I Error Rate", labcex = 0.8, method = "edge", axes = FALSE, cex.lab = 1.25)
contour(x, y, w_full, levels = c(0.05), col = "firebrick", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, z_full, levels = c(0.8), col = "seagreen", add = TRUE, labcex = 0.8, labels = "", method = "edge")
contour(x, y, w_full, levels = c(0.055), col = "black", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, w_full, levels = c(0.06), col = "black", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, w_full, levels = c(0.04), col = "black", add = TRUE, labcex = 0.8, method = "edge")
axis(side = 1, at = seq(25, 45, 5), cex.axis = 1.15)
axis(side = 2, at = seq(0.93, 0.97, 0.01), cex.axis = 1.15)
box()

z_full2 <- matrix(unlist(read.csv("z_full_lin2.csv")), nrow = 21, ncol = 50)
w_full2 <- matrix(unlist(read.csv("w_full_lin2.csv")), nrow = 21, ncol = 50)

contour(x, y, w_full2, levels = c(seq(0.02, 0.035, 0.005), 0.045, seq(0.065, 0.08, 0.005)), 
        xlab = expression(italic("n")), xlim = c(25,45),
        ylab = expression(gamma),  main = "Type I Error Rate", labcex = 0.8, method = "edge", axes = FALSE, cex.lab = 1.25)
contour(x, y, w_full2, levels = c(0.05), col = "firebrick", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, z_full2, levels = c(0.8), col = "seagreen", add = TRUE, labcex = 0.8, labels = "", method = "edge")
contour(x, y, w_full2, levels = c(0.055), col = "black", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, w_full2, levels = c(0.06), col = "black", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, w_full2, levels = c(0.04), col = "black", add = TRUE, labcex = 0.8, method = "edge")
axis(side = 1, at = seq(25, 45, 5), cex.axis = 1.15)
axis(side = 2, at = seq(0.93, 0.97, 0.01), cex.axis = 1.15)
box()

gammas <- seq(log(0.93) - log(0.07), log(0.97) - log(0.03), length.out = 50)
gammas <- sort(c(opt_gamma, gammas))
y <- 1/(1 + exp(-gammas))
x <- seq(25, 45, 1)

contour(x, y, z, levels = c(0.68, 0.725, 0.77, 0.83, 0.845, 0.875), 
        xlab = expression(italic("n")), xlim = c(25,45),axes = FALSE, cex.lab = 1.25,
        ylab = expression(gamma),  main = "Power", labcex = 0.8, method = "edge")
contour(x, y, w, levels = c(0.05), col = "firebrick", add = TRUE, labcex = 0.8,labels = "", method = "edge") 
contour(x, y, z, levels = c(0.8), col = "seagreen", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, z, levels = c(0.65, 0.665, 0.695, 0.71, 0.74, 0.755, 0.785, 0.815, 0.86, 0.89, 0.905), 
        col = "black", add = TRUE, method = "edge", labcex = 0.8, labels = " ")
contour(x, y, z, levels = c(0.785), col = "black", add = TRUE, method = "edge", labcex = 0.8, labels = " ")
axis(side = 1, at = seq(25,45, 5), cex.axis = 1.15)
axis(side = 2, at = seq(0.93, 0.97, 0.01), cex.axis = 1.15)
box()
points(x = first_rep[2], y = 1/(1 + exp(-opt_gamma)), pch = 19, col = adjustcolor("grey50", 0.75))

gammas <- seq(log(0.93) - log(0.07), log(0.97) - log(0.03), length.out = 50)
y <- 1/(1 + exp(-gammas))

contour(x, y, z_full, levels = c(seq(0.65, 0.77, 0.015), 0.80, seq(0.83, 0.92, 0.015)), 
        xlab = expression(italic("n")), xlim = c(25,45),axes = FALSE, cex.lab = 1.25,
        ylab = expression(gamma),  main = "Power", labcex = 0.8, method = "edge")
contour(x, y, w_full, levels = c(0.05), col = "firebrick", add = TRUE, labcex = 0.8,labels = "", method = "edge") 
contour(x, y, z_full, levels = c(0.8), col = "seagreen", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, z_full, levels = c(0.815), col = "black", add = TRUE, method = "edge", labcex = 0.8, labels = " ")
contour(x, y, z_full, levels = c(0.785), col = "black", add = TRUE, method = "edge", labcex = 0.8, labels = " ")
axis(side = 1, at = seq(25,45, 5), cex.axis = 1.15)
axis(side = 2, at = seq(0.93, 0.97, 0.01), cex.axis = 1.15)
box()

contour(x, y, z_full2, levels = c(seq(0.65, 0.77, 0.015), 0.80, seq(0.83, 0.92, 0.015)), 
        xlab = expression(italic("n")), xlim = c(25,45),axes = FALSE, cex.lab = 1.25,
        ylab = expression(gamma),  main = "Power", labcex = 0.8, method = "edge")
contour(x, y, w_full2, levels = c(0.05), col = "firebrick", add = TRUE, labcex = 0.8,labels = "", method = "edge") 
contour(x, y, z_full2, levels = c(0.8), col = "seagreen", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, z_full2, levels = c(0.815), col = "black", add = TRUE, method = "edge", labcex = 0.8, labels = " ")
contour(x, y, z_full2, levels = c(0.785), col = "black", add = TRUE, method = "edge", labcex = 0.8, labels = " ")
axis(side = 1, at = seq(25,45, 5), cex.axis = 1.15)
axis(side = 2, at = seq(0.93, 0.97, 0.01), cex.axis = 1.15)
box()

par(mfrow=c(1,1))
dev.off()