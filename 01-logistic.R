## BEGIN SETUP ##

## load necessary packages

require(qrng)
require(foreach)
require(doParallel)
require(doSNOW)
require(rjags)
require(coda)
require(R.utils)

## define regression coefficients for Psi_1
eta_plus1 = c(-2.71, log(1.25), 0.25)
## regression coefficients for Psi_0
eta_plus0 = c(-2.71, log(2), 0.25) 
## define criteria for power and the type I error rate
pwr = 0.75
typeI = 0.4 
## define interval hypothesis and q
deltas = c(-Inf, log(2))
q = 2 

## get a vector of hyperparameters for the regression coefficients
## the mean and precision for beta_0's normal prior are the first two
## elements, followed by the normal parameters for beta_1 and beta_2
hypers = c(-2.71, 1, 0, 0.01, 0, 0.01) 

## mcmc settings consist of the number of chains, number of burnin iterations,
## the number of retained draws, and the thinning parameter
mcmc_settings = c(1, 500, 1000, 1)

## define these functions for later use
expit <- function(x){1/(1 + exp(-x))}
logit <- function(x){log(x) - log(1-x)}

## find initial sample size using code that is specific to this example

## use a Sobol' sequence to implement numerical integration
sob <- sobol(32768, 2, seed = 1, randomize = "digital.shift")

## generate a large low-discrepancy sample where there is a 2/3 probability
## of someone being assigned to the treatment (x1 covariate) and the covariate
## x2 follows a standard normal distribution (after centering and scaling)
x1 <- qbinom(sob[,1], 1, q/(q+1))
x2 <- qnorm(sob[,2])

## take the coefficients from Psi_1 to calculate sample size based on power
beta0 <- eta_plus1[1]
beta1 <- eta_plus1[2]
beta2 <- eta_plus1[3]

## estimate the components of the Fisher information matrix using empirical averages
W00 <- mean(expit(beta0 + beta1*x1 + beta2*x2)*(1-expit(beta0 + beta1*x1 + beta2*x2)))
W01 <- mean(x1*expit(beta0 + beta1*x1 + beta2*x2)*(1-expit(beta0 + beta1*x1 + beta2*x2)))
W02 <- mean(x2*expit(beta0 + beta1*x1 + beta2*x2)*(1-expit(beta0 + beta1*x1 + beta2*x2)))

W11 <- mean(x1^2*expit(beta0 + beta1*x1 + beta2*x2)*(1-expit(beta0 + beta1*x1 + beta2*x2)))
W12 <- mean(x1*x2*expit(beta0 + beta1*x1 + beta2*x2)*(1-expit(beta0 + beta1*x1 + beta2*x2)))
W22 <- mean(x2^2*expit(beta0 + beta1*x1 + beta2*x2)*(1-expit(beta0 + beta1*x1 + beta2*x2)))

W_mat <- rbind(c(W00, W01, W02), c(W01, W11, W12),
               c(W02, W12, W22))

## this equation returns the initial sample size for group B
n0 <- ceiling((qnorm(pwr) + qnorm(1 -typeI))^2*solve(W_mat)[2,2]/(beta1-deltas[2])^2/(q+1))
n <- n0

## set up parallelization with 10000 simulation repetitions
cores=detectCores()
cl <- makeSOCKcluster(cores[1]-1)

m <- 100000
registerDoSNOW(cl)
pb <- txtProgressBar(max = m, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## estimate sampling distribution under H1 at n0
H1_vec <- foreach(i=1:m, .packages=c('rjags', 'coda'), .combine=rbind,
                  .options.snow=opts) %dopar% {
                    log.wd <- paste(getwd(), '/JAGS_logistic.txt', sep='')
                    
                    beta0 <- eta_plus1[1]
                    beta1 <- eta_plus1[2]
                    beta2 <- eta_plus1[3]
                    
                    ## generate data
                    set.seed(i)
                    x1 <- rep(c(0,rep(1, q)), each = n)
                    x2 <- rnorm((q+1)*n)
                    eta <- beta0 + beta1*x1 + beta2*x2
                    y <- rbinom((q+1)*n, 1, expit(eta))
                    
                    ## extract MCMC settings
                    n.chains = mcmc_settings[1]
                    n.burnin = mcmc_settings[2]
                    n.draws = mcmc_settings[3]
                    n.thin = mcmc_settings[4]
                    
                    ## fit model and obtain posterior draws
                    model1.fit <- jags.model(file=log.wd,
                                             data=list(N=(q+1)*n, Y=y, X1=x1, X2=x2,
                                                       m0 = hypers[1], p0 = hypers[2],
                                                       m1 = hypers[2], p1 = hypers[4],
                                                       m2 = hypers[5], p2 = hypers[6]),
                                             n.chains = n.chains)
                    
                    update(model1.fit, n.burnin)
                   
                    model1.samples <- coda.samples(model1.fit, c("beta1"), 
                                                   n.iter=n.draws, thin=n.thin)
                    
                    beta1.post <- unlist(model1.samples[1])
                    
                    ## estimate posterior of beta1 using kernel density estimation      
                    kd <- density(beta1.post)
                    ## get the complementary of the probability H1 is true (more stable to
                    ## take the logit of this probability and multiply by -1)
                    comp.prob <- mean(pnorm(deltas[2], beta1.post, kd$bw, lower.tail = FALSE))
                          
                    res.i <- -1*logit(comp.prob)
                    c(i, res.i)
   
}
write.csv(H1_vec, paste0("H1_vec_", n,".csv"), row.names = FALSE)

## repeat the process above to estimate the sampling distribution under H0 at n1
H0_vec <- foreach(i=1:m, .packages=c('rjags', 'coda'), .combine=rbind,
                  .options.snow=opts) %dopar% {
                    log.wd <- paste(getwd(), '/JAGS_logistic.txt', sep='')
                    
                    beta0 <- eta_plus0[1]
                    beta1 <- eta_plus0[2]
                    beta2 <- eta_plus0[3]
                    
                    set.seed(m + i)
                    x1 <- rep(c(0,rep(1, q)), each = n)
                    x2 <- rnorm((q+1)*n)
                    eta <- beta0 + beta1*x1 + beta2*x2
                    y <- rbinom((q+1)*n, 1, expit(eta))
                    
                    n.chains = mcmc_settings[1]
                    n.burnin = mcmc_settings[2]
                    n.draws = mcmc_settings[3]
                    n.thin = mcmc_settings[4]
                    
                    model1.fit <- jags.model(file=log.wd,
                                             data=list(N=(q+1)*n, Y=y, X1=x1, X2=x2,
                                                       m0 = hypers[1], p0 = hypers[2],
                                                       m1 = hypers[2], p1 = hypers[4],
                                                       m2 = hypers[5], p2 = hypers[6]),
                                             n.chains = n.chains)
                    
                    update(model1.fit, n.burnin)
                    
                    model1.samples <- coda.samples(model1.fit, c("beta1"), 
                                                   n.iter=n.draws, thin=n.thin)
                    
                    beta1.post <- unlist(model1.samples[1])
                    
                    kd <- density(beta1.post)
                    comp.prob <- mean(pnorm(deltas[2], beta1.post, kd$bw, lower.tail = FALSE))
                    
                    res.i <- -1*logit(comp.prob)
                    c(i, res.i)
}
write.csv(H0_vec, paste0("H0_vec_", n,".csv"), row.names = FALSE)

## check the operating characteristics at this sample size
## save the logits of the probabilities to a vector
H1_vec_n0 <- H1_vec
H0_vec_n0 <- H0_vec

## we initialize the indices of the order statistics
mid1 <- ceiling((1 - typeI)*m) # index for xi_0
mid2 <- floor((1-pwr)*m) # index for xi_1

## show the sample size n0 is sufficiently large in this case
sort(H0_vec_n0[,2])[mid1] > sort(H1_vec_n0[,2])[mid2]

## initialize first round of binary search in the sample size space
H1_high <- as.numeric(H1_vec_n0[,2]); H0_high <- as.numeric(H0_vec_n0[,2])
upper <- n
n_high <- n

## choose a smaller sample size n_1 using linear approximations
upper_temp <- n0
lower_temp <- floor(n0/2)

## get the limiting slope from Theorem 1
slopes_H1_an <- (q+1)*0.5*(eta_plus1[2] - deltas[2])^2/solve(W_mat)[2,2]

while (upper_temp - lower_temp > 1){
  n_temp <- ceiling(0.5*(lower_temp + upper_temp))
  print(n_temp)
  
  ## approximate sampling distribution using linear approximations
  H1_vec_temp <- H1_vec_n0[,2] + slopes_H1_an*(n_temp - n0)
  
  ## check if xi_0 <= xi_1 (i.e., if sample size is large enough)
  if (sort(H0_vec_n0[,2])[mid1] <= sort(H1_vec_temp)[mid2]){
    upper_temp <- n_temp
  } else{
    lower_temp <- n_temp
  }
}

## extract the next sample size n1 from binary search
n <- lower_temp

## estimate sampling distribution under H1 at n1
H1_vec <- foreach(i=1:m, .packages=c('rjags', 'coda'), .combine=rbind,
                  .options.snow=opts) %dopar% {
                    log.wd <- paste(getwd(), '/JAGS_logistic.txt', sep='')
                    
                    beta0 <- eta_plus1[1]
                    beta1 <- eta_plus1[2]
                    beta2 <- eta_plus1[3]
                    
                    ## generate data
                    set.seed(2*m + i)
                    x1 <- rep(c(0,rep(1, q)), each = n)
                    x2 <- rnorm((q+1)*n)
                    eta <- beta0 + beta1*x1 + beta2*x2
                    y <- rbinom((q+1)*n, 1, expit(eta))
                    
                    ## extract MCMC settings
                    n.chains = mcmc_settings[1]
                    n.burnin = mcmc_settings[2]
                    n.draws = mcmc_settings[3]
                    n.thin = mcmc_settings[4]
                    
                    ## fit model and obtain posterior draws
                    model1.fit <- jags.model(file=log.wd,
                                             data=list(N=(q+1)*n, Y=y, X1=x1, X2=x2,
                                                       m0 = hypers[1], p0 = hypers[2],
                                                       m1 = hypers[2], p1 = hypers[4],
                                                       m2 = hypers[5], p2 = hypers[6]),
                                             n.chains = n.chains)
                    
                    update(model1.fit, n.burnin)
                    
                    model1.samples <- coda.samples(model1.fit, c("beta1"), 
                                                   n.iter=n.draws, thin=n.thin)
                    
                    beta1.post <- unlist(model1.samples[1])
                    
                    ## estimate posterior of beta1 using kernel density estimation      
                    kd <- density(beta1.post)
                    ## get the complementary of the probability H1 is true (more stable to
                    ## take the logit of this probability and multiply by -1)
                    comp.prob <- mean(pnorm(deltas[2], beta1.post, kd$bw, lower.tail = FALSE))
                    
                    res.i <- -1*logit(comp.prob)
                    c(i, res.i)
                    
                  }
write.csv(H1_vec, paste0("H1_vec_", n,".csv"), row.names = FALSE)

## repeat the process above to estimate the sampling distribution under H0 at n1
H0_vec <- foreach(i=1:m, .packages=c('rjags', 'coda'), .combine=rbind,
                  .options.snow=opts) %dopar% {
                    log.wd <- paste(getwd(), '/JAGS_logistic.txt', sep='')
                    
                    beta0 <- eta_plus0[1]
                    beta1 <- eta_plus0[2]
                    beta2 <- eta_plus0[3]
                    
                    set.seed(3*m + i)
                    x1 <- rep(c(0,rep(1, q)), each = n)
                    x2 <- rnorm((q+1)*n)
                    eta <- beta0 + beta1*x1 + beta2*x2
                    y <- rbinom((q+1)*n, 1, expit(eta))
                    
                    n.chains = mcmc_settings[1]
                    n.burnin = mcmc_settings[2]
                    n.draws = mcmc_settings[3]
                    n.thin = mcmc_settings[4]
                    
                    model1.fit <- jags.model(file=log.wd,
                                             data=list(N=(q+1)*n, Y=y, X1=x1, X2=x2,
                                                       m0 = hypers[1], p0 = hypers[2],
                                                       m1 = hypers[2], p1 = hypers[4],
                                                       m2 = hypers[5], p2 = hypers[6]),
                                             n.chains = n.chains)
                    
                    update(model1.fit, n.burnin)
                    
                    model1.samples <- coda.samples(model1.fit, c("beta1"), 
                                                   n.iter=n.draws, thin=n.thin)
                    
                    beta1.post <- unlist(model1.samples[1])
                    
                    kd <- density(beta1.post)
                    comp.prob <- mean(pnorm(deltas[2], beta1.post, kd$bw, lower.tail = FALSE))
                    
                    res.i <- -1*logit(comp.prob)
                    c(i, res.i)
                  }
write.csv(H0_vec, paste0("H0_vec_", n,".csv"), row.names = FALSE)

## check the operating characteristics at this sample size
## save the logits of the probabilities to a vector
H1_vec_n1 <- H1_vec
H0_vec_n1 <- H0_vec

## show that the sample size n1 is not large enough
sort(H0_vec_n1[,2])[mid1] > sort(H1_vec_n1[,2])[mid2]

## initialize binary search with better linear approximations
n1 <- n
n_low <- n
H1_low <- H1_vec_n1[,2]; H0_low <- H0_vec_n1[,2]

## obtain linear approximations using order statistics
H1_slope <- (sort(H1_high) - sort(H1_low))/(n_high-n_low)
H1_int <- sort(H1_low) - H1_slope*n_low

H0_slope <- (sort(H0_high) - sort(H0_low))/(n_high-n_low)
H0_int <- sort(H0_low) - H0_slope*n_low

lower <- n
## implement binary search given these bounds for the sample size
while ((upper - lower) > 1){
  n <- ceiling(0.5*(upper + lower))
  print(n)
  
  ## use the linear approximations to approximate sampling distributions
  ## at new sample sizes
  H1_vec <- H1_int + H1_slope*n
  H0_vec <- H0_int + H0_slope*n
  
  if (sort(H0_vec)[mid1] <= sort(H1_vec)[mid2]){
    upper <- n
    H1_final <- H1_vec; H0_final <- H0_vec
  } else{
    lower <- n
  }
}

## extract the final sample size recommendation
n <- upper

## save the results
results <- list(c(n, 1/(1 + exp(-as.numeric(sort(H0_final)[mid1]))),
                  1/(1 + exp(-as.numeric(sort(H1_final)[mid2]))),
                  as.numeric(mean(H0_final > sort(H0_final)[mid1])), 
                  as.numeric(mean(H1_final > sort(H0_final)[mid1])),
                  sort(c(n0, n1))), H0_slope, H0_int, H1_slope, H1_int)

## code to implement the bootstrap procedure
H1_high <- read.csv("H1_vec_93.csv")[,2]; H1_low <- read.csv("H1_vec_63.csv")[,2]
H0_high <- read.csv("H0_vec_93.csv")[,2]; H0_low <- read.csv("H0_vec_63.csv")[,2]

## save the vectors of estimated posterior probabilities to vectors that will
## be sampled from
H1_high_back <- H1_high; H1_low_back <- H1_low
H0_high_back <- H0_high; H0_low_back <- H0_low
n_high <- 93; n_low <- 63

## mm is m^* in the manuscript (the number of resamples)
mm <- c(5000, 10000, 25000, 50000, 75000, 100000)
for (kk in 1:length(mm)){
  set.seed(kk)
  
  ## create matrix to save the results and update the indices of 
  ## the order statistics to reflect new number of resamples
  res_matrix <- NULL
  mid1 <- ceiling((1 - typeI)*mm[kk])
  mid2 <- floor((1-pwr)*mm[kk])
  
  for (ii in 1:10000){
    ## resample from each estimated sampling distribution
    H1_high <- H1_high_back[sample(1:m, mm[kk], replace = TRUE)]
    H1_low <- H1_low_back[sample(1:m, mm[kk], replace = TRUE)]
    H0_high <- H0_high_back[sample(1:m, mm[kk], replace = TRUE)]
    H0_low <- H0_low_back[sample(1:m, mm[kk], replace = TRUE)]
  
    ## create new linear approximations
    H1_slope <- (sort(H1_high) - sort(H1_low))/(n_high-n_low)
    H1_int <- sort(H1_low) - H1_slope*n_low
  
    H0_slope <- (sort(H0_high) - sort(H0_low))/(n_high-n_low)
    H0_int <- sort(H0_low) - H0_slope*n_low
  
    ## implement binary search with new linear approximations
    lower <- 63; upper <- 93
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
  
    ## save results for each simulation repetition
    res_matrix <- rbind(res_matrix, c(n, 1/(1 + exp(-as.numeric(sort(H0_final)[mid1]))),
      1/(1 + exp(-as.numeric(sort(H1_final)[mid2]))),
      as.numeric(mean(H0_final > sort(H0_final)[mid1])), as.numeric(mean(H1_final > sort(H0_final)[mid1])),
      sort(c(n0, n1))))
  }
  ## write results to a .csv file
  write.csv(res_matrix, paste0("bootstrap_cis_", mm[kk], ".csv"), row.names = FALSE)
}

## construct table with the bootstrap confidence intervals (Table 1)
CI_n <- NULL
CI_gamma <- NULL
for (kk in 1:length(mm)){
  temp_dat <- read.csv(paste0("bootstrap_cis_", mm[kk], ".csv"))
  CI_n <- rbind(CI_n, quantile(temp_dat[,1], c(0.025, 0.975)))
  CI_gamma <- rbind(CI_gamma, quantile(temp_dat[,2], c(0.025, 0.975)))
}

nn <- c(68, 73, 78, 83, 88)
## confirmatory estimates of operating characteristics, which 
## are used to construct the contour plot
for (jj in 1:length(nn)){
  n <- nn[jj]
  
  ## estimate sampling distribution under H1 at n
  H1_vec <- foreach(i=1:m, .packages=c('rjags', 'coda'), .combine=rbind,
                    .options.snow=opts) %dopar% {
                      log.wd <- paste(getwd(), '/JAGS_logistic.txt', sep='')
                      
                      beta0 <- eta_plus1[1]
                      beta1 <- eta_plus1[2]
                      beta2 <- eta_plus1[3]
                      
                      ## generate data
                      set.seed(4*m + (2*jj-2)*m + i)
                      x1 <- rep(c(0,rep(1, q)), each = n)
                      x2 <- rnorm((q+1)*n)
                      eta <- beta0 + beta1*x1 + beta2*x2
                      y <- rbinom((q+1)*n, 1, expit(eta))
                      
                      ## extract MCMC settings
                      n.chains = mcmc_settings[1]
                      n.burnin = mcmc_settings[2]
                      n.draws = mcmc_settings[3]
                      n.thin = mcmc_settings[4]
                      
                      ## fit model and obtain posterior draws
                      model1.fit <- jags.model(file=log.wd,
                                               data=list(N=(q+1)*n, Y=y, X1=x1, X2=x2,
                                                         m0 = hypers[1], p0 = hypers[2],
                                                         m1 = hypers[2], p1 = hypers[4],
                                                         m2 = hypers[5], p2 = hypers[6]),
                                               n.chains = n.chains)
                      
                      update(model1.fit, n.burnin)
                      
                      model1.samples <- coda.samples(model1.fit, c("beta1"), 
                                                     n.iter=n.draws, thin=n.thin)
                      
                      beta1.post <- unlist(model1.samples[1])
                      
                      ## estimate posterior of beta1 using kernel density estimation      
                      kd <- density(beta1.post)
                      ## get the complementary of the probability H1 is true (more stable to
                      ## take the logit of this probability and multiply by -1)
                      comp.prob <- mean(pnorm(deltas[2], beta1.post, kd$bw, lower.tail = FALSE))
                      
                      res.i <- -1*logit(comp.prob)
                      c(i, res.i)
                      
                    }
  write.csv(H1_vec, paste0("H1_vec_", n,".csv"), row.names = FALSE)
  
  ## repeat the process above to estimate the sampling distribution under H0 at n1
  H0_vec <- foreach(i=1:m, .packages=c('rjags', 'coda'), .combine=rbind,
                    .options.snow=opts) %dopar% {
                      log.wd <- paste(getwd(), '/JAGS_logistic.txt', sep='')
                      
                      beta0 <- eta_plus0[1]
                      beta1 <- eta_plus0[2]
                      beta2 <- eta_plus0[3]
                      
                      set.seed(4*m + (2*jj-1)*m + i)
                      x1 <- rep(c(0,rep(1, q)), each = n)
                      x2 <- rnorm((q+1)*n)
                      eta <- beta0 + beta1*x1 + beta2*x2
                      y <- rbinom((q+1)*n, 1, expit(eta))
                      
                      n.chains = mcmc_settings[1]
                      n.burnin = mcmc_settings[2]
                      n.draws = mcmc_settings[3]
                      n.thin = mcmc_settings[4]
                      
                      model1.fit <- jags.model(file=log.wd,
                                               data=list(N=(q+1)*n, Y=y, X1=x1, X2=x2,
                                                         m0 = hypers[1], p0 = hypers[2],
                                                         m1 = hypers[2], p1 = hypers[4],
                                                         m2 = hypers[5], p2 = hypers[6]),
                                               n.chains = n.chains)
                      
                      update(model1.fit, n.burnin)
                      
                      model1.samples <- coda.samples(model1.fit, c("beta1"), 
                                                     n.iter=n.draws, thin=n.thin)
                      
                      beta1.post <- unlist(model1.samples[1])
                      
                      kd <- density(beta1.post)
                      comp.prob <- mean(pnorm(deltas[2], beta1.post, kd$bw, lower.tail = FALSE))
                      
                      res.i <- -1*logit(comp.prob)
                      c(i, res.i)
                    }
  write.csv(H0_vec, paste0("H0_vec_", n,".csv"), row.names = FALSE)
}

## we now obtain a sample size recommendation when 
## fixing the critical value at gamma = 0.6
H1_high <- H1_vec_n0[,2]; H0_high <- H0_vec_n0[,2]
H1_low <- H1_vec_n1[,2]; H0_low <- H0_vec_n1[,2]

## now we only require the linear approximations for H1
H1_slope <- (sort(H1_high) - sort(H1_low))/(n_high-n_low)
H1_int <- sort(H1_low) - H1_slope*n_low

lower <- n0; upper <- n0 + 30
## implement binary search given these bounds for the sample size
while ((upper - lower) > 1){
  n <- ceiling(0.5*(upper + lower))
  print(n)
  
  ## use linear approximations to estimate the sampling distribution under H1
  H1_vec <- H1_int + H1_slope*n
  
  ## we compare the order statistic of the sampling distribution under H1
  ## to logit(0.6) (i.e., a fixed cutoff)
  if (logit(0.6) <= sort(H1_vec)[mid2]){
    upper <- n
    H1_final <- H1_vec; H0_final <- H0_vec
  } else{
    lower <- n
  }
}

## the final sample size that we explore is n = 97
nn <- upper
## confirmatory estimates of operating characteristics
for (jj in 1:length(nn)){
  n <- nn[jj]
  ## estimate sampling distribution under H1 at n
  H1_vec <- foreach(i=1:m, .packages=c('rjags', 'coda'), .combine=rbind,
                    .options.snow=opts) %dopar% {
                      log.wd <- paste(getwd(), '/JAGS_logistic.txt', sep='')
                      
                      beta0 <- eta_plus1[1]
                      beta1 <- eta_plus1[2]
                      beta2 <- eta_plus1[3]
                      
                      ## generate data
                      set.seed(14*m + (2*jj-2)*m + i)
                      x1 <- rep(c(0,rep(1, q)), each = n)
                      x2 <- rnorm((q+1)*n)
                      eta <- beta0 + beta1*x1 + beta2*x2
                      y <- rbinom((q+1)*n, 1, expit(eta))
                      
                      ## extract MCMC settings
                      n.chains = mcmc_settings[1]
                      n.burnin = mcmc_settings[2]
                      n.draws = mcmc_settings[3]
                      n.thin = mcmc_settings[4]
                      
                      ## fit model and obtain posterior draws
                      model1.fit <- jags.model(file=log.wd,
                                               data=list(N=(q+1)*n, Y=y, X1=x1, X2=x2,
                                                         m0 = hypers[1], p0 = hypers[2],
                                                         m1 = hypers[2], p1 = hypers[4],
                                                         m2 = hypers[5], p2 = hypers[6]),
                                               n.chains = n.chains)
                      
                      update(model1.fit, n.burnin)
                      
                      model1.samples <- coda.samples(model1.fit, c("beta1"), 
                                                     n.iter=n.draws, thin=n.thin)
                      
                      beta1.post <- unlist(model1.samples[1])
                      
                      ## estimate posterior of beta1 using kernel density estimation      
                      kd <- density(beta1.post)
                      ## get the complementary of the probability H1 is true (more stable to
                      ## take the logit of this probability and multiply by -1)
                      comp.prob <- mean(pnorm(deltas[2], beta1.post, kd$bw, lower.tail = FALSE))
                      
                      res.i <- -1*logit(comp.prob)
                      c(i, res.i)
                      
                    }
  write.csv(H1_vec, paste0("H1_vec_", n,".csv"), row.names = FALSE)
  
  ## repeat the process above to estimate the sampling distribution under H0 at n1
  H0_vec <- foreach(i=1:m, .packages=c('rjags', 'coda'), .combine=rbind,
                    .options.snow=opts) %dopar% {
                      log.wd <- paste(getwd(), '/JAGS_logistic.txt', sep='')
                      
                      beta0 <- eta_plus0[1]
                      beta1 <- eta_plus0[2]
                      beta2 <- eta_plus0[3]
                      
                      set.seed(15*m + (2*jj-1)*m + i)
                      x1 <- rep(c(0,rep(1, q)), each = n)
                      x2 <- rnorm((q+1)*n)
                      eta <- beta0 + beta1*x1 + beta2*x2
                      y <- rbinom((q+1)*n, 1, expit(eta))
                      
                      n.chains = mcmc_settings[1]
                      n.burnin = mcmc_settings[2]
                      n.draws = mcmc_settings[3]
                      n.thin = mcmc_settings[4]
                      
                      model1.fit <- jags.model(file=log.wd,
                                               data=list(N=(q+1)*n, Y=y, X1=x1, X2=x2,
                                                         m0 = hypers[1], p0 = hypers[2],
                                                         m1 = hypers[2], p1 = hypers[4],
                                                         m2 = hypers[5], p2 = hypers[6]),
                                               n.chains = n.chains)
                      
                      update(model1.fit, n.burnin)
                      
                      model1.samples <- coda.samples(model1.fit, c("beta1"), 
                                                     n.iter=n.draws, thin=n.thin)
                      
                      beta1.post <- unlist(model1.samples[1])
                      
                      kd <- density(beta1.post)
                      comp.prob <- mean(pnorm(deltas[2], beta1.post, kd$bw, lower.tail = FALSE))
                      
                      res.i <- -1*logit(comp.prob)
                      c(i, res.i)
                    }
  write.csv(H0_vec, paste0("H0_vec_", n,".csv"), row.names = FALSE)
}

## we now construct the contour plots for Appendix C.1
## create matrices for left contour plot in Figure 1 (based on single
## sample size calculation)
first_rep <- results[[1]]
n_low <- first_rep[6]; n_high <- first_rep[7]

## read in the posterior probabilities corresponding to n0 and n1
H0_slope <- results[[2]]
H0_int <- results[[3]]
H1_slope <- results[[4]]
H1_int <- results[[5]]

## approximate the sampling distributions of posterior probabilities
## on the logit scale using these approximations
for (i in seq(63, 97, 1)){
  assign(paste0("H1_vec_", i), H1_int + H1_slope*i)
  assign(paste0("H0_vec_", i), H0_int + H0_slope*i)
}

## create a vector of gamma values on the logit scale to compute power
## and type I error rate estimates
opt_gamma <- as.numeric(first_rep[2])
opt_gamma <- log(opt_gamma) - log(1 - opt_gamma)

gammas <- seq(log(0.5) - log(0.5), log(0.7) - log(0.3), length.out = 50)
gammas <- sort(c(opt_gamma, gammas))

x <- seq(63, 97, 1)
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

## create matrices for right contour plot in Figure 1 (based on
## simulating data and repeatedly estimating the sampling distribution)
z_full2 <- matrix(0, nrow = 8, ncol = 50)
w_full2 <- matrix(0, nrow = 8, ncol = 50)
## convert the posterior probabilities for each approximated sampling
## distribution to the logit scale (with error checking to ensure
## to logits are finite)
for (i in c(seq(63, 93,5), 97)){
  assign(paste0("H1_vec_", i), 
         as.numeric(unlist(read.csv(paste0("H1_vec_", i,".csv"))[,2])))
  assign(paste0("H0_vec_", i), 
         as.numeric(unlist(read.csv(paste0("H0_vec_", i,".csv"))[,2])))
}

## this process mirrors what was done to create the z and w matrices in 
## the previous two plots but with the estimates obtained by simulating data
gammas <- seq(log(0.5) - log(0.5), log(0.7) - log(0.3), length.out = 50)

x <- c(seq(63, 93,5), 97)
y <- 1/(1 + exp(-gammas))

z_mat <- NULL
for (i in 1:length(x)){
  z_mat <- rbind(z_mat, get(paste0("H1_vec_", x[i])))
}
z <- NULL
for (j in 1:length(y)){
  z <- cbind(z, rowMeans(z_mat > gammas[j]))
}

z_full2 <- z

w_mat <- NULL
for (i in 1:length(x)){
  w_mat <- rbind(w_mat, get(paste0("H0_vec_", x[i])))
}
w <- NULL
for (j in 1:length(y)){
  w <- cbind(w, rowMeans(w_mat > gammas[j]))
}

w_full2 <- w

## write output to a .csv file  
write.csv(z_full2, "z_full_mat2.csv", row.names = FALSE)
write.csv(w_full2, "w_full_mat2.csv", row.names = FALSE)

## create the three contour plots and output as .pdf file for the article
pdf(file = "FigureC1.pdf",   # The directory you want to save the file in
    width = 6, 
    height = 6) 

par(mfrow=c(2,2), mar = c(3.75, 3.75, 2, 0.35) + 0.1, mgp=c(2.25,1,0))

## read in matrices for left plot
z <- matrix(unlist(read.csv("z_mat1.csv")), nrow = 35, ncol = 51)
w <- matrix(unlist(read.csv("w_mat1.csv")), nrow =35, ncol = 51)
gammas <- seq(log(0.5) - log(0.5), log(0.7) - log(0.3), length.out = 50)
gammas <- sort(c(opt_gamma, gammas))
y <- 1/(1 + exp(-gammas))
x <- seq(63, 97, 1)

contour(x, y, w, levels = c(seq(0.24, 0.38, 0.02), seq(0.42, 0.46, 0.02)), 
        xlab = expression(italic("n")['B']), xlim = c(63,97),
        ylab = expression(gamma),  main = "Type I Error Rate", labcex = 0.8, method = "edge", axes = FALSE, cex.lab = 1.25)
contour(x, y, w, levels = c(0.4), col = "firebrick", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, z, levels = c(0.75), col = "seagreen", add = TRUE, labcex = 0.8, labels = "", method = "edge")
points(x = first_rep[1], y = 1/(1 + exp(-opt_gamma)), pch = 19, col = adjustcolor("grey50", 0.75))
axis(side = 1, at = seq(65, 95, 10), cex.axis = 1.15)
axis(side = 2, at = seq(0.5, 0.7, 0.05), cex.axis = 1.15)
box()

gammas <- seq(log(0.5) - log(0.5), log(0.7) - log(0.3), length.out = 50)
y <- 1/(1 + exp(-gammas))
x <- c(seq(63, 93,5), 97)

z_full2 <- matrix(unlist(read.csv("z_full_mat2.csv")), nrow = 8, ncol = 50)
w_full2 <- matrix(unlist(read.csv("w_full_mat2.csv")), nrow = 8, ncol = 50)

contour(x, y, w_full2, levels = c(seq(0.24, 0.38, 0.02), seq(0.42, 0.46, 0.02)), 
        xlab = expression(italic("n")['B']), xlim = c(63,97),
        ylab = expression(gamma),  main = "Type I Error Rate", labcex = 0.8, method = "edge", axes = FALSE, cex.lab = 1.25)
contour(x, y, w_full2, levels = c(0.4), col = "firebrick", add = TRUE, labcex = 0.8, method = "edge")
contour(x, y, z_full2, levels = c(0.75), col = "seagreen", add = TRUE, labcex = 0.8, labels = "", method = "edge")
axis(side = 1, at = seq(65, 95, 10), cex.axis = 1.15)
axis(side = 2, at = seq(0.5, 0.7, 0.05), cex.axis = 1.15)
box()

gammas <- seq(log(0.5) - log(0.5), log(0.7) - log(0.3), length.out = 50)
gammas <- sort(c(opt_gamma, gammas))
y <- 1/(1 + exp(-gammas))
x <- seq(63, 97, 1)

contour(x, y, z, levels = c(seq(0.61, 0.73, 0.02), seq(0.77, 0.83, 0.02)), 
        xlab = expression(italic("n")['B']), xlim = c(63,97),
        ylab = expression(gamma),  main = "Power", labcex = 0.8, method = "edge", axes = FALSE, cex.lab = 1.25)
contour(x, y, w, levels = c(0.4), col = "firebrick", add = TRUE, labcex = 0.8,labels = "", method = "edge") 
contour(x, y, z, levels = c(0.75), col = "seagreen", add = TRUE, labcex = 0.8, method = "edge")
axis(side = 1, at = seq(65, 95, 10), cex.axis = 1.15)
axis(side = 2, at = seq(0.5, 0.7, 0.05), cex.axis = 1.15)
box()
points(x = first_rep[1], y = 1/(1 + exp(-opt_gamma)), pch = 19, col = adjustcolor("grey50", 0.75))

gammas <- seq(log(0.5) - log(0.5), log(0.7) - log(0.3), length.out = 50)
y <- 1/(1 + exp(-gammas))
x <- c(seq(63, 93,5), 97)

contour(x, y, z_full2, levels = c(seq(0.61, 0.73, 0.02), seq(0.77, 0.83, 0.02)), 
        xlab = expression(italic("n")['B']), xlim = c(63,97),
        ylab = expression(gamma),  main = "Power", labcex = 0.8, method = "edge", axes = FALSE, cex.lab = 1.25)
contour(x, y, w_full2, levels = c(0.4), col = "firebrick", add = TRUE, labcex = 0.8,labels = "", method = "edge") 
contour(x, y, z_full2, levels = c(0.75), col = "seagreen", add = TRUE, labcex = 0.8, method = "edge")
axis(side = 1, at = seq(65, 95, 10), cex.axis = 1.15)
axis(side = 2, at = seq(0.5, 0.7, 0.05), cex.axis = 1.15)
box()

par(mfrow=c(1,1))
dev.off()
