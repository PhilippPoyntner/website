rm(list = ls())

compute_elbo <- function(y, x, beta_mu, beta_sd, nu, tau2, nr_samples = 1e4) {
  n <- length(y)
  sum_y2 <- sum(y^2)
  sum_x2 <- sum(x^2)
  sum_yx <- sum(x*y)
  
  # Takes a function and computes its expectation with respect to q(\beta)
  E_q_beta <- function(fn) {
    integrate(function(beta) {
      dnorm(beta, beta_mu, beta_sd) * fn(beta)
    }, -Inf, Inf)$value
  }
  
  # Takes a function and computes its expectation with respect to q(\sigma^2)
  E_q_sigma2 <- function(fn) {
    integrate(function(sigma) {
      MCMCpack::dinvgamma(sigma^2, (n + 1)/2, nu) * fn(sigma)
    }, 0, Inf)$value
  }
  
  # Compute expectations of log p(\sigma^2)
  E_log_p_sigma2 <- E_q_sigma2(function(sigma) log(1/sigma^2))
  
  # Compute expectations of log p(\beta | \sigma^2)
  E_log_p_beta <- (
    log(tau2 / beta_sd^2) * E_q_sigma2(function(sigma) log(sigma^2)) +
      (beta_sd^2 + tau2) / (tau2) * E_q_sigma2(function(sigma) 1/sigma^2)
  )
  
  # Compute expectations of the log variational densities q(\beta)
  E_log_q_beta <- E_q_beta(function(beta) dnorm(beta, beta_mu, beta_sd, log = TRUE))
  # E_log_q_sigma2 <- E_q_sigma2(function(x) log(MCMCpack::dinvgamma(x, (n + 1)/2, nu))) # fails
  
  # Compute expectations of the log variational densities q(\sigma^2)
  sigma2 <- MCMCpack::rinvgamma(nr_samples, (n + 1)/2, nu)
  E_log_q_sigma2 <- mean(log(MCMCpack::dinvgamma(sigma2, (n + 1)/2, nu)))
  
  # Compute the expected log likelihood
  E_log_y_b <- sum_y2 - 2*sum_yx*beta_mu + (beta_sd^2 + beta_mu^2)*sum_x2
  E_log_y_sigma2 <- E_q_sigma2(function(sigma) log(sigma^2) * 1/sigma^2)
  E_log_y <- n/4 * log(2*pi) * E_log_y_b * E_log_y_sigma2
  
  # Compute and return the ELBO
  ELBO <- E_log_y + E_log_p_beta + E_log_p_sigma2 - E_log_q_beta - E_log_q_sigma2
  ELBO
}

lmcavi <- function(y, x, tau2, nr_samples = 1e5, epsilon = 1e-4) {
  n <- length(y)
  sum_y2 <- sum(y^2)
  sum_x2 <- sum(x^2)
  sum_yx <- sum(x*y)
  
  # beta mean (independent of optimization)
  beta_mu <- sum_yx / (sum_x2 + 1/tau2)
  # sigma^2 shape parameter 
  cN <- (n+1)/2
  
  # starting values
  res <- list()
  res[['nu']] <- 5
  res[['cN']] <- cN
  res[['beta_mu']] <- beta_mu
  res[['beta_sd']] <- 1
  res[['E_sigma2']] <- NA
  res[['sigma2_sd']]<- NA
  res[['ELBO']] <- 0

  ELBO <- compute_elbo(y = y,
                       x = x,
                       beta_mu = beta_mu,
                       beta_sd = res[['beta_sd']],
                       nu = res[['nu']],
                       tau2 = tau2,
                       nr_samples = nr_samples)
  
  # while the ELBO has not converged
  j <- 1 # init counter
  while(abs(res[['ELBO']][j] - ELBO) > epsilon){
    
    nu_prev <- res[['nu']][j]
    beta_sd_prev <- res[['beta_sd']][j]
    
    # used in the update of beta_sd and nu
    E_qA <- sum_y2 - 2*sum_yx*beta_mu + (beta_sd_prev^2 + beta_mu^2)*(sum_x2 + 1/tau2)
    
    # update the variational parameters for sigma2 and beta
    nu <- 1/2 * E_qA
    beta_sd <- sqrt(((n + 1) / E_qA) / (sum_x2 + 1/tau2))
    
    E_sigma2 <- nu/(cN-1)
    sigma2_sd<- nu/((cN-1)*(cN-2)^0.5)
    
    # update results object
    res[['nu']] <- c(res[['nu']], nu)
    res[['beta_sd']] <- c(res[['beta_sd']], beta_sd)
    res[['E_sigma2']] <- E_sigma2
    res[['sigma2_sd']]<- sigma2_sd
    res[['ELBO']] <- c(res[['ELBO']], ELBO)
    
    # compute new ELBO
    ELBO <- compute_elbo(y, x, beta_mu, beta_sd, nu, tau2, nr_samples = nr_samples)
    j <- j + 1
  }
  res[['summary']]<- cbind('mean' = c('beta' = beta_mu, 'sigma2' = E_sigma2), 'sd' = c(beta_sd, sigma2_sd))
  
  res
}

lmgibbs <- function(y, x, tau2, Nstore = 10000, Nburn = 2000){
  
  n <- length(y)
  # Prior Parameters
  tau2 = tau2
  
  # Sampler Setup
  Nrep   <- Nstore+Nburn
  
  store  <- matrix(NA, ncol = 2, nrow = Nstore)
  colnames(store) <- c("beta_i", "sigma2_i")
  
  # Starting Value
  beta_i <- 1
  
  for(i in -(Nburn-1):Nstore){
    
    # draw sigma|beta
    cN      <- (n+1)/2
    CN      <- 1/2 * sum((y - x * beta_i)^2 + beta_i^2 / tau2)
    sigma2_i<- 1 / rgamma(1,shape = cN, rate = CN)
    
    # draw beta|sigma
    bN      <- 1 / (sum(x^2) + 1 / tau2) * sum(x * y)
    BN      <- sigma2_i * 1 / (sum(x^2) + 1 / tau2)
    beta_i  <- rnorm(1,bN,BN^0.5)
    
    # store
    if(i > 0){
      store[i,] <- c(beta_i, sigma2_i)
    }
  }
  res <- list()
  res[['draws']] <- store
  res[['summary']] <- cbind('mean'=colMeans(store), 'sd'=diag(var(store))^0.5)
  res
}

gen_dat <- function(n, beta, sigma) {
  x <- rnorm(n)
  y <- 0 + beta*x + rnorm(n, 0, sigma)
  data.frame(x = x, y = y)
}

set.seed(1)
dat <- gen_dat(100, 0.30, 1)

tau2 <- 10
VI_res <- lmcavi(y = dat$y, x = dat$x, tau2 = tau2)
GBS_res<- lmgibbs(dat$y, dat$x, tau2 = tau2, Nstore = 10000, Nburn = 1000)

VI_res$summary
GBS_res$summary


plot(density(GBS_res$draws[,1]),lwd=2, main = "Gibbs vs. VI")
curve(dnorm(x,VI_res$beta_mu,VI_res$beta_sd[length(VI_res$beta_sd)]),add=T,lwd=2,col="red")
leg.txt <- c("Gibbs","VI")
legend("topright",legend = leg.txt, col = c("black","red"),lty=c(1,1),lwd=3)

plot(density(GBS_res$draws[,2]),lwd=2, main = "Gibbs vs. VI")
curve(MCMCpack::dinvgamma(x,VI_res$cN,VI_res$nu[length(VI_res$nu)]),add=T,lwd=2,col="red")
leg.txt <- c("Gibbs","VI")
legend("topright",legend = leg.txt, col = c("black","red"),lty=c(1,1),lwd=3)

# set.seed(1)
# dat <- gen_dat(100, 0.30, 1)
# x100 <- microbenchmark::microbenchmark(
#   "VI" = lmcavi(y = dat$y, x = dat$x, tau2 = tau2),
#   "Gibbs" = lmgibbs(dat$y, dat$x, tau2 = tau2, Nstore = 10000, Nburn = 1000), 
#   times = 20)
# 
# set.seed(1)
# dat <- gen_dat(1000, 0.30, 1)
# x1000 <- microbenchmark::microbenchmark(
#   "VI" = lmcavi(y = dat$y, x = dat$x, tau2 = tau2),
#   "Gibbs" = lmgibbs(dat$y, dat$x, tau2 = tau2, Nstore = 10000, Nburn = 1000), 
#   times = 20)
# 
# set.seed(1)
# dat <- gen_dat(10000, 0.30, 1)
# x10000 <- microbenchmark::microbenchmark(
#   "VI" = lmcavi(y = dat$y, x = dat$x, tau2 = tau2),
#   "Gibbs" = lmgibbs(dat$y, dat$x, tau2 = tau2, Nstore = 10000, Nburn = 1000), 
#   times = 20)
# 
# list("n=100"=x100,"n=1k"=x1000,"n=10k"=x10000)
