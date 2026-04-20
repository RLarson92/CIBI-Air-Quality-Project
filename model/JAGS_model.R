model{
  ### PRIORS
  ## Hyperpriors
  mu.omega ~ dnorm(0,0.1)
  tau.omega ~ dgamma(1,1)
  sig.omega <- 1 / sqrt(tau.omega)
  # for abundance
  mu.mu.beta0 ~ dnorm(0, 0.1)
  mu.tau.beta0 ~ dgamma(1, 1)
  mu.sig.beta0 <- 1 / sqrt(mu.tau.beta0)
  mu.mu.beta1 ~ dnorm(0, 0.1)
  mu.tau.beta1 ~ dgamma(1, 1)
  mu.sig.beta1 <- 1 / sqrt(mu.tau.beta1)

  # for captures
  mu.mu.alpha0 ~ dnorm(0, 0.1)
  mu.tau.alpha0 ~ dgamma(1, 1)
  mu.sig.alpha0 <- 1 / sqrt(mu.tau.alpha0)
  mu.mu.alpha1 ~ dnorm(0, 0.1)
  mu.tau.alpha1 ~ dgamma(1,1)
  mu.sig.alpha1 <- 1 / sqrt(mu.tau.alpha1)

  ## Order-Specific Priors
  for (q in 1:nOrder){
    mu.beta0[q] ~ dnorm(mu.mu.beta0, mu.tau.beta0)
    mu.beta1[q] ~ dnorm(mu.mu.beta1, mu.tau.beta1)
    mu.alpha0[q] ~ dnorm(mu.mu.alpha0, mu.tau.alpha0)
    mu.alpha1[q] ~ dnorm(mu.mu.alpha1, mu.tau.alpha1)
  }
  
  tau.beta0 ~ dgamma(1,1)
  sd.beta0 <- 1 / sqrt(tau.beta0)
  tau.beta1 ~ dgamma(1,1)
  sd.beta1 <- 1 / sqrt(tau.beta1)
  
  tau.alpha0 ~ dgamma(1,1)
  sd.alpha0 <- 1 / sqrt(tau.alpha0)
  tau.alpha1 ~ dgamma(1,1)
  sd.alpha1 <- 1 / sqrt(tau.alpha1)
  
  ## Species-Specific Priors
  # species-specific coefficients
  for (i in 1:nTaxa) {
    alpha0[i] ~ dnorm(mu.alpha0[order[i]], tau.alpha0)
    alpha1[i] ~ dnorm(mu.alpha1[order[i]], tau.alpha1)
    for (j in 1:nSite){
      beta0[i,j] ~ dnorm(mu.beta0[order[i]], tau.beta0)
      beta1[i,j] ~ dnorm(mu.beta1[order[i]], tau.beta1)
      omega[i,j] ~ dnorm(mu.omega, tau.omega)
    }
  }
  
  ### MODEL
  # First sampling month
  # State Process
  for (i in 1:nTaxa) {
    for (j in 1:nSite) {
      w[i,j] ~ dbern(omega[i,j])
      for (t in 1:nMonth){
        N[i,j,t] ~ dpois(lambda[i,j,t]*w[i,j])
        log(lambda[i,j,t]) <- beta0[i,j] + beta1[i,j]*PM[j,t]
  # Observation Process
        y[i,j,t] ~ dbin(p[i,j,t]*N[i,j,t], J[j,t])
        logit(p[i,j,t]) <- alpha0[i] + alpha1[i]*wind[j,t]
      }
    }
  }
  for (j in 1:nSite){
    for (t in 1:nMonth){
      Nsite[j,t] <- sum(N[,j,t])
    }
  }
}