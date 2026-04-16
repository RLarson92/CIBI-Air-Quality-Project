model{
  ### PRIORS
  ## Hyperpriors
  omega ~ dunif(0,1)
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

  ## Order-Specific Priors
  for (q in 1:nOrder){
    mu.beta0[q] ~ dnorm(mu.mu.beta0, mu.tau.beta0)
    mu.beta1[q] ~ dnorm(mu.mu.beta1, mu.tau.beta1)
    mu.alpha0[q] ~ dnorm(mu.mu.alpha0, mu.tau.alpha0)
  }
  
  tau.beta0 ~ dgamma(1,1)
  sd.beta0 <- 1 / sqrt(tau.beta0)
  tau.beta1 ~ dgamma(1,1)
  sd.beta1 <- 1 / sqrt(tau.beta1)
  mu.phi ~ dnorm(0, 0.1)
  tau.phi ~ dgamma(1,1)
  sd.phi <- 1 / sqrt(tau.phi)
  
  tau.alpha0 ~ dgamma(1,1)
  sd.alpha0 <- 1 / sqrt(tau.alpha0)
  
  ## Species-Specific Priors
  # species-specific coefficients
  for (i in 1:nTaxa) {
    alpha0[i] ~ dnorm(mu.alpha0[order[i]], tau.alpha0)
    for (j in 1:nSite){
      beta0[i,j] ~ dnorm(mu.beta0[order[i]], tau.beta0)
      beta1[i,j] ~ dnorm(mu.beta1[order[i]], tau.beta1)
      phi[i,j] ~ dnorm(mu.phi, tau.phi)
    }
  }
  
  ### MODEL
  # First sampling month
  # State Process
  for (i in 1:nTaxa) {
    w[i] ~ dbern(omega)
    for (j in 1:nSite) {
      N[i,j,1] ~ dpois(lambda[i,j,1])
      log(lambda[i,j,1]) <- beta0[i,j]
  # Observation Process
      y[i,j,1] ~ dbin(p[i,j,1], N[i,j,1])
      logit(p[i,j,1]) <- alpha0[i]
  # Subsequent Sampling months
  # State Process
      for (t in 2:nMonth){
        N[i,j,t] ~ dpois(lambda[i,j,t])
        log(lambda[i,j,t]) <- beta0[i,j] + phi[i,j]
  # Observation Process
        y[i,j,t] ~ dbin(p[i,j,t], N[i,j,t])
        logit(p[i,j,t]) <- alpha0[i]
      }
    }
  }
  for (j in 1:nSite){
    for (t in 1:nMonth){
      Nsite[j,t] <- sum(N[,j,t])
    }
  }
}