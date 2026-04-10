model{
  ### PRIORS
  ## Hyperpriors
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
  mu.mu.alpha2 ~ dnorm(0,0.1)
  mu.tau.alpha2 ~ dgamma(1,1)
  mu.sig.alpha2 <- 1 / sqrt(mu.tau.alpha2)

  ## Order-Specific Priors
  for (q in 1:nOrder){
    mu.beta0[q] ~ dnorm(mu.mu.beta0, mu.tau.beta0)
    mu.beta1[q] ~ dnorm(mu.mu.beta1, mu.tau.beta1)
    mu.alpha0[q] ~ dnorm(mu.mu.alpha0, mu.tau.alpha0)
    mu.alpha2[q] ~ dnorm(mu.mu.alpha2, mu.tau.alpha2)
  }
  
  tau.beta0 ~ dgamma(1,1)
  sd.beta0 <- 1 / sqrt(tau.beta0)
  tau.beta1 ~ dgamma(1,1)
  sd.beta1 <- 1 / sqrt(tau.beta1)
  
  tau.alpha0 ~ dgamma(1,1)
  sd.alpha0 <- 1 / sqrt(tau.alpha0)
  alpha1 ~ dnorm(0, 0.1)
  tau.alpha2 ~ dgamma(1,1)
  sd.alpha2 <- 1 / sqrt(tau.alpha2)
  
  ## Species-Specific Priors
  # species-specific coefficients
  for (i in 1:nTaxa) {
    beta0[i] ~ dnorm(mu.beta0[order[i]], tau.beta0)
    alpha0[i] ~ dnorm(mu.alpha0[order[i]], tau.alpha0)
    alpha2[i] ~ dnorm(mu.alpha2[order[i]], tau.alpha2)
    for (j in 1:nSite){
      beta1[i,j] ~ dnorm(mu.beta1[order[i]], tau.beta1)
    }
  }
  
  ### MODEL
  # State Process
  for (i in 1:nTaxa) {
    for (j in 1:nSite) {
      for (t in 1:nMonth){
        N[i,j,t] ~ dpois(lambda[i,j,t])
        log(lambda[i,j,t]) <- beta0[i] + beta1[i,j]
  # Observation Process
        y[i,j,t] ~ dbin(p[i,j,t], N[i,j,t])
        logit(p[i,j,t]) <- alpha0[i] + alpha1*J[j,t] + alpha2[i]*wind[j,t]
      }
    }
  }
  for (j in 1:nSite){
    for (t in 1:nMonth){
      Nsite[j,t] <- sum(N[,j,t])
    }
  }
}