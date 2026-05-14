model{
  ### PRIORS
  ## Hyperpriors
  mu.omega ~ dnorm(0, 0.1)
  tau.omega ~ dgamma(1,1)
  sig.omega <- 1 / sqrt(tau.omega)
  # for abundance
  mu.mu.beta0 ~ dnorm(0, 0.1)
  mu.tau.beta0 ~ dgamma(1, 1)
  mu.sig.beta0 <- 1 / sqrt(mu.tau.beta0)
  mu.mu.beta1 ~ dnorm(0, 0.1)
  mu.tau.beta1 ~ dgamma(1, 1)
  mu.sig.beta1 <- 1 / sqrt(mu.tau.beta1)
  mu.mu.beta2 ~ dnorm(0, 0.1)
  mu.tau.beta2 ~ dgamma(1, 1)
  mu.sig.beta2 <- 1 / sqrt(mu.tau.beta2)
  mu.mu.beta3 ~ dnorm(0, 0.1)
  mu.tau.beta3 ~ dgamma(1, 1)
  mu.sig.beta3 <- 1 / sqrt(mu.tau.beta3)
  mu.phi ~ dnorm(0, 0.1)
  tau.phi ~ dgamma(1,1)
  sig.phi <- 1 / sqrt(tau.phi)

  # for captures
  mu.mu.alpha0 ~ dnorm(0, 0.1)
  mu.tau.alpha0 ~ dgamma(1, 1)
  mu.sig.alpha0 <- 1 / sqrt(mu.tau.alpha0)
  alpha1 ~ dnorm(0, 0.1)
  mu.mu.alpha2 ~ dnorm(0, 0.1)
  mu.tau.alpha2 ~ dgamma(1,1)
  mu.sig.alpha2 <- 1 / sqrt(mu.tau.alpha2)
  mu.mu.alpha3 ~ dnorm(0, 0.1)
  mu.tau.alpha3 ~ dgamma(1,1)
  mu.sig.alpha3 <- 1 / sqrt(mu.tau.alpha3)

  ## Order-Specific Priors
  for (q in 1:nOrder){
    mu.beta0[q] ~ dnorm(mu.mu.beta0, mu.tau.beta0)
    mu.beta1[q] ~ dnorm(mu.mu.beta1, mu.tau.beta1)
    mu.beta2[q] ~ dnorm(mu.mu.beta2, mu.tau.beta2)
    mu.beta3[q] ~ dnorm(mu.mu.beta3, mu.tau.beta3)
    mu.alpha0[q] ~ dnorm(mu.mu.alpha0, mu.tau.alpha0)
    mu.alpha2[q] ~ dnorm(mu.mu.alpha2, mu.tau.alpha2)
    mu.alpha3[q] ~ dnorm(mu.mu.alpha3, mu.tau.alpha3)
  }
  
  tau.beta0 ~ dgamma(1,1)
  sd.beta0 <- 1 / sqrt(tau.beta0)
  tau.beta1 ~ dgamma(1,1)
  sd.beta1 <- 1 / sqrt(tau.beta1)
  tau.beta2 ~ dgamma(1,1)
  sd.beta2 <- 1 / sqrt(tau.beta2)
  tau.beta3 ~ dgamma(1,1)
  sd.beta3 <- 1 / sqrt(tau.beta3)
  
  tau.alpha0 ~ dgamma(1,1)
  sd.alpha0 <- 1 / sqrt(tau.alpha0)
  tau.alpha2 ~ dgamma(1,1)
  sd.alpha2 <- 1 / sqrt(tau.alpha2)
  tau.alpha3 ~ dgamma(1,1)
  sd.alpha3 <- 1 / sqrt(tau.alpha3)
  
  ## Species-Specific Priors
  for (i in 1:nTaxa) {
    alpha0[i] ~ dnorm(mu.alpha0[order[i]], tau.alpha0)
    alpha2[i] ~ dnorm(mu.alpha2[order[i]], tau.alpha2)
    alpha3[i] ~ dnorm(mu.alpha3[order[i]], tau.alpha3)
    phi[i] ~ dnorm(mu.phi, tau.phi)
    for (j in 1:nSite){
      beta0[i,j] ~ dnorm(mu.beta0[order[i]], tau.beta0)
      beta1[i,j] ~ dnorm(mu.beta1[order[i]], tau.beta1)
      beta2[i,j] ~ dnorm(mu.beta2[order[i]], tau.beta2)
      beta3[i,j] ~ dnorm(mu.beta3[order[i]], tau.beta3)
      omega[i,j] ~ dnorm(mu.omega, tau.omega)
    }
  }
  
  ### MODEL
  # First sampling month
  # State Process
  for (i in 1:nTaxa) {
    for (j in 1:nSite) {
      w[i,j] ~ dbern(1-omega[i,j])
      N[i,j,1] ~ dpois(lambda[i,j,1]*w[i,j])
      log(lambda[i,j,1]) <- beta0[i,j] + beta1[i,j]*PM[j,1] + 
        beta2[i,j]*smoke[j,1] + beta3[i,j]*PC1[j,1]
      # Observation Process
      y[i,j,1] ~ dbinom(p[i,j,1], N[i,j,1])
      logit(p[i,j,1]) <- alpha0[i] + alpha1*J[j,1] + alpha2[i]*PC1[j,1] + 
        alpha3[i]*PC2[j,1]
  # Subsequent Sampling Months
  # State Process
      for (t in 2:nMonth){
        N[i,j,t] ~ dpois(lambda[i,j,t]*w[i,j])
        log(lambda[i,j,t]) <- beta0[i,j] + beta1[i,j]*PM[j,t] + 
          beta2[i,j]*smoke[j,t] + beta3[i,j]*PC1[j,t] + phi[i]*log(N[i,j,t-1]+1)
        # Observation Process
        y[i,j,t] ~ dbinom(p[i,j,t], N[i,j,t])
        logit(p[i,j,t]) <- alpha0[i] + alpha1*J[j,t] + alpha2[i]*PC1[j,t] + 
          alpha3[i]*PC2[j,t]
      }
    }
  }
  # Derived Parameters
  for(i in 1:nTaxa){
    Ntotal[i]<-sum(N[i,,]) # total abundance of each species
  }
  for(j in 1:nSite){
    for(t in 1:nMonth){
      Nsite[j,t]<-sum(N[,j,t]) # abundance of species at each site each month
    }
  }
}