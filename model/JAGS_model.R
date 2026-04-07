model{
  ### PRIORS
  ## Hyperpriors (community average distributions to pull species-specific responses from)
  # for abundance
  mu.beta0 ~ dnorm(0, 0.1)
  tau.beta0 ~ dgamma(1, 1)
  sig.beta0 <- 1 / sqrt(tau.beta0)
  mu.beta1 ~ dnorm(0, 0.1)
  tau.beta1 ~ dgamma(1, 1)
  sig.beta1 <- 1 / sqrt(tau.beta1)
  # for captures
  mu.alpha0 ~ dnorm(0, 0.1)
  tau.alpha0 ~ dgamma(1, 1)
  sig.alpha0 <- 1 / sqrt(tau.alpha0)
  alpha1 ~ dnorm(0, 0.1)                   # technically not hyperpriors but here for convenience
  alpha2 ~ dnorm(0, 0.1)
  alpha3 ~ dnorm(0, 0.1)
  # for survival
  mu.phi0 ~ dnorm(0, 0.1)
  tau.phi0 ~ dgamma(1, 1)
  sig.phi0 <- 1 / sqrt(tau.phi0)
  mu.phi1 ~ dnorm(0, 0.1)
  tau.phi1 ~ dgamma(1, 1)
  sig.phi1 <- 1 / sqrt(tau.phi1)
  mu.phi2 ~ dnorm(0, 0.1)
  tau.phi2 ~ dgamma(1, 1)
  sig.phi2 <- 1 / sqrt(tau.phi2)
  # for recruitment
  mu.gamma0 ~ dnorm(0, 0.1)
  tau.gamma0 ~ dgamma(1, 1)
  sig.gamma0 <- 1 / sqrt(tau.gamma0)
  mu.gamma1 ~ dnorm(0, 0.1)
  tau.gamma1 ~ dgamma(1, 1)
  sig.gamma1 <- 1 / sqrt(tau.gamma1)
  # because season is dummy-coded, we have a different coefficient for each season
  # although season isn't a hyper prior, it's here for convenience
  for (tmpk in 1:2){
    tmp_gamma[tmpk] ~ dnorm(0,0.1)
  }
  gamma2[1] <- 0
  for (z in 2:nseason){
    gamma2[z] <- tmp_gamma[z-1]
  }
  ## Species-Specific Priors
  # species-specific coefficients
  for (i in 1:nspec) {
    # for abundance
    beta0[i] ~ dnorm(mu.beta0, tau.beta0)
    beta1[i] ~ dnorm(mu.beta1, tau.beta1)
    # for captures
    alpha0[i] ~ dnorm(mu.alpha0, tau.alpha0)
    # for survival
    phi0[i] ~ dnorm(mu.phi0, tau.phi0)
    phi1[i] ~ dnorm(mu.phi1, tau.phi1)
    phi2[i] ~ dnorm(mu.phi2, tau.phi2)
    # for recruitment
    gamma0[i] ~ dnorm(mu.gamma0, tau.gamma0)
    gamma1[i] ~ dnorm(mu.gamma1, tau.gamma1)
  }
  
  ### MODEL
  # State Process (initial abundance)
  for (i in 1:nTaxa) {
    for (j in 1:nSite) {
      N[i,j,1] ~ dpois(lambda[i,j,1])
      log(lambda[i,j,1]) <- beta0[i] + beta1[i]*temp[j,1]
      # Observation Process
      y[i,j,1] ~ dbin(p[i,j,1], N[i,j,1])
      logit(p[i,j,1]) <- alpha0[i] + alpha1[i]*wind[j,1] + alpha2*temp[j,1] + alpha3*effort[j,1]
    }
      
      # State Process (transition model)
      for (t in 2:nMonth) {
        N[i,j,t] <- S[i,j,t] + R[i,j,t]  # 'S' for survivors, 'R' for recruits
        S[i,j,t] ~ dbin(phi[i,j,t], N[i,j,t-1])
        logit(phi[i,j,t]) <- phi0[i] + phi1[i]*PM2.5[j] + phi2[i]*temp[j]
        R[i,j,t] ~ dpois(gamma[i,j,t])
        log(gamma[i,j,t]) <- gamma0[i] + gamma1[i]*contag[j] + gamma2[season[t]]
        # Observation Process
        logit(p[i,j,t]) <- alpha0[i] + alpha1[i]*wind[j,t] + alpha2*temp[j,t] + alpha3*effort[j,t]
        y[i,j,t] ~ dbin(p[i,j,t], N[i,j,t])
      }
    }
  }