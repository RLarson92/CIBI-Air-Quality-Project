model{
  # Priors for Abundance
  for(l in 1:ncov_abundance){
    b_abundance[l] ~ dt(0, 2.5, 1)
    tau_abundance[l] ~ dgamma(1,1)
    tau_shape[l] ~ dunif(0.01,10)  # hyperpriors for species random effects
    tau_rate[l] ~ dunif (0.01, 10)
    for(j in 1:nTaxa){
      b_taxa[l,j] ~ dnorm(b_abundance[l], tau_abundance[l])
      tau_taxa[l,j] ~ dgamma(tau_shape[l], tau_rate[l])
      for(k in 1:nSite){
        b_taxa_site[l,j,k] ~ dnorm(b_taxa[l,j], tau_taxa[l,j])
      }
    }
  }
  # Setting up monthly variation now that we have a structure set up for a taxon
  # at a site
  c_shape_lambda ~ dunif(0.001, 10)
  c_rate_lambda ~ dunif(0.001, 10)
  for(k in 1:nSite){
    c_tau_lambda[k] ~ dgamma(c_shape_lambda, c_rate_lambda)
  }
  # Reference the correct month_site_tau for the month random effect
  for(m in 1:nMonth_Params){
    # taxa, month, site seasonal variation
    ssc_lambda[m] ~ dnorm(
      b_taxa_site[1,combo_taxa_idx[m],combo_site_idx[m]],
      c_tau_lambda[combo_site_idx[m]]
    )
  }
  # Priors for Capture Rates
  for(l in 1:ncov_capture){
    c_cap[l] ~ dt(0,2.5,1)
    tau_cap[l] ~ dgamma(1,1)
    tau_shape_rho[l] ~ dunif(0.01,10) # hyperpriors for species random effects
    tau_rate_rho[l] ~ dunif(0.01,10)
    for(j in 1:nTaxa){
      c_taxa_cap[l,j] ~ dnorm(c_cap[l], tau_cap[l])
      tau_taxa_cap[l,j] ~ dgamma(tau_shape_rho[l], tau_rate_rho[l])
      for(k in 1:nSite){
        c_taxa_site[l,j,k] ~ dnorm(c_taxa_cap[l,j], tau_taxa_cap[l,j])
      }
    }
  }
  # Setting up monthly variation now that we have a structure set up for a taxon
  # at a site
  c_shape_rho ~ dunif(0.001, 10)
  c_rate_rho ~ dunif(0.001, 10)
  for(k in 1:nSite){
    c_tau_rho[k] ~ dgamma(c_shape_rho, c_rate_rho)
  }
  # Reference the correct month_site_tau for the month random effect
  for(m in 1:nMonth_Params){
    # taxa, month, site seasonal variation
    ssc_rho[m] ~ dnorm(
      c_taxa_site[1,combo_taxa_idx[m],combo_site_idx[m]],
      c_tau_rho[combo_site_idx[m]]
    )
  }
  # Abundance (latent state)
  for(s in 1:nSamples){
    log(lambda[s]) <- inprod(b_taxa_site[2:ncov_abundance,taxa_idx[s],site_idx[s]], 
                             lambda_covs[s,2:ncov_abundance]) + ssc_lambda[combo_idx[s]]
    N[s] ~ dpois(lambda[s])
  }
  # Captures
  for(s in 1:nSamples){
    logit(rho[s]) <- inprod(c_taxa_site[2:ncov_capture,taxa_idx[s],site_idx[s]],
                            rho_covs[s,2:ncov_capture]) + ssc_rho[combo_idx[s]]
    y[s] ~ dbinom(rho[s], N[s])
  }
}