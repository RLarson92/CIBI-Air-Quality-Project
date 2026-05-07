inits_Indexed <- function(chain){
  gen_list <- function(chain = chain){
    list(
      N = rep(1150, nrow(tmp)),
      b_abundance = rnorm(data_list$ncov_abundance),
      tau_abundance = rgamma(data_list$ncov_abundance,1,1),
      b_taxa = matrix(
        rnorm(data_list$nTaxa * data_list$ncov_abundance),
        ncol = data_list$nTaxa,
        nrow = data_list$ncov_abundance),
      tau_taxa = matrix(
        rgamma(data_list$nTaxa * data_list$ncov_abundance,1,1),
        ncol = data_list$nTaxa,
        nrow = data_list$ncov_abundance),
      tau_shape = runif(data_list$ncov_abundance, 0.5, 2),
      tau_rate = runif(data_list$ncov_abundance, 0.5, 2),
      b_taxa_site = array(
        rnorm(data_list$ncov_abundance * data_list$nTaxa * data_list$nSite),
        dim = c(data_list$ncov_abundance, data_list$nTaxa, data_list$nSite)
      ),
      c_shape_lambda = runif(1),
      c_rate_lambda = runif(1),
      c_tau_lambda = rgamma(data_list$nSite,1,1),
      ssc_lambda = rnorm(data_list$nMonth_Params),
      c_shape_rho = runif(1),
      c_rate_rho = runif(1),
      c_tau_rho = rgamma(data_list$nSite,1,1),
      ssc_rho = rnorm(data_list$nMonth_Params),
      .RNG.name = switch(chain,
                         "1" = "base::Wichmann-Hill",
                         "2" = "base::Marsaglia-Multicarry",
                         "3" = "base::Super-Duper",
                         "4" = "base::Mersenne-Twister",
                         "5" = "base::Wichmann-Hill",
                         "6" = "base::Marsaglia-Multicarry",
                         "7" = "base::Super-Duper",
                         "8" = "base::Mersenne-Twister"),
      .RNG.seed = sample(1:1e+06, 1)
    )
  }
  return(
    switch(
      chain,
      "1" = gen_list(chain),
      "2" = gen_list(chain),
      "3" = gen_list(chain),
      "4" = gen_list(chain),
      "5" = gen_list(chain),
      "6" = gen_list(chain),
      "7" = gen_list(chain),
      "8" = gen_list(chain)
    )
  )
}