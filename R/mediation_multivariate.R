
#' BayesMediation
#' 
#' Bayes approach for FLOW MR multivariate model.
#' 
#' @param Gamma_hat The estimated exposure and outcome effects gamma_hat
#' @param Sd_hat The standarded errors of gamma_hat
#' @param init Starting value for gibbs sampler. Either "Random" or "EM"
#' @param iter The number of iterations for gibbs sampler
#' @param warmup The length of warm-up periods
#' @param second Whether to run the second stage. Default is False.
#' @param inv When inv = False, we are estiming B in Gamma = (I + B) alpha;
#' inv = True, we are estimating B in Gamma = B Gamma + alpha.
#' @param cor The correlation matrix of noise.
#' @param Raw Whether to include the unprocessed raw outputs. Default
#' is F
#' @param total Whether to include the total effects. Default is F
#' @param indirect Whether to include the indirect effects. Default is F.
#' @param n_traits number of traits at each time stage
#' @param alpha_0 hyper-parameter
#' @param alpha_1 hyper-parameter
#'  
#' @return A list with elements
#' \item{summary_first} The summary table of first stage
#' \item{summary_second} The summary table of second stage (included if second  = T)
#' \item{total_effect_first} The total effect of each exposure on the outcome,
#' computed by first stage (included if total = T)
#' \item{total_effect_second} The total effect of each exposure on the outcome,
#' computed by second stage (included if total = T and second = T)
#' \item{indirect_effect_first} The indirect effect of each exposure on the outcome,
#' computed by first stage (included if indirect = T)
#' \item{indirect_effect_second} The indirect effect of each exposure on the outcome,
#' computed by second stage (included if indirect = T and  second = T)
#' \item{raw_first} The unprocessed output from the first stage (included if raw= T)
#' \item{raw_second} The unprocessed output from the second stage (included if raw = T and second = T)
#' 
#' @export
MultivariateMediation = function(Gamma_hat, Sd_hat, init = "Random", iter = 6000,
                                 warmup = 3000, second = F,  inv = T,
                                 cor = NULL, Raw = T, total = F, indirect = F, n_traits = NULL, 
                                 alpha_0 = 3, alpha_1 = 3) {
  registerDoParallel(cores=4)
  P = length(Gamma_hat[1,])
  K = length(Gamma_hat[, 1])
  
  #Estimate pleiotropy for exposures and outcome
  exps = list()
  for(i in 1: (K-1)) {
    exps[[i]] = zero.centered.em(Gamma_hat[(i + 1), ], 
                                 Sd_hat[(i + 1),], show_outputs = FALSE)
  }
  out = zero.centered.em(Gamma_hat[(1), ], 
                         Sd_hat[(1),], show_outputs = FALSE)
  
  sigma1 = rep(0, K)
  sigma0 = rep(0, K)
  pi = rep(0, K)
  for(i in 1:(K - 1)){
    sigma1[i + 1] = exps[[i]]$S1
    sigma0[i + 1] = exps[[i]]$S2
    pi[i + 1] = exps[[i]]$Pi
  }
  
  sigma1[1] = out$S1 
  sigma0[1] = out$S2
  pi[1] = out$Pi 
  
  #First Stage
  result1 = gibbs_wrapper_multivariate(iter, warmup, chains = 4,
                                       Gamma_hat, Sd_hat, sigma = 0.1, 
                                       sigma_1 = sigma1,
                                       sigma_0 = sigma0,
                                       p = pi,
                                       parallel = TRUE, init = init, cor_mat = cor, num_traits = n_traits, alpha_0 = alpha_0,
                                       alpha_1 = alpha_1)
  
  
  
  summarys1  = summary_gibbs(result1, c("B", 
                                        "sigma", 
                                        'sigma1', "sigma0", 
                                        "p"), K, T = warmup, inv = inv)
  
  if(max(summarys1[, 8][! is.na(summarys1[, 8])]) > 1.1){
    warning("The maximum R_hat is greater than 1.1. The gibbs sampler might fail
            to converge in the first stage")
  }
  # second stage  (if necessary)
  if(second == TRUE){
    B_hat_vec = summarys1[1:((K^2-K)/2), 1]
    B_hat = matrix(0, K, K)
    B_hat[lower.tri(B_hat, diag = FALSE)] = B_hat_vec
    B_hat = t(B_hat)
    if(inv == FALSE){
      B_hat = convert(B_hat)
    }
    
    #Note (I - B_hat) Gamma_hat = alpha + (I - B_hat) epsilon
    alpha_hat = (diag(K) - B_hat) %*% Gamma_hat
    cor_hat = matrix(nrow = K, ncol = P)
    for(c in 1:P){
      cor_hat[, c] = diag( (diag(K) - B_hat) %*% 
                             diag(Sd_hat[,c])^2 %*% t(diag(K) - B_hat) )
    }
    
    exps = list()
    for(i in 1: (K-1)) {
      exps[[i]] = zero.centered.em(alpha_hat[(i + 1), ], 
                                   sqrt( cor_hat[(i + 1),] ), show_outputs = FALSE)
    }
    out = zero.centered.em(alpha_hat[(1), ], 
                           sqrt( cor_hat[(1),] ), show_outputs = FALSE)
    
    sigma1 = rep(0, K)
    sigma0 = rep(0, K)
    pi = rep(0, K)
    for(i in 1:(K - 1)){
      sigma1[i + 1] = exps[[i]]$S1
      sigma0[i + 1] = exps[[i]]$S2
      pi[i + 1] = exps[[i]]$Pi
    }
    
    sigma1[1] = out$S1 
    sigma0[1] = out$S2
    pi[1] = out$Pi 
    
    result2 = gibbs_wrapper(iter, warmup, chains = 4,
                            Gamma_hat, Sd_hat, sigma = 0.1, 
                            sigma_1 = sigma1,
                            sigma_0 = sigma0,
                            p = pi,
                            parallel = TRUE, init = init, cor_mat = cor, num_traits = n_traits, alpha_0 = alpha_0,
                            alpha_1 = alpha_1)
    
    
    summarys2  = summary_gibbs(result2, c("B", 
                                          "sigma", 
                                          'sigma1', "sigma0", 
                                          "p"), K, T = warmup, inv = inv)
    
    if(max(summarys1[, 8][! is.na(summarys1[, 8])]) > 1.1){
      warning("The maximum R_hat is greater than 1.1. The gibbs sampler might fail
            to converge in the first stage")
    }
    
    ls = list(summary_first = summarys1, summary_second = summarys2)
    if(total){
      ls$total_effect_first = total_effect(result1, K, warmup)
      ls$total_effect_second = total_effect(result2, K, warmup)
    }
    if(indirect){
      ls$indirect_effect_first = indirect_effect(result1, K, warmup)
      ls$indirect_effect_second = indirect_effect(result2, K, warmup)
      
    }
    
    if(Raw == T){
      ls$raw_first = result1
      ls$raw_second = result2
      
    }
    return(ls)
    
  }
  ls = list(summary= summarys1)
  if(total){
    ls$total_effect = total_effect(result1, K, warmup)
  }
  if(indirect){
    ls$indirect_effect= indirect_effect(result1, K, warmup)
  }
  
  if(Raw == T){
    ls$raw = result1
  }
  return(ls)
  
  
}




#' Gibbs sampler
#' 
#' Given the observed effect and standard error, run gibbs sampler to draw 
#' samples from the posterior density
#' 
#' @param N Number of total iterations
#' @param warmup Number of warm-up iterations.
#' @param chains Number of Chains.
#' @param Gamma_hat The observed effects
#' @param Sd_hat The standard deviation of Gamma_hat
#' @param sigma The prior mean for sigma
#' @param sigma_1 The prior mean for sigma_1
#' @param sigma_0 The prior mean for sigma_0
#' @param p The prior mean for p
#' @param parallel Indicating whether to parallelly run all the chains. 
#' Default is True.
#' @param init Method of Initialization. Must be one of "EM" or "Random". Default
#' is "Random"
#' @param ratio The ratio of initial sigma1 and sigma0 when the initialization method
#' is "Random". Default is 10.
#' @param cor_mat The correlation matrix for noise. Default is null (no correlation)
#' @param Lambda The block diagonal matrix such that each block is the covariance 
#' matrix of the noise for a SNP. Default is NULL (not provided)
#' @param Lambda_inv The block diagonal matrix such that each block is the inverse of 
#' covariance matrix of the noise for a SNP. Default is NULL (not provided).
#' @param alpha_0 hyper-parameter
#' @param alpha_1 hyper-parameter
#' 
#' 
#' @return A list of which the elements are outputs from each chain. Each of these
#' outputs is a list with attributes
#' \item{B}
#' \item{p}
#' \item{sigma1}
#' \item{sigma0}
#' \item{pi}
#' \item{sigma}
#' 
#' @import foreach
#' @import doParallel
#' @import coda
#' @import stats
#' @export

gibbs_wrapper_multivariate = function(N, warmup, chains = 4,
                                      Gamma_hat, Sd_hat, sigma, sigma_1, sigma_0, p,
                                      parallel = TRUE, init = "Random", ratio = 10, cor_mat = NULL, Lambda = NULL,
                                      Lambda_inv = NULL, num_traits = NULL, alpha_0 = 3, alpha_1 = 3) {
  P = length(Gamma_hat[1,])
  K = length(Gamma_hat[, 1])
  
  #Set Hyperparameters
  alpha_B = 2
  beta_B = 0.5
  alpha_0 = rep(alpha_0, K)
  alpha_1 = rep(alpha_1, K)
  stopifnot(length(sigma_0) == K)
  stopifnot(length(sigma_1) == K)
  beta_0 = sigma_0^2 * (alpha_0[1] - 1)
  beta_1 = sigma_1^2 * (alpha_1[1] - 1)
  a = rep(2, K)
  stopifnot(length(p) == K)
  b = a/p - a
  
  A = matrix(nrow=K, ncol=P);
  Z = matrix(nrow=K, ncol = P);
  
  
  outputs = list()
  
  #If correlation is non-zero, compute the covariance and inverse covariance
  #matrices, if not provided.
  if(! is.null(cor_mat) | ! is.null(num_traits)){
    if (is.null(cor_mat)){
      cor_mat = diag(K)
    }
    
    if(is.null(Lambda_inv)){
      require(magic)
      
      new_Lambda = diag( Sd_hat[,1] )%*% cor_mat %*% diag( Sd_hat[,1] )
      
      Lambda = new_Lambda
      Lambda_inv = solve(new_Lambda )
      
      for(i in 2:length(Sd_hat[1, ])){
        
        new_Lambda = diag( Sd_hat[,i] )%*% cor_mat %*% diag( Sd_hat[,i] )
        
        Lambda = adiag(Lambda, new_Lambda)
        Lambda_inv = adiag(Lambda_inv, solve( new_Lambda ))
        
      }
      
      
    }
    
  }
  print(alpha_0)
  each_chain_task = function(i) {
    #Set initial values for parameters.
    
    sigma = runif(1)
    if(init == "Random"){
      sigma1 = runif(K)
      sigma0 = sigma1 / ratio
      p = runif(K, min = 0, max = 0.1)
    }
    else{
      sigma1 = sigma_1
      sigma0 = sigma_0 
      p = p
    }
    B = random_upper_tri(K)
    
    
    
    if(is.null(num_traits)){
      if( all(cor_mat == diag(K))){
        this_output = gibbs_sampler(Gamma_hat, Sd_hat, N, B, sigma, sigma1,
                                    sigma0, p, A, Z,
                                    alpha_B, beta_B, alpha_0, alpha_1, beta_0, beta_1, a, b)
      }
      else{
        this_output = gibbs_sampler_with_corr(Gamma_hat,
                                              Sd_hat, cor_mat, N, B,
                                              sigma, sigma1,
                                              sigma0, p, A, Z,
                                              alpha_B, beta_B, 
                                              alpha_0, alpha_1, beta_0,
                                              beta_1, a, b, Lambda, Lambda_inv)
      }
      
    }
    else{
      this_output = gibbs_sampler_new(Gamma_hat,
                                      Sd_hat, cor_mat, N, B,
                                      sigma, sigma1,
                                      sigma0, p, A, Z,
                                      alpha_B, beta_B, 
                                      alpha_0, alpha_1, beta_0,
                                      beta_1, a, b, Lambda, Lambda_inv, num_traits)
    }
    
    
    this_output
    
  }
  
  
  print(Sys.time())
  if (parallel){
    outputs = foreach(i=1:chains) %dopar% each_chain_task(i)
  } else {
    outputs = lapply(1:1, each_chain_task)
  }
  print(Sys.time())
  
  return (outputs)
}
