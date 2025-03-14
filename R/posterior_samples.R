#' Posterior samples
#' 
#' Estimate the total effect of an exposure on the outcome
#' 
#' @param raw raw results 
#' @param par parameter of interets
#' @param ind index of the parameter
#' @param warump The length of warmup period. Default is 3000
#'  
#' @return A vector of posterior samples from all chains
#' 
#' @export

get_posterior_samples  = function(raw, par, ind = NULL, warmup = 3000) {
  m = length(raw)
  raw = result_process(raw, K)
  
  len = utils::tail(dim(raw[[1]][["B"]]), n = 1)
  
  
  if (par == "B"){
    if (length(ind) != 2){
      warning( paste( "incorrect index dimension for", par))
      return(NULL)
    }
    samples = c()
    for(i in 1:m){
      samples = c(samples, raw[[i]][["B"]][ind[1], ind[2], ((warmup + 1):len)])
    }
    
  } else if (par == "sigma1" | par == "sigma0" | par == "p"){
    if (length(ind) != 1){
      warning( paste( "incorrect index dimension for", par))
      return(NULL)
    }
    samples = c()
    for(i in 1:m){
      samples = c(samples, raw[[i]][[par]][ind, ((warmup + 1):len)])
    }
  } else if (par == "sigma") {
    if (length(ind) != 0){
      warning( paste( "incorrect index dimension for", par))
      return(NULL)
    }
    samples = c()
    for(i in 1:m){
      samples = c(samples, raw[[i]][[par]][((warmup + 1):len)])
    }
  } else{
    warning("Parameter must be one of B, sigma1, sigma0, sigma, p")
  }
  
  return(samples)
  
  
}
















#' Posterior samples of total effects
#' 
#' Estimate the total effect of an exposure on the outcome
#' 
#' @param results results from gibbs_wrapper
#' @param K number of phenotypes
#' @param ind index of interet
#' @param warump The length of warmup period. Default is 3000
#'  
#' @return A K by K by 4 array of summary of total effect.
#' 
#' @export

get_total_samples = function(results, K, ind, warmup = 3000) {
  
  if (! length(ind) == 2){
    warning("incorrect ind dimension")
    return(NULL)
  }
  
  results = result_process(results, K)
  m = length(results)
  len = utils::tail(dim(results[[1]][["B"]]), n = 1)
  effects = array( dim = c(K , K, m * (len - warmup) ))
  for(i in 1:m){
    for(j in 1:(len - warmup)){
      temp_B = results[[i]][["B"]][,,(j + warmup)]
      effects[,, ( (i - 1) * (len - warmup) + j )] = (compute_total(temp_B))
    }
  }
  
  return(effects[ind[1], ind[2], ])
  
  
}


#' Posterior samples of indirect effects
#' 
#' Estimate the indirrect effect of an exposure on the outcome
#' 
#' @param results results from gibbs_wrapper
#' @param K number of phenotypes
#' @param warmup The length of warmup period. Default is 3000
#' @param path A vector representing the pathway of interests
#' @param ind index of indirect effect of interest if path is NULL
#' 
#' @return A vector of posterior samples
#' 
#' @export
get_indirect_samples = function(results, K,  warmup = 3000, path = NULL, ind = NULL) {
  
  results = result_process(results, K)
  m = length(results)
  len = utils::tail(dim(results[[1]][["B"]]), n = 1)
  
  if (all(path == "all")) {
    
    effects = array( dim = c(K , K, m * (len - warmup) ))
    for(i in 1:m){
      for(j in 1:(len - warmup)){
        temp_B = results[[i]][["B"]][,,(j + warmup)]
        effects[,, ( (i - 1) * (len - warmup) + j )] = (compute_total(temp_B) - temp_B)
      }
    }
    
    return(effects[ind[1], ind[2], ])
    
  } else{
    effects = NULL
    for(i in 1:m){
      for(j in 1:(len - warmup)){
        temp_B = results[[i]][["B"]][,,(j + warmup)]
        cur = 1
        for (ele in 1:(length(path) - 1)){
          cur = cur * temp_B[path[ele + 1], path[ele]]
        }
        effects = c(effects, cur)
      }
    }
    
    
  }
  
  return(effects)
  
}
