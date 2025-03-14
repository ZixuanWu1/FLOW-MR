#' Total effect
#' 
#' Estimate the total effect of an exposure on the outcome
#' 
#' @param results results from gibbs_wrapper
#' @param K number of phenotypes
#' @param warump The length of warmup period. Default is 3000
#' @param quantiles Quantiles of interests
#'  
#' @return A K by K by 4 array of summary of total effect.
#' 
#' @export

total_effect = function(results, K, warmup = 3000, quantiles = c(0.025, 0.5, 0.975)) {
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
  r = length(quantiles)
  df = array(dim = c(K, K, r + 1))
  for(i in 1:K){
    for(j in 1:K){
      df[i, j, 1] = mean(effects[i,j, ])
      
      for(rr in 1:r){
        df[i, j, rr + 1] = quantile(effects[i, j, ], quantiles[rr])
        
      }
      
      
    }
    
  }
  dimnames(df)[[3]] = c("Mean", paste("Quantile", quantiles, sep="_"))
  
  return(df)
  
  
}


#' Indirrect effedt
#' 
#' Estimate the indirrect effect of an exposure on the outcome
#' 
#' @param results results from gibbs_wrapper
#' @param K number of phenotypes
#‘ @param path A vector representing the pathway of interests
#' @param quantiles A vector of quantiles of interests
#' @param warump The length of warmup period. Default is 3000
#' 
#' @return A K by K by 4 array of summary of direct effect.
#' 
#' @export
indirect_effect = function(results, K,  warmup = 3000, path = NULL, quantiles =  c(0.025, .5, 0.975)) {
  
  results = result_process(results, K)
  m = length(results)
  len = utils::tail(dim(results[[1]][["B"]]), n = 1)
  
  if(is.null(path)){
    effects = matrix(0, K, m * (len - warmup))
    for(i in 1:m){
      for(j in 1:(len - warmup)){
        temp_B = results[[i]][["B"]][,,(j + warmup)]
        effects[,( (i - 1) * (len - warmup) + j )] = (compute_total(temp_B) - temp_B)[1,]
      }
    }
    df = matrix(nrow = K, ncol = 4)
    for(i in 1:K){
      df[i, 1] = mean(effects[i,])
      df[i, 2:4] = quantile(effects[i,], quantiles)
      
    }
    df = as.data.frame(df)
    for(i in 1:K){
      row.names(df)[i] = paste("Exp", toString(i))
    }
    
    
  }  else if (all(path == "all")) {
    
    effects = array( dim = c(K , K, m * (len - warmup) ))
    for(i in 1:m){
      for(j in 1:(len - warmup)){
        temp_B = results[[i]][["B"]][,,(j + warmup)]
        effects[,, ( (i - 1) * (len - warmup) + j )] = (compute_total(temp_B) - temp_B)
      }
    }
    r = length(quantiles)
    df = array(dim = c(K, K, r + 1))
    for(i in 1:K){
      for(j in 1:K){
        df[i, j, 1] = mean(effects[i,j, ])
        
        for(rr in 1:r){
          df[i, j, rr + 1] = quantile(effects[i, j, ], quantiles[rr])
          
        }

        
      }
      
    }
    dimnames(df)[[3]] = c("Mean", paste("Quantile", quantiles, sep="_"))
    
    return(df)
    
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
    df = c(mean(effects), quantile(effects, quantiles)  )
    df = as.data.frame(matrix(df, nrow= 1))
    
  }
  colnames(df)[1] = "Mean"
  for(j in 2:ncol(df)){
    colnames(df)[j] = paste(toString(quantiles[j - 1] * 100), "%", sep = "")
  }
  
  return(df)
  
}

#' @keywords internal
#'
compute_total <- function(B){
  K = length(B[,1])
  return( solve(diag(K) - B)  - diag(K) )
}
