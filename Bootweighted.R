source("Bootstrapping.r")

## Bootstrapping weighted matrix
## with restriction and without restriction

library(compiler)

Bootweighted = function(weighted_mtx){
  
  M = matrix(0,nrow(weighted_mtx),ncol(weighted_mtx))
  mcat = sort(unique(as.vector(weighted_mtx)))
  if(mcat[1]!=0){mcat = c(0,mcat)}
  
  for(j in 2:length(mcat)){
    
    c = which(weighted_mtx>=mcat[j])
    B = matrix(0,nrow(weighted_mtx),ncol(weighted_mtx))
    B[c] = 1
    boot = Bootbinary(B)$Matrix*(mcat[j]-mcat[j-1])
    M = M+boot
    
  }#end of j loop
  return(M)
}

Bootweighted_restrict = function(weighted_mtx,
                                 restr_etry){
  
  M = matrix(0,nrow(weighted_mtx),ncol(weighted_mtx))
  for(k in 1:nrow(restr_etry)){
    weighted_mtx[restr_etry[k,1],restr_etry[k,2]] = 0
  }
  mcat = sort(unique(as.vector(weighted_mtx)))
  if(mcat[1]!=0){mcat = c(0,mcat)}
  
  for(j in 2:length(mcat)){
    
    c = which(weighted_mtx>=mcat[j])
    B = matrix(0,nrow(weighted_mtx),ncol(weighted_mtx))
    B[c] = 1
    boot = Bootbinary_restrict(B,restr_etry)$Matrix*(mcat[j]-mcat[j-1])
    M = M+boot
    
  }#end of j loop
  return(M)
  
}

Bootweighted = cmpfun(Bootweighted)
Bootweighted_restrict = cmpfun(Bootweighted_restrict)

