source("EnergyOptim.r")
library(compiler)

## this is a function for bootstrapping an adjancency matrix 
## with constraint of its row and column sums.
## return: (1) the energy of the bootstrapping adj matrix (2) the new Adj matrix

pij = function(da,db,dahat,dbhat,m,alpha){
  p = (dahat*dbhat)^alpha*(1-((da*db)/(4*m)))
  return(p)
}

Bootbinary = function(Adj){
  
  Adj = as.matrix(Adj)
  
  m = mean(Adj)
  if(m==0 || m==1){
    Bootmatrix = Adj
  }else{
    
    Bootmatrix = matrix(0,nrow(Adj),ncol(Adj))
    ii=1
    
    while(sum(abs(rowSums(Bootmatrix)-rowSums(Adj)))>0 & 
            sum(abs(colSums(Bootmatrix)-colSums(Adj)))>0){
      
      ######
	#print(ii)

      if(ii==1000){Bootmatrix = Adj; break}
      da = rowSums(Adj)
      db = colSums(Adj)
      
      Na = nrow(Adj)
      Nb = ncol(Adj)
      
      m = (sum(da)+sum(db))/2
      
      E = c()
      dahat = da
      dbhat = db
      P = 1
      C = expand.grid(1:Na,1:Nb)
      #alpha = 1
      
      for (i in 1:(Na*Nb)){
        
        diff = abs(dahat[C[,1]]-dbhat[C[,2]])
        
        alpha = rep(1,nrow(C)) #initial of alpha
        alpha[which(diff>quantile(diff,0.7))] = 2
        #print(table(alpha))
        
        # calculate probability pij
        Pr = pij(da[C[,1]],db[C[,2]],dahat[C[,1]],dbhat[C[,2]],m,alpha)
        ## weird...
        if(sum(Pr)==0) break
        
        Pr = Pr/sum(Pr)
        
        idx = which(rmultinom(1,1,Pr)==1)
        P = P*Pr[idx]
        dahat[C[idx,1]] = dahat[C[idx,1]]-1
        dbhat[C[idx,2]] = dbhat[C[idx,2]]-1
        
        E = rbind(E,C[idx,])
        
        C = C[-idx,]
        
        # break the for loop when there is no edges and add in
        if((sum(dahat)+sum(dbhat))==0) break
        
      }
      
      Bootmatrix = matrix(0,nrow(Adj),ncol(Adj))
      for(j in 1:nrow(E)){
        Bootmatrix[E[j,1],E[j,2]] = 1
      }
      
      ii = ii+1
    }# end of while loop
    
    
  } # End of if
  
  #print(ii)
  Energy = GetBipEnergy(Bootmatrix)
  
  return(list(Energy=Energy, Matrix=Bootmatrix))
}

Bootbinary = cmpfun(Bootbinary)
