#-------------------------------------------------------------------------------
# PLS-SEM : procedure hanafi-wold
#-------------------------------------------------------------------------------
woldProc <- function(X, A, mode, wscheme, nb_LV, nb_MV, tol = 1e-5)
{
  #-------------------------------------------------------------------------------
  # input X, A, mode, wscheme
  # mode = A or B
  # wscheme = centroid, factor, path
  # nb_LV = nb of LV
  # nb_MV = vector with nb of MV in each bloc
  #-------------------------------------------------------------------------------
  # check param
  if(!(mode %in% c("A","B"))){
    return("Erreur : mode must be one of the following: \"A\", \"B\" ")
  }
  if(!(wscheme %in% c("centroid","factor","path"))){
    return("Erreur : wscheme must be one of the following: \"centroid\", \"factor\", \"path\" ")
  }
  if(length(nb_MV) != nb_LV)
  {
    return("Erreur : length of nb_MV mist be equal to the nb of LV")
  }
  if(sum(nb_MV) != dim(X)[2])
  {
    return("Erreur : the total nb of nb of MV is not equal to the nb of columns of X")  
  }
  conv = FALSE
  
  #-------------------------------------------------------------------------------
  # centrer X
  mean.X = apply(X, 2, mean)
  sd.X = sqrt((nrow(X)-1)/nrow(X)) * apply(X, 2, sd)
  X = scale(X, center = mean.X, scale=sd.X)
  
  #-------------------------------------------------------------------------------
  # initialization 
  # init w
  MVsum = c(0,cumsum(nb_MV))+1
  new_WM = matrix(0, nrow = sum(nb_MV), ncol = nb_LV)
  for(i in 1:nb_LV)
  {
    new_WM[MVsum[i]:(MVsum[i+1]-1),i] = rep(1,nb_MV[i])
  }
  # standardize w and score
  old_score = scale(X%*%new_WM)
  new_score = old_score
  
  #-------------------------------------------------------------------------------
  # PLS loop
  i = 0
  R = matrix(1,ncol=nb_LV, nrow = nb_LV)
  C = t(A) + A
  
  #-------------------------------------------------------------------------------
  # inner estimation 
  while(!conv)
  {
    i = i+1
    # estimation for the block X_k0
    for(k in 1:nb_LV)
    {
      ## compute correlation
      if(k!=1)
      {
        for(k1 in 1:(k-1))
        {
          R[k,k1] = cor(old_score[,k],new_score[,k1])
        }
      }
      
      if(k!=nb_LV)
      {
        for(k2 in (k+1):nb_LV)
        {
          R[k,k2] = cor(old_score[,k],old_score[,k2])
        }
      }
      
      ## compute theta
      if(wscheme == "centroid")
      {
        theta = sign(R)*C
      }
      else if(wscheme == "factor")
      {
        theta = R*C
      }
      # else if(wscheme == "path")
      # {
      #   # compute P
      #   P = path_coeff(old_score,nb_LV,A)
      #   theta = P + t(A)*R
      # }
      
      #-------------------------------------------------------------------------------
      ## compute score for block k0
      z_tilde = 0
      if(k!=1)
      {
        z_tilde = new_score[,1:(k-1)]%*%as.matrix(theta[1:(k-1),k])
      }
      if(k!=nb_LV)
      {
        z_tilde = z_tilde + new_score[,(k+1):nb_LV]%*%as.matrix(theta[(k+1):nb_LV,k])
      }
      
      #-------------------------------------------------------------------------------
      # outer estimation 
      # compute w_tilde
      sfact_k = as.matrix(z_tilde)
      mf_k = X[,MVsum[k]:(MVsum[k+1]-1)]

      new_WM[MVsum[k]:(MVsum[k+1]-1),k] = t(mf_k)%*%sfact_k
      if(mode == "B")
      {
        new_WM[MVsum[k]:(MVsum[k+1]-1),k] = solve(cor(mf_k))%*%new_WM[MVsum[k]:(MVsum[k+1]-1),k]
      }
      
      new_score[,k] = scale(mf_k%*%new_WM[MVsum[k]:(MVsum[k+1]-1),k])
      
    }
    
    #-------------------------------------------------------------------------------
    # stop criterion
    crit = sum((new_score - old_score)^2)
    conv = (crit < tol)
    
    old_score = new_score
  }
  
  #-------------------------------------------------------------------------------
  # path coeff
  pC = path_coeff(new_score, nb_LV, t(A))
  
  # outer weights
  outerW = matrix(0, ncol = nb_LV, nrow = sum(nb_MV))
  for(k in 1:nb_LV)
  {
    outerW[MVsum[k]:(MVsum[k+1]-1),k] = solve(cor(X[,MVsum[k]:(MVsum[k+1]-1)]))%*%
      cor(X[,MVsum[k]:(MVsum[k+1]-1)],new_score[,k])
  }
  
  cross_loading = cor(X, new_score)
  loading = cross_loading
  loading[outerW == 0] = 0
  
  #-------------------------------------------------------------------------------
  return(list(nb_iter = i, score = new_score, outerW = outerW, innerW = theta, 
              data = X, pathCoeff = pC, cross_loading = cross_loading, loading = loading))
  
}


#-------------------------------------------------------------------------------
# compute path coefficient with OLS reg
#-------------------------------------------------------------------------------
path_coeff = function(score, nb_LV, A)
{
  pC = matrix(0, ncol = nb_LV, nrow = nb_LV)
  fun_pred = function(x) which(x == 1)
  pred = apply(A, 2, fun_pred)
  for(i in 1:nb_LV)
  {
    # path coeff
    if(length(pred[[i]]) == 0) next
    pC[i,pred[[i]]] = solve(cor(as.matrix(score[,pred[[i]]])))%*%
      cor(score[,pred[[i]]],score[,i])
  }
  
  return(pC)
  
}