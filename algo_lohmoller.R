#-------------------------------------------------------------------------------
# PLS-SEM : procedure lohmoller
#-------------------------------------------------------------------------------
lohmProc <- function(X, A, mode, wscheme, nb_LV, nb_MV, tol = 1e-5, iter = 100)
{
  #-----------------------------------------------------------------------------
  # input X, A, mode, wscheme
  # A = adjacency matrix
  # mode = A or B
  # wscheme = centroid, factor, path
  # nb_LV = nb of LV
  # nb_MV = vector with nb of MV in each bloc
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  # check param
  if(!(mode %in% c("A","B"))){
    return("Erreur : mode must be one of the following: \"A\", \"B\" ")
  }
  if(!(wscheme %in% c("centroid","factorial","path"))){
    return("Erreur : wscheme must be one of the following: \"centroid\", \"factorial\", \"path\" ")
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
    eigVect = eigen(cor(as.matrix(X[,MVsum[i]:(MVsum[i+1]-1)])))$vectors[,1]
    new_WM[MVsum[i]:(MVsum[i+1]-1),i] = eigVect
  }

  # standardize w and score
  old_score = scale(X%*%new_WM)
  new_score = old_score

  #-------------------------------------------------------------------------------
  # PLS loop
  #-------------------------------------------------------------------------------
  # inner estimation 
  for(i in 1:iter)
  {
    # compute R*C
    R = cor(old_score)
    C = A+t(A)

    if(wscheme == "centroid")
    {
      theta = sign(R)*C
    }
    else if(wscheme == "factorial")
    {
      theta = R*C
    }
    else if(wscheme == "path")
    {
      # compute P
      P = path_coeff(old_score,nb_LV,A)
      theta = t(P) + t(A)*R
    }
    # compute z_tilde
    z_tilde = old_score%*%theta
    
    #-------------------------------------------------------------------------------
    # outer estimation
    
    # compute w_tilde
    for(k in 1:nb_LV)
    {
      sfact_k = as.matrix(z_tilde[,k])
      mf_k = as.matrix(X[,MVsum[k]:(MVsum[k+1]-1)])

      new_WM[MVsum[k]:(MVsum[k+1]-1),k] = cor(mf_k,sfact_k)
      if(mode == "B")
      {
        new_WM[MVsum[k]:(MVsum[k+1]-1),k] = solve(cor(mf_k))%*%new_WM[MVsum[k]:(MVsum[k+1]-1),k]
      }
    }

    #-------------------------------------------------------------------------------
    # standardize w_tilde
    new_score = scale(X%*%new_WM)

    #-------------------------------------------------------------------------------
    # stop criterion
    crit = sum((new_score - old_score)^2)
    conv = (crit < tol)

    if(conv || i == iter)
    {
      break
    }

    old_score = new_score
  }
  
  #-------------------------------------------------------------------------------
  # restimation of param
  # path coeff
  pC = path_coeff(new_score, nb_LV, t(A))
  
  # outer weights
  outerW = matrix(0, ncol = nb_LV, nrow = sum(nb_MV))
  for(k in 1:nb_LV)
  {
    outerW[MVsum[k]:(MVsum[k+1]-1),k] = solve(cor(as.matrix(X[,MVsum[k]:(MVsum[k+1]-1)])))%*%
      cor(as.matrix(X[,MVsum[k]:(MVsum[k+1]-1)]),new_score[,k])
  }
  
  cross_loading = cor(X, new_score)
  loading = cross_loading
  loading[outerW == 0] = 0

  #-------------------------------------------------------------------------------
  return(list(nb_iter = i, score = new_score, last_w = new_WM, outerW = outerW, innerW = theta,
              data = X, pathCoeff = pC, cross_loading = cross_loading, loading = loading, path_model = A))
  
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
  colnames(pC) = colnames(A)
  rownames(pC) = rownames(A)
  
  #-------------------------------------------------------------------------------
  return(pC)
  
}