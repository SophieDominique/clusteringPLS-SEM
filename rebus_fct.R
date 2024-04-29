#-------------------------------------------------------------------------------
# REBUS-PLS
#-------------------------------------------------------------------------------
rebusPLS <- function(data, pls_matrix, wscheme, nb_LV, nb_MV, tol = 1e-5, iter.pls = 100, 
                     iter = 100, stop.crit = 0.05, nb_cluster = -1)
{
  #-----------------------------------------------------------------------------
  # pls_matrix must be a lower triangular matrix
  # there is no choice of the mode because rebus-pls method is only available 
  # when all blocks are in mode A
  # nb_cluster = -1 : automatic selection of the number of classes from the CAH
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  # PLSPM on all individuals
  # source("algo_lohmoller.R") #own function
  
  N = nrow(X)
  res_pls = lohmProc(X, pls_matrix, "A", wscheme, nb_LV, nb_MV, tol, iter.pls)
  
  #-----------------------------------------------------------------------------
  rebus_score = res_pls$data %*% res_pls$outerW
  cor.XY = cor(res_pls$data, rebus_score)
  
  # sign ambiguity
  ODM = res_pls$outerW
  ODM[res_pls$outerW != 0] = 1
  w_sign = sign(colSums(sign((cor.XY * ODM))))
  if (any(w_sign <= 0)) {
    w_sign[w_sign == 0] = -1
    # scores
    rebus_score = rebus_score %*% diag(w_sign, ncol(res_pls$outerW), ncol(res_pls$outerW))
  }
  rebus_crossload = cor(res_pls$data, rebus_score)
  rebus_load = rebus_crossload
  rebus_load[res_pls$outerW == 0] = 0
  rebus_path = path_coeff(rebus_score, nb_LV, t(pls_matrix))
  
  
  #-----------------------------------------------------------------------------
  # compute res
  mres = res_mesure(res_pls$data, rebus_load, rebus_score)
  sres = res_struct(rebus_score, pls_matrix, rebus_path)
  
  #-----------------------------------------------------------------------------
  # compute CAH
  res = cbind(mres, sres)
  res_clus = HCPC(as.data.frame(res), nb.clust=nb_cluster, graph = FALSE, consol = FALSE)
  nb_class = length(levels(res_clus$data.clust[,"clust"]))
  CM = matrix(0, nrow = N, ncol = nb_class)
  old.class = as.numeric(as.character(res_clus$data.clust[,"clust"]))

  #-----------------------------------------------------------------------------
  # initialisation matrix
  new.class = rep(NA, N)
  gof_grp = matrix(NA, nrow = iter, ncol = 3*nb_class)
  all_gqi = rep(NA, iter)
  class1 = rep(1,N)
  gqi1 = computeGQI(X, pls_matrix, 1, class1, res_pls, nb_LV, nb_MV)$GQI
  gqi_mes = rep(NA, iter)
  gqi_struc = rep(NA, iter)
  
  #-----------------------------------------------------------------------------
  # loop REBUS
  for(i in 1:iter)
  {
    for(k in 1:nb_class)
    {
      ind_k = which(old.class == k)
      res_pls_k = lohmProc(X[ind_k,], pls_matrix, "A", wscheme, nb_LV, nb_MV, tol, iter.pls)
      
      # compute res on all units
      nk = nrow(X[ind_k,])
      mean.k = apply(X[ind_k,], 2, mean)
      sd.k = sqrt((nk-1)/nk) * apply(X[ind_k,], 2, sd)
      X.k = scale(X, center=mean.k, scale=sd.k)
      allScore = X.k %*% res_pls_k$outerW
      
      #-----------------------------------------------------------------------------
      mres_k = res_mesure(X.k, res_pls_k$loading, allScore)
      sres_k = res_struct(allScore, pls_matrix, res_pls_k$pathCoeff)
      
      # compute CM
      CM[,k] = compute_CM(N, nk, mres_k, sres_k, res_pls_k$loading, allScore[ind_k,], pls_matrix)
    }
    nbIterconv = i
    
    #-----------------------------------------------------------------------------
    # set groups
    new.class = apply(CM, 1, which.min)

    # check convergence
    diff_class = new.class - old.class 
    unit_change = length(which(diff_class != 0))
    
    # rate of unit change
    if (unit_change/N <= stop.crit)
      break
    
    old.class = new.class
    
    if (any(table(new.class) <= 5) || length(table(new.class)) != nb_class) 
      stop("\nToo few units: a class with less than 6 units was detected")
    
    #-----------------------------------------------------------------------------
    # compute quality index model
    nb_endo = length(which(rowSums(pls_matrix) != 0))
    model_k = as.list(1:nb_class)
    com_model = matrix(NA, nrow = nb_LV, ncol = nb_class)
    rownames(com_model) = paste("Com",c(1:nb_LV), sep="")
    R2_model = matrix(NA, nrow = nb_endo, ncol = nb_class)
    rownames(R2_model) = paste("R2", which(rowSums(pls_matrix) != 0), sep=".")
    red_model = matrix(NA, nrow = nb_endo, ncol = nb_class)
    rownames(red_model) = paste("Red", which(rowSums(pls_matrix) != 0), sep="")
    
    for(k in 1:nb_class)
    {
      ind_model = which(new.class == k)
      
      model_k[[k]] = lohmProc(X[ind_model,], pls_matrix, "A", wscheme, nb_LV, nb_MV, tol, iter.pls)
      
      res_quality = index_quality(model_k[[k]], nb_LV, nb_MV)
      com_model[,k] = res_quality$commu
      R2_model[,k] = res_quality$rsq
      red_model[,k] = res_quality$redun
      
      gof_grp[i,3*k-2] = weighted.mean(com_model[,k], nb_MV)
      gof_grp[i,3*k-1] = mean(R2_model[,k])
      gof_grp[i,3*k] = sqrt(weighted.mean(com_model[,k], nb_MV) * mean(R2_model[,k]))
      
    }
    res_gqi = computeGQI(X, pls_matrix, nb_class, new.class, model_k, nb_LV, nb_MV)
    all_gqi[i] = res_gqi$GQI
    gqi_mes[i] = res_gqi$commu
    gqi_struc[i] = res_gqi$rsq
    
  }
  
  #-----------------------------------------------------------------------------
  # final compute quality index model
  nb_endo = length(which(rowSums(pls_matrix) != 0))
  model_k = as.list(1:nb_class)
  com_model = matrix(NA, nrow = nb_LV, ncol = nb_class)
  rownames(com_model) = paste("Com",c(1:nb_LV), sep="")
  R2_model = matrix(NA, nrow = nb_endo, ncol = nb_class)
  rownames(R2_model) = paste("R2", which(rowSums(pls_matrix) != 0), sep=".")
  red_model = matrix(NA, nrow = nb_endo, ncol = nb_class)
  rownames(red_model) = paste("Red", which(rowSums(pls_matrix) != 0), sep="")
  GoF = rep(NA, nb_class)
  
  for(k in 1:nb_class)
  {
    ind_model = which(new.class == k)
    model_k[[k]] = lohmProc(X[ind_model,], pls_matrix, "A", wscheme, nb_LV, nb_MV, tol, iter.pls)
    res_quality = index_quality(model_k[[k]], nb_LV, nb_MV)
    com_model[,k] = res_quality$commu
    R2_model[,k] = res_quality$rsq
    red_model[,k] = res_quality$redun
    GoF[k] = sqrt(weighted.mean(com_model[,k], nb_MV) * mean(R2_model[,k]))
    
    gof_grp[i,3*k-2] = weighted.mean(com_model[,k], nb_MV)
    gof_grp[i,3*k-1] = mean(R2_model[,k])
    gof_grp[i,3*k] = GoF[k]
  }
  
  res_gqi = computeGQI(X, pls_matrix, nb_class, new.class, model_k, nb_LV, nb_MV)
  gqi = res_gqi$GQI
  all_gqi[i] = res_gqi$GQI
  gqi_mes[i] = res_gqi$commu
  gqi_struc[i] = res_gqi$rsq
  
  allIndex = rbind(com_model, red_model, R2_model, GoF)
  
  #-----------------------------------------------------------------------------
  return(list(model = model_k, quality = allIndex, segment = new.class, GQI = gqi, 
              nbIterReb = nbIterconv, all.gqi =all_gqi, all.gof = gof_grp,
              GQI_mes = gqi_mes, GQI_struc = gqi_struc, nb_class = nb_class, CM_index = CM) )
  
  #-----------------------------------------------------------------------------
}



#-----------------------------------------------------------------------------
# functions for index quality 
#-----------------------------------------------------------------------------
computeGQI= function(X, pls_matrix, nb_class, segment, model_k, nb_LV, nb_MV)
{
  nb_endo = length(which(rowSums(pls_matrix) != 0))
  gqi.local = rep(NA, nb_class)
  prop.unit = table(segment) / length(segment)
  tempC = rep(NA, nb_class)
  tempR = rep(NA, nb_class)
  
  for(k in 1:nb_class)
  {
    if(nb_class == 1)
    {
      res_quality = index_quality(model_k, nb_LV, nb_MV)
    }
    else
    {
      res_quality = index_quality(model_k[[k]], nb_LV, nb_MV)
    }
    
    com_model = weighted.mean(res_quality$commu, nb_MV)
    R2_model = mean(res_quality$rsq)
    tempC[k] = com_model*prop.unit[k]
    tempR[k] = R2_model*prop.unit[k]
    
  }
  
  gqi = sqrt(sum(tempC)*sum(tempR))
  
  return(list(commu = sum(tempC), rsq = sum(tempR), GQI = gqi))
  
}

#-----------------------------------------------------------------------------
# commulnality, redundandy, R2
#-----------------------------------------------------------------------------
index_quality = function(model_pls, nb_LV, nb_MV)
{
  isEndo = which(rowSums(model_pls$path_model) != 0)
  com = rep(NA, nb_LV)
  sumMV = cumsum(nb_MV)
  mv1 = 1
  
  # communality
  for(j in 1:nb_LV)
  {
    mv2 = sumMV[j]
    
    if(mv1 == mv2) # cas o? une seule MV
    {
      load_com = t(as.matrix(model_pls$loading[mv2,]))
    }
    else
    {
      load_com = model_pls$loading[mv1:mv2,]
    }
    
    com[j] = mean(load_com^2 %*% matrix(1, nrow = ncol(load_com), ncol = 1))
    
    mv1 = mv2 +1
  }
  # com = model_pls$loading^2 %*% matrix(1, nrow = ncol(model_pls$loading), ncol = 1)
  # R-square
  R2 = computeR2(model_pls$score, model_pls$path_model)
  # redundancy
  red = com[isEndo] * R2
  
  return(list(commu = com, rsq = R2, redun = red))
  
}

#-----------------------------------------------------------------------------
# compute measurement model residual
#-----------------------------------------------------------------------------
res_mesure = function(X, loading, score)
{
  X_chap = score%*%t(loading)
  return(X-X_chap)
}

#-----------------------------------------------------------------------------
# compute structural model residual
#-----------------------------------------------------------------------------
res_struct = function(score, pls_matrix, pathCoeff)
{
  isEndo = as.logical(rowSums(pls_matrix))
  
  if(sum(isEndo) != 1)
  {
    nu = score%*%t(pathCoeff[isEndo == 1,])
  }
  else
  {
    nu = score%*%pathCoeff[isEndo == 1,]
  }

  return(score[,isEndo==1]-nu)
}

#-----------------------------------------------------------------------------
# compte R2 on LV endo
#-----------------------------------------------------------------------------
computeR2 = function(score, pls_matrix)
{
  isEndo = as.logical(rowSums(pls_matrix))
  nbEndo = sum(isEndo)
  rsq = rep(0, nrow(pls_matrix))
  
  for(j in 1:nbEndo)
  {
    # index endo LV
    j1 = which(isEndo)[j]
    # index pred LV
    j2 = which(pls_matrix[j1,] == 1)
    
    regEndo = summary(lm(score[,j1]~score[,j2]))
    rsq[j1] = regEndo$r.squared
  }
  
  return(rsq[isEndo])
}


#-----------------------------------------------------------------------------
# compute clossness measure
#-----------------------------------------------------------------------------
compute_CM = function(N, nk, mres, sres, loading, score, pls_matrix)
{
  # mesure part
  # communality
  com = (loading^2) %*% matrix(1, nrow = ncol(loading), ncol = 1) %*% matrix(1, nrow = 1, ncol = N)
  num_mes = rowSums(mres^2 / t(com))
  denom_mes = sum(num_mes) / (N-2)
  part_mes = num_mes / denom_mes
  
  # struct part
  # R2
  rsq = computeR2(score, pls_matrix)

  if(length(rsq) == 1)
  {
    num_struct = rowSums(sres^2%*%(1/rsq))
  }
  else
  {
    num_struct = rowSums(sres^2%*%diag(1/rsq))
  }
  
  denom_struct = sum(num_struct) / (N-2)
  part_struct = num_struct / denom_struct
  
  # CM for all i in class k
  return(sqrt(part_mes * part_struct))
}

#-----------------------------------------------------------------------------
# compute path coefficient
#-----------------------------------------------------------------------------
path_coeff = function(score, nb_LV, pls_matrix)
{
  pC = matrix(0, ncol = nb_LV, nrow = nb_LV)
  fun_pred = function(x) which(x == 1)
  pred = apply(pls_matrix, 2, fun_pred)
  for(i in 1:nb_LV)
  {
    # path coeff
    if(length(pred[[i]]) == 0) next
    pC[i,pred[[i]]] = solve(cor(as.matrix(score[,pred[[i]]])))%*%
      cor(score[,pred[[i]]],score[,i])
  }
  
  return(pC)
  
}

