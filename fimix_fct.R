fimixPLS <- function(data, pathMat, blocksPLS, modesPLS, wscheme = "path", nbCluster, itFIM = 500, 
                            tolFIM = 1e-7,
                            tolPLS = 1e-10, itPLS = 300,
                            pathMatLOC = NULL,
                            withConst = FALSE, init_pik = NULL, init_scores = NULL)
{
  #----------------------------------------------------------------------------------------
  # pathMat = adjacent matrix of the global model
  # pathMatLOC = adjacent matrix of the local model
  # wscheme = c("centroid", "factorial", "path")
  library(plspm)
  library(MASS)
  #----------------------------------------------------------------------------------------
  # PLS-SEM model to compute scores and common path coefficients
  resPLS = plspm(data, path_matrix = pathMat, blocks = blocksPLS, modes = modesPLS,
                 scheme = wscheme, tol = tolPLS, maxiter = itPLS)
 
  library(FactoMineR)
  if(is.null(init_scores))
  {
    scoresPLSPM = resPLS$scores
  } else {
    scoresPLSPM = init_scores
  }
  #----------------------------------------------------------------------------------------
  # intialisation partition
  N = nrow(data)
  
  if(is.null(init_pik))
  {
    pik = matrix(runif(N*nbCluster), nrow = N, ncol = nbCluster)
    pik = sweep(pik, 1, rowSums(pik), "/")
  } else {
    pik = init_pik
  }
  seg = apply(pik, 1, which.max)
  
  #----------------------------------------------------------------------------------------
  # param
  nbArrow = sum(pathMat)
  
  if(!is.null(pathMatLOC))
  {
    if(all(pathMatLOC == 0))
    {
      print("Error : the model matrix is null.")
      return(NA)
    }
    
    clustLocal = TRUE
    segLocal = as.vector(t(pathMatLOC))[which(as.vector(t(pathMat)) == 1)]
  }  else  {
    clustLocal = FALSE
    pathMatLOC = pathMat
    segLocal = rep(1,nbArrow)
  }
  
  nbArrowLoc = sum(pathMatLOC)
  
  #----------------------------------------------------------------------------------------
  tauGlobal = c(t(resPLS$path_coefs))
  tauGlobal = tauGlobal[-which(tauGlobal==0)]
  
  #----------------------------------------------------------------------------------------
  isEndo = which(rowSums(pathMat) != 0)
  nbEndo = length(isEndo)
  isEndoLoc = which(apply(pathMatLOC,1,sum) != 0)
  nbEndoLoc = length(isEndoLoc)
  isExoPur = which(rowSums(pathMat) == 0)
  nbExoPur = length(isExoPur)
  nbLV = ncol(pathMat)
  
  #----------------------------------------------------------------------------------------
  oldLn = -1e10
  oldLnC = -1e10
  devLnL = matrix(NA, nrow = itFIM, ncol = 2)
  devLnLC = matrix(NA, nrow = itFIM, ncol = 2)
  tauLoc = matrix(0, nrow = nbArrowLoc, ncol = nbCluster)
  tauFinal = matrix(0, nrow = nbArrowLoc, ncol = nbCluster)
  psyAll = matrix(1, nrow = nbEndo, ncol = nbCluster)
  constAll = matrix(0, nrow = nbEndo, ncol = nbCluster)
  
  #----------------------------------------------------------------------------------------
  nbGamma = sum(pathMat[,isExoPur])
  indGamma = rep(0,nbGamma)
  nbBeta = sum(pathMat[,isEndo])
  indBeta = rep(0,nbBeta)
  
  #----------------------------------------------------------------------------------------
  oldTauLoc = tauLoc
  oldPsi = psyAll
  oldPik = pik
  oldConst = constAll
  
  #----------------------------------------------------------------------------------------
  # start loop
  for(iter in 1:itFIM)
  {
    #---------------------------------------------------------------------------
    ## M-step
    #---------------------------------------------------------------------------
    pk = apply(pik, 2, mean)
    
    # path coefficients calculated by regression
    for(k in 1:nbCluster)
    {
      compt_Reg = 1
      tot_Reg = 0
      namesVect = c()
      nbReg = length(which(apply(pathMatLOC,1,sum) != 0)) # =nbEndo
      
      for(m in 1:nbReg)
      {
        
        indEndo = isEndo[m]

        if(indEndo %in% isEndoLoc)
        {
          predEndoLoc = which(pathMatLOC[indEndo,] != 0)
          nbPredLoc = length(predEndoLoc)
          
          tot_Reg = tot_Reg + nbPredLoc
          
          Ym = scoresPLSPM[,indEndo, drop = FALSE] # N x 1
          if(withConst)
          {
            Xm = cbind(1,scoresPLSPM[,predEndoLoc, drop = FALSE]) # N x nbpred
          } else {
            Xm = scoresPLSPM[,predEndoLoc, drop = FALSE] # N x nbpred
          }
          diagPik = diag(pik[,k])

          tau_mk = ginv(t(Xm)%*%diagPik%*%Xm) %*% (t(Xm)%*%diagPik%*%Ym)
          
          tempPik = pik[,k]
          tempPik[which(seg == k)] = 1
          tempPik[which(seg != k)] = 0
          diagPik2 = diag(tempPik)

          tau_mk_final = ginv(t(Xm)%*%diagPik2%*%Xm) %*% (t(Xm)%*%diagPik2%*%Ym)

          if(withConst)
          {
            tauLoc[compt_Reg:tot_Reg,k] = tau_mk[-1]
            tauFinal[compt_Reg:tot_Reg,k] = tau_mk_final[-1]
            constAll[m,k] = tau_mk[1]
          } else {
            tauLoc[compt_Reg:tot_Reg,k] = tau_mk
            tauFinal[compt_Reg:tot_Reg,k] = tau_mk_final
          }
          
          compt_Reg = compt_Reg + nbPredLoc
          
          tempPsy = (t(Ym-Xm%*%tau_mk)%*%diagPik%*%(Ym-Xm%*%tau_mk))/(sum(N*pk[k]))
          if(tempPsy < 0.001)
          {
            psyAll[m,k] = 0.001
          } else {
            psyAll[m,k] = tempPsy
          }
          
        }
        
      }
      rownames(psyAll) = colnames(pathMat)[isEndo]
    }
    
    #----------------------------------------------------------------------------------------
    # likelihood calculation for the specific part
    newLn = computeML(scoresPLSPM, pathMatLOC, tauLoc, constAll, psyAll, nbCluster, pk, seg, 
                       clustLocal, pathMat)
    
    newLnC = computeMLC(scoresPLSPM, pathMatLOC, tauLoc, constAll, psyAll, nbCluster, pk, seg, pik, 
                         clustLocal, pathMat, FALSE)
    
    #---------------------------------------------------------------------------
    ## E-step
    #---------------------------------------------------------------------------
    # compute new proba
    pik = computePik(scoresPLSPM, pathMatLOC, tauLoc, constAll, psyAll, nbCluster, pk, pathMat)
    seg = apply(pik, 1, which.max)
    
    #----------------------------------------------------------------------------------------
    # convergence 
    if(abs(newLnC - oldLnC) < tolFIM || iter == itFIM)
    {
      devLnL[iter,1] = newLn
      devLnL[iter,2] = newLn - oldLn
      devLnLC[iter,1] = newLnC
      devLnLC[iter,2] = newLnC - oldLnC
      
      tauAll = matrix(0, nrow = nbArrow, ncol = nbCluster)
      
      cmpLoc = 1
      cmpAll = 1
      for(indRow in nbExoPur:nbLV)
      {
        for(indCol in 1:nbLV)
        {
          if(pathMat[indRow, indCol] == 1)
          {
            if(pathMatLOC[indRow, indCol] == 1)
            {
              tauAll[cmpAll,] = tauLoc[cmpLoc,]
              cmpLoc = cmpLoc + 1
            } else {
              tauAll[cmpAll,] = rep(tauGlobal[cmpAll],nbCluster)
            }
            cmpAll = cmpAll + 1
          }
        }
      }
      
      pk = apply(pik, 2, mean)
      
      #----------------------------------------------------------------------------------------
      # information criteria
      infoCrit = matrix(NA,11,1)
      rownames(infoCrit) = c("AIC", "AIC3", "AIC4", "BIC", "CAIC", "HQ", "MDL5",
                             "LnL", "EN", "NFI", "NEC")
      
      Ns = (nbCluster - 1) + nbCluster*(nbArrowLoc + nbEndoLoc)
      
      if(withConst)
      {
        Ns = Ns + nbEndoLoc
      }
      
      infoCrit[1,] = (-2*newLn) + 2*Ns # AIC
      infoCrit[2,] = (-2*newLn) + 3*Ns # AIC3
      infoCrit[3,] = (-2*newLn) + 4*Ns # AIC4
      infoCrit[4,] = (-2*newLn) + log(N)*Ns # BIC
      infoCrit[5,] = (-2*newLn) + (log(N)+1)*Ns# CAIC
      infoCrit[6,] = (-2*newLn) + 2*log(log(N))*Ns # HQ
      infoCrit[7,] = (-2*newLn) + 5*log(N)*Ns # MDL5
      infoCrit[8,] = newLn # LnL
      Es = (-1)*sum(pik*log(pik)) # estimated entropy
      infoCrit[9,] = 1 - (Es / (N*log(nbCluster))) # EN
      infoCrit[10,] = (nbCluster*sum(pik*pik)-N) / (N*(nbCluster-1)) # NFI
      infoCrit[11,] = Es / log(nbCluster) # NEC
      
      
      return(list(proba = pik, coeffProb = tauAll, constProb = constAll, psy = psyAll, 
                  segments = seg, segSize = pk, lnL = devLnL, lnLC = devLnLC, 
                  scores = scoresPLSPM, nbIter = iter, critMod = infoCrit)) 
    }
    
    #----------------------------------------------------------------------------------------
    devLnL[iter,1] = newLn
    devLnL[iter,2] = newLn - oldLn
    devLnLC[iter,1] = newLnC
    devLnLC[iter,2] = newLnC - oldLnC
    
    oldLn = newLn
    oldLnC = newLnC
    
    oldPik = pik
    oldPsi = psyAll
    oldTauLoc = tauLoc
    oldConst = constAll
  }
  
}

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# compute log likelihood LnL
computeML <- function(scores, pathMatLOC, tauAll, constAll, psyAll, nbCluster, pk, seg, 
                       clusterLOC = FALSE, pathMatGlo)
{
  # param
  N = nrow(scores)
  nbLV = ncol(pathMatGlo)
  
  isExo = which(colSums(pathMatGlo) != 0)
  nbExo = length(isExo)
  isExoLoc = which(colSums(pathMatLOC) != 0)
  nbExoLoc = length(isExoLoc)
  
  isEndo = which(rowSums(pathMatGlo) != 0)
  nbEndo = length(isEndo)
  isEndoLoc = which(rowSums(pathMatLOC) != 0)
  nbEndoLoc = length(isEndoLoc)
  
  isExoPur = which(rowSums(pathMatGlo) == 0)
  nbExoPur = nbLV - nbEndo
  isExoPurLoc = which(rowSums(pathMatLOC) == 0)
  nbExoPurLoc = nbLV - nbEndoLoc
  
  lnL = 0
  
  for(i in 1:N)
  {
    sumF = 0
    
    for(k in 1:nbCluster)
    {
      
      tempCoeff = pathMatLOC
      cmp = 1
      for(indRow in 1:nbLV)
      {
        for(indCol in 1:nbLV)
        {
          if(tempCoeff[indRow, indCol] == 1)
          {
            tempCoeff[indRow, indCol] = tauAll[cmp,k]
            cmp = cmp + 1
          }
        }
      }
      
      tempCoeff = tempCoeff[isEndo,, drop = FALSE]
      gammaK = tempCoeff[,-isEndo, drop = FALSE]
      betaK = tempCoeff[,(nbExoPur+1):nbLV, drop = FALSE]
      constK = constAll[,k,drop = FALSE]
      
      # compute fik
      # constant 
      const = 1 / ( (2*pi)^(nbEndoLoc/2) * sqrt(prod(psyAll[,k])) )
      # expo part
      temp = (diag(nbEndo)-betaK)%*%t(scores[i,(nbExoPur+1):nbLV, drop = FALSE]) +
        (-gammaK)%*%t(scores[i,-isEndo, drop = FALSE]) + (-constK)
      
      if(nbEndo == 1){expPart = t(temp)%*%(1/psyAll[,k])%*%temp}
      else{expPart = t(temp)%*%diag(as.vector(1/psyAll[,k]))%*%temp}
      
      # all
      fik = exp((-0.5)*expPart)
      sumF = sumF + pk[k] * const * fik
    }
    
    lnL = lnL + log(sumF)
    
  }
  
  return(as.numeric(lnL))
}

# -------------------------------------------------------------------------------
# compute log likekihood complete LnLC
computeMLC <- function(scores, pathMatLOC, tauAll, constAll, psyAll, nbCluster, pk, seg, pik, 
                        clusterLOC = FALSE, pathMatGlo, display = FALSE)
{
  # param
  N = nrow(scores)
  nbLV = ncol(pathMatGlo)
  isExo = which(colSums(pathMatGlo) != 0)
  nbExo = length(isExo)
  isEndo = which(rowSums(pathMatGlo) != 0)
  nbEndo = length(isEndo)
  isEndoLoc = which(rowSums(pathMatLOC) != 0)
  nbEndoLoc = length(isEndoLoc)
  isExoPur = which(rowSums(pathMatGlo) == 0)
  nbExoPur = nbLV - nbEndo
  
  lnLC = 0
  
  for(i in 1:N)
  {
    sumF = 0
    
    for(k in 1:nbCluster)
    {
      tempCoeff = pathMatLOC
      cmp = 1
      for(indRow in 1:nbLV)
      {
        for(indCol in 1:nbLV)
        {
          if(tempCoeff[indRow, indCol] == 1)
          {
            tempCoeff[indRow, indCol] = tauAll[cmp,k]
            cmp = cmp + 1
          }
        }
      }
      
      tempCoeff = tempCoeff[isEndo,, drop = FALSE]
      gammaK = tempCoeff[,-isEndo, drop = FALSE]
      betaK = tempCoeff[,(nbExoPur+1):nbLV, drop = FALSE]
      constK = constAll[,k,drop = FALSE]
      
      # compute fik
      # log const
      const = 1 / ( (2*pi)^(nbEndoLoc/2) * sqrt(prod(psyAll[,k])) )
      # constLn = log(const)
      # constLn = (-1)*((nbEndoLoc/2)*log(2*pi)+0.5*log(prod(psyAll[,k])))
      # expo part
      temp = (diag(nbEndo)-betaK)%*%t(scores[i,(nbExoPur+1):nbLV, drop = FALSE]) +
        (-gammaK)%*%t(scores[i,-isEndo, drop = FALSE]) + (-constK)
      
      if(nbEndo == 1){expPart = t(temp)%*%(1/psyAll[,k])%*%temp}
      else{expPart = t(temp)%*%diag(as.vector(1/psyAll[,k]))%*%temp}
      
      # all
      fik = const * exp((-0.5)*expPart)
      zik = pik[i,k]
      
      if(fik == 0)
      {
        sumF = sumF + 0
      } else {
        sumF = sumF + zik * (log(fik)+log(pk[k]))
      }
      
    }
    
    lnLC = lnLC + sumF
  }
  
  return(lnLC)
  
}

# -------------------------------------------------------------------------------
# compute proba
computePik <- function(scores, pathMatLOC, tauAll, constAll, psyAll, nbCluster, pk, pathMatGlo)
{
  # param / init
  N = nrow(scores)
  nbLV = nrow(pathMatGlo)
  
  isExo = which(colSums(pathMatGlo) != 0)
  nbExo = length(isExo)
  isEndo = which(rowSums(pathMatGlo) != 0)
  nbEndo = length(isEndo)
  isEndoLoc = which(rowSums(pathMatLOC) != 0)
  nbEndoLoc = length(isEndoLoc)
  nbExoPur = nbLV - nbEndo
  nbEndoPur = nbLV - nbExo
  
  pik = matrix(0, nrow = N, ncol= nbCluster)
  sumK = rep(0, N)
  
  for(i in 1:N)
  {
    for(k in 1:nbCluster)
    {
      tempCoeff = pathMatLOC
      cmp = 1
      for(indRow in 1:nbLV)
      {
        for(indCol in 1:nbLV)
        {
          if(tempCoeff[indRow, indCol] == 1)
          {
            tempCoeff[indRow, indCol] = tauAll[cmp,k]
            cmp = cmp + 1
          }
        }
      }
      
      tempCoeff = tempCoeff[isEndo,, drop = FALSE]
      gammaK = tempCoeff[,-isEndo, drop = FALSE]
      betaK = tempCoeff[,(nbExoPur+1):nbLV, drop = FALSE]
      constK = constAll[,k,drop=FALSE]
      
      # compute fik
      # constant
      const = 1 / ( (2*pi)^(nbEndoLoc/2) * sqrt(prod(psyAll[,k])) )
      # expo part
      temp = (diag(nbEndo)-betaK)%*%t(scores[i,(nbExoPur+1):nbLV, drop = FALSE]) +
        (-gammaK)%*%t(scores[i,-isEndo, drop = FALSE]) + (-constK)
      
      if(nbEndo == 1){expPart = t(temp)%*%(1/psyAll[,k])%*%temp}
      else{expPart = t(temp)%*%diag(as.vector(1/psyAll[,k]))%*%temp}
      
      # all
      fik = exp((-0.5)*expPart)
      tempPik = pk[k] * const * fik
      pik[i,k] = tempPik
      sumK[i] = sumK[i] + tempPik
      
    }
    
    
  }
  
  pik = sweep(pik, 1, sumK, "/")
  
  return(pik)
}

