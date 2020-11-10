# Use requires attribution to Author Giacomo Bignardi giacomo.bignardi@mrc-cbu.cam.ac.uk

# Required Packages -----------------------------------------------------------------------------------------

if (!require("caret")) install.packages("caret")
if (!require("patchwork")) install.packages("patchwork")
if (!require("jtools")) install.packages("jtools")
if (!require("MCMCpack")) install.packages("MCMCpack")
if (!require("dataPreparation")) install.packages("dataPreparation")
if (!require("matrixStats")) install.packages("matrixStats")
if (!require("Rfast")) install.packages("Rfast")

library(MCMCpack)
library(jtools)
library(patchwork)
library(caret)


# CCA with training/testing datasets ------------------------------------------------------------------------
    # Simplifies the use of the gb_CCA function by automatically estimating confidence intervals & p-values

# X_FIT = df0_imputed_group1[,envvar]
# Y_FIT = df0_imputed_group1[,outcomes]
# X_PRED = df0_imputed_group2[,envvar]
# Y_PRED = df0_imputed_group2[,outcomes]
# ncomp=NULL
# alpha = .05


gb_CCA_splithalf = function(X_FIT,Y_FIT,X_PRED,Y_PRED, 
                            ncomp=NULL, alpha = 0.05){
  
  model_results = gb_CCA(X_FIT=X_FIT,Y_FIT=Y_FIT,X_PRED=X_PRED,Y_PRED=Y_PRED, 
                         ncomp=ncomp, 
                         ProcrustX = NULL, ProcrustY = NULL,
                         SafetyChecks=TRUE)
  
  # Estimate Confidence Intervals
  
  samplesize = nrow(X_PRED)
  predicted_cc = model_results$cc_pred
  confint_cc = sapply(predicted_cc, function(x) CorrelationCIEstimator(x, samplesize, alpha=alpha)$CI)
    rownames(confint_cc) = c("LB","UB") # 'lower bound' and 'upper bound' of confidence interval 
    colnames(confint_cc) = paste0("cc",1:ncol(confint_cc))
  pvalue_cc = sapply(predicted_cc, function(x) CorrelationCIEstimator(x, samplesize, alpha=alpha)$p)
  combined_cc = rbind.data.frame(predicted_cc, confint_cc, pvalue_cc)
    rownames(combined_cc) = c("cc", "LB", "UB", "p")
  
  return(list(
      model_results = model_results,
      predicted_cc  = predicted_cc,
      confint_cc    = confint_cc,
      pvalue_cc     = pvalue_cc,
      combined_cc   = combined_cc
  ))
  
} 



# Bootstrap CCA Code  ---------------------------------------------------------------------------------------
# This function runs the gb_CCA Canonical Correlation Analysis function multiple times to assess variability in the CCA loadings and canonical correlations
    # Because bootstrap resampling can change the order of canonical variates that are extracted, or sign flipping can occur 
    # in some cases (i.e. a very similar latent variable is extracted but on some occasions the loadings are mostly positive or negative), we rotate the loadings
    # in each bootstrap resample to map onto the loadings generated from the full, raw input datsets. 

# X_FIT - matrix including 'predictor' variables
# Y_FIT - matrix inluding 'outcome' variables
# ncomp - number of components to extract (must be lower than number of variables in X or Y matrix )
# Nboot - number of bootstrap resamples to perform (ideally 1000-3000 as minimum)

gb_CCA_boot = function(X_FIT,Y_FIT, ncomp=10, Nboot=30,ProcrustX = NULL, ProcrustY = NULL){
  # browser()
  pb <- utils::txtProgressBar(min = 1, max = Nboot, style = 3)
  
  #Run model on full dataset 
  CCA_OriginalData = gb_CCA(X_FIT=X_FIT, Y_FIT=Y_FIT, X_PRED=NULL, Y_PRED=NULL, 
                            ProcrustX = ProcrustX, ProcrustY = ProcrustY,
                            ncomp=ncomp)

  #Lists to store bootstrap data in
  YLoadings_ROTATED = list()  #Store rotated loadings from bootstrap resamples in here 
  XLoadings_ROTATED = list()
  cc = list()                #Store canonical correlations from bootstrap resamples in hre 
  
  for(i in 1:Nboot){
    # Randomly bootstrap resample data and fit CCA model. 
    BootResample = sample(nrow(X_FIT), replace = TRUE)
  
    
    #Run CCA on bootstrap resampled data, rotating loadings using procrustes to map  onto the original model
    CCA_BootData = gb_CCA(X_FIT=X_FIT[BootResample,], Y_FIT=Y_FIT[BootResample,], 
                          ProcrustX = CCA_OriginalData$XLoadings, ProcrustY = CCA_OriginalData$YLoadings, 
                          X_PRED=NULL, Y_PRED=NULL, 
                          ncomp=ncomp)
    
    XLoadings_ROTATED[[i]] = CCA_BootData$XLoadings
    YLoadings_ROTATED[[i]] = CCA_BootData$YLoadings
    cc[[i]]  = CCA_BootData$cc_pred
    
  setTxtProgressBar(pb, value=i)
    
  }
#Estimate quantiles from boostrap distribution
  
  XLoadings_Quantiles = apply(base::simplify2array(XLoadings_ROTATED), 1:2, quantile, prob = c(0.025, .5, 0.975))
  XLoadings_Quantiles = lapply(1:ncomp, function(i) XLoadings_Quantiles[,,i])
  for(i in 1:ncomp){
    colnames(XLoadings_Quantiles[[i]]) = colnames(X_FIT)
    XLoadings_Quantiles[[i]] = data.frame(t(XLoadings_Quantiles[[i]]))
    XLoadings_Quantiles[[i]]$original = CCA_OriginalData$XLoadings[,i]
  }

  YLoadings_Quantiles = apply(base::simplify2array(YLoadings_ROTATED), 1:2, quantile, prob = c(0.025, .5, 0.975))
  YLoadings_Quantiles = lapply(1:ncomp, function(i) YLoadings_Quantiles[,,i])
  for(i in 1:ncomp){
    colnames(YLoadings_Quantiles[[i]]) = colnames(Y_FIT)
    YLoadings_Quantiles[[i]] = data.frame(t(YLoadings_Quantiles[[i]]))
    YLoadings_Quantiles[[i]]$original = CCA_OriginalData$YLoadings[,i]
  }
  
  cc_quantiles = apply(base::simplify2array(cc), 1, quantile, prob = c(0.025, .5, 0.975))
  
  
# Ouput Data 
  
  out = list(
    XLoadings_Quantiles=XLoadings_Quantiles,
    YLoadings_Quantiles=YLoadings_Quantiles,
    cc_quantiles,
    cc=cc
  )
  
  return(out)
  
}



# Cross-Validtation Bootstrap CCA code  ---------------------------------------------------------------------------------------
      # This function runs the gb_CCA Canonical Correlation Analysis function multiple times to assess variability and bias in canonical correlations
        # On each iteration, N-fold (default is 10 fold) cross-validation is used to generated predicted canonical variates for the complete sample. Following this,
        # the predicted variates are bootstrap resampled and canonical correlations are estimated from them. 
      # Because bootstrap resampling can change the order of canonical variates that are extracted, or sign flipping can occur 
        # in some cases (i.e. a very similar latent variable is extracted but on some occasions the loadings are mostly positive or negative), we rotate the loadings
        # in each during cross-validation to map onto the loadings generated from the full, raw input datsets. 

# NOTE; 
# NBoot - Number of times to repeat the cross-validation + bootstrap procedure. In each iteration a brand new set of cross-validated predicted variates are generated and
    # ... bootstrap resampled from. 

gb_CCA_CVboot = function(X_FIT,Y_FIT, ncomp=10, Nboot=30, Nfolds=10,ProcrustX = NULL, ProcrustY = NULL, UseProgressBar=TRUE, UseProcrustes=TRUE){
  # browser()
  if (UseProgressBar){
    pb <- utils::txtProgressBar(min = 0, max = Nboot, style = 3)
  }
  
  CCA_OriginalData = gb_CCA(X_FIT=X_FIT, Y_FIT=Y_FIT, X_PRED=NULL, Y_PRED=NULL, ncomp=ncomp, ProcrustX = ProcrustX, ProcrustY = ProcrustY)
  
  if (UseProcrustes==FALSE){
    CCA_OriginalData$XLoadings = NULL # By setting this to NULL, the gb_CCA function called below will not use procrustes rotations
    CCA_OriginalData$YLoadings = NULL
  }
  
  cc_CVboot = list()
  cc_CV = list()
  R2_matrix = list()
  
  for(b in 1:Nboot){
    #Divide data into folds... 
    df_folds = caret::createFolds(1:nrow(X_FIT), k=Nfolds)
    # df_folds = cut(sample(nrow(X_FIT)), breaks=Nfolds, labels=FALSE)
    
    Variate_predictions = list()
    for (f in 1:Nfolds){
      
      # Fit_Index = (df_folds!=f) #Rows to select to fit data to 
      # Pred_Index = (df_folds==f) #Rows to select to make predictions for
      # 
      Fit_Index = .Internal(unlist(df_folds[-f], FALSE, FALSE)) #Row numbers - Training Data
      Pred_Index = .Internal(unlist(df_folds[f], FALSE, FALSE)) #Row numbers - Hold-Out Data

      #Estimate CCA in trainning dataset and generate list of predictions for hold-out data - and append to Variate_predictions list
      Variate_predictions =  c(Variate_predictions,
                               list(gb_CCA(X_FIT  = X_FIT[Fit_Index,],  Y_FIT  = Y_FIT[Fit_Index,], 
                                           X_PRED = X_FIT[Pred_Index,], Y_PRED = Y_FIT[Pred_Index,], ncomp=ncomp,
                                           ProcrustX = CCA_OriginalData$XLoadings, ProcrustY = CCA_OriginalData$YLoadings)$Variates)
                              )
      
    }
    
    #Put cross-validated predictions back in original order and in a data frame 
    Variates_CrossValidated = data.table::rbindlist(Variate_predictions)
    Variates_CrossValidated = Variates_CrossValidated[order(.Internal(unlist(df_folds, FALSE, FALSE))),]
    
    #Bootstrap cross-validation predictions
    boot_i = .Internal(sample(nrow(Variates_CrossValidated),nrow(Variates_CrossValidated), TRUE, NULL))
    Variates_CrossValidated_b = as.matrix(Variates_CrossValidated[boot_i,])
    
    #Estimate canonical correlations
    cc_CVboot[[b]] = Rfast::corpairs(as.matrix(Variates_CrossValidated_b[,1:ncomp]),as.matrix(Variates_CrossValidated_b[,(ncomp+1):(2*ncomp)]))
    # cc_CVboot = c(cc_CVboot , list(Rfast::corpairs(as.matrix(Variates_CrossValidated_b[,1:ncomp]),as.matrix(Variates_CrossValidated_b[,(ncomp+1):(2*ncomp)])))) #Not necessary 
    
    #Estimate R2 for all outcome variables (with boot)
    Variates_CV_Scaled = as.matrix(dataPreparation::fastScale(Variates_CrossValidated_b, verbose = FALSE))
    Y_FIT_Scaled =       as.matrix(dataPreparation::fastScale(Y_FIT[boot_i,],verbose=FALSE))

    R2_matrix[[b]] =  
        sapply(1:ncomp, function(ncomp_i)
          sapply(1:ncol(Y_FIT), function(y_i) 
            R2quickcalc(X=Variates_CV_Scaled[,1:ncomp_i],Y=Y_FIT_Scaled[,y_i])
          )
        )
    
    #Standard Cross-Validation canonical correlation
    cc_CV[[b]] = Rfast::corpairs(as.matrix(Variates_CrossValidated[,1:ncomp]),as.matrix(Variates_CrossValidated[,(ncomp+1):(2*ncomp)]))
    
    
    if (UseProgressBar){
      setTxtProgressBar(pb, value=b)
    }
  }
  
  
  # R2 Means 
  R2_matrix = apply(base::simplify2array(R2_matrix), 1:2, mean)
  #Prettify output 
  R2_matrix = data.frame(R2_matrix)
  rownames(R2_matrix) = colnames(Y_FIT)
  colnames(R2_matrix) = paste0("NVariates_",1:ncomp)

  # Cross-Validation Quantiles
  cc_CV_quantiles = do.call("rbind.data.frame", cc_CV)
  colnames(cc_CV_quantiles) = paste0("cc",1:ncomp)
  cc_CV_quantiles = apply(cc_CV_quantiles, 2, function(x) quantile(x, probs = c(0.025,.5,.975), type=6))
  # cc_CV_pval = apply(cc_CV_quantiles, 2, function(x) boot_pval(x))
  
  # Cross-validation + Bootstrap Quantiles
  cc_CVBoot_quantiles = do.call("rbind.data.frame", cc_CVboot)
  colnames(cc_CVBoot_quantiles) = paste0("cc",1:ncomp)
  cc_CVBoot_quantiles2 = apply(cc_CVBoot_quantiles, 2, function(x) quantile(x, probs = c(0.025,.5,.975), type=6))
  cc_CVBoot_pval = apply(cc_CVBoot_quantiles, 2, function(x) boot_pval(x))
  
  
  return(list(
    R2_matrix = R2_matrix,
    CrossValidationQuantiles = cc_CV_quantiles,
    CrossValidationBootstrapQuantiles = cc_CVBoot_quantiles2,
    CrossValidationBootstrapPvalues = cc_CVBoot_pval
    
  ))
  
  
}

# Cross-Validtation Bootstrap CCA code  ---------------------------------------------------------------------------------------
      # This function runs the gb_CCA Canonical Correlation Analysis function multiple times to assess variability and bias in canonical correlations
        # On each iteration, N-fold (default is 10 fold) cross-validation is used to generated predicted canonical variates for the complete sample. Following this,
        # the predicted variates are bootstrap resampled and canonical correlations are estimated from them. 
      # Because bootstrap resampling can change the order of canonical variates that are extracted, or sign flipping can occur 
        # in some cases (i.e. a very similar latent variable is extracted but on some occasions the loadings are mostly positive or negative), we rotate the loadings
        # in each during cross-validation to map onto the loadings generated from the full, raw input datsets. 

# NOTE; 
# NBoot - Number of times to repeat the cross-validation + bootstrap procedure. In each iteration a brand new set of cross-validated predicted variates are generated and
    # ... bootstrap resampled from. 

gb_CCA_CVboot = function(X_FIT,Y_FIT, ncomp=10, Nboot=30, Nfolds=10,ProcrustX = NULL, ProcrustY = NULL, UseProgressBar=TRUE, UseProcrustes=TRUE){
  # browser()
  if (UseProgressBar){
    pb <- utils::txtProgressBar(min = 0, max = Nboot, style = 3)
  }
  
  CCA_OriginalData = gb_CCA(X_FIT=X_FIT, Y_FIT=Y_FIT, X_PRED=NULL, Y_PRED=NULL, ncomp=ncomp, ProcrustX = ProcrustX, ProcrustY = ProcrustY)
  
  if (UseProcrustes==FALSE){
    CCA_OriginalData$XLoadings = NULL # By setting this to NULL, the gb_CCA function called below will not use procrustes rotations
    CCA_OriginalData$YLoadings = NULL
  }
  
  cc_CVboot = list()
  cc_CV = list()
  R2_matrix = list()
  
  for(b in 1:Nboot){
    #Divide data into folds... 
    df_folds = caret::createFolds(1:nrow(X_FIT), k=Nfolds)
    # df_folds = cut(sample(nrow(X_FIT)), breaks=Nfolds, labels=FALSE)
    
    Variate_predictions = list()
    for (f in 1:Nfolds){
      
      # Fit_Index = (df_folds!=f) #Rows to select to fit data to 
      # Pred_Index = (df_folds==f) #Rows to select to make predictions for
      # 
      Fit_Index = .Internal(unlist(df_folds[-f], FALSE, FALSE)) #Row numbers - Training Data
      Pred_Index = .Internal(unlist(df_folds[f], FALSE, FALSE)) #Row numbers - Hold-Out Data

      #Estimate CCA in trainning dataset and generate list of predictions for hold-out data - and append to Variate_predictions list
      Variate_predictions =  c(Variate_predictions,
                               list(gb_CCA(X_FIT  = X_FIT[Fit_Index,],  Y_FIT  = Y_FIT[Fit_Index,], 
                                           X_PRED = X_FIT[Pred_Index,], Y_PRED = Y_FIT[Pred_Index,], ncomp=ncomp,
                                           ProcrustX = CCA_OriginalData$XLoadings, ProcrustY = CCA_OriginalData$YLoadings)$Variates)
                              )
      
    }
    
    #Put cross-validated predictions back in original order and in a data frame 
    Variates_CrossValidated = data.table::rbindlist(Variate_predictions)
    Variates_CrossValidated = Variates_CrossValidated[order(.Internal(unlist(df_folds, FALSE, FALSE))),]
    
    #Bootstrap cross-validation predictions
    boot_i = .Internal(sample(nrow(Variates_CrossValidated),nrow(Variates_CrossValidated), TRUE, NULL))
    Variates_CrossValidated_b = as.matrix(Variates_CrossValidated[boot_i,])
    
    #Estimate canonical correlations
    cc_CVboot[[b]] = Rfast::corpairs(as.matrix(Variates_CrossValidated_b[,1:ncomp]),as.matrix(Variates_CrossValidated_b[,(ncomp+1):(2*ncomp)]))
    # cc_CVboot = c(cc_CVboot , list(Rfast::corpairs(as.matrix(Variates_CrossValidated_b[,1:ncomp]),as.matrix(Variates_CrossValidated_b[,(ncomp+1):(2*ncomp)])))) #Not necessary 
    
    #Estimate R2 for all outcome variables (with boot)
    Variates_CV_Scaled = as.matrix(dataPreparation::fastScale(Variates_CrossValidated_b, verbose = FALSE))
    Y_FIT_Scaled =       as.matrix(dataPreparation::fastScale(Y_FIT[boot_i,],verbose=FALSE))

    R2_matrix[[b]] =  
        sapply(1:ncomp, function(ncomp_i)
          sapply(1:ncol(Y_FIT), function(y_i) 
            R2quickcalc(X=Variates_CV_Scaled[,1:ncomp_i],Y=Y_FIT_Scaled[,y_i])
          )
        )
    
    #Standard Cross-Validation canonical correlation
    cc_CV[[b]] = Rfast::corpairs(as.matrix(Variates_CrossValidated[,1:ncomp]),as.matrix(Variates_CrossValidated[,(ncomp+1):(2*ncomp)]))
    
    
    if (UseProgressBar){
      setTxtProgressBar(pb, value=b)
    }
  }
  
  
  # R2 Means 
  R2_matrix = apply(base::simplify2array(R2_matrix), 1:2, mean)
  #Prettify output 
  R2_matrix = data.frame(R2_matrix)
  rownames(R2_matrix) = colnames(Y_FIT)
  colnames(R2_matrix) = paste0("NVariates_",1:ncomp)

  # Cross-Validation Quantiles
  cc_CV_quantiles = do.call("rbind.data.frame", cc_CV)
  colnames(cc_CV_quantiles) = paste0("cc",1:ncomp)
  cc_CV_quantiles = apply(cc_CV_quantiles, 2, function(x) quantile(x, probs = c(0.025,.5,.975), type=6))
  # cc_CV_pval = apply(cc_CV_quantiles, 2, function(x) boot_pval(x))
  
  # Cross-validation + Bootstrap Quantiles
  cc_CVBoot_quantiles = do.call("rbind.data.frame", cc_CVboot)
  colnames(cc_CVBoot_quantiles) = paste0("cc",1:ncomp)
  cc_CVBoot_quantiles2 = apply(cc_CVBoot_quantiles, 2, function(x) quantile(x, probs = c(0.025,.5,.975), type=6))
  cc_CVBoot_pval = apply(cc_CVBoot_quantiles, 2, function(x) boot_pval(x))
  
  
  return(list(
    R2_matrix = R2_matrix,
    CrossValidationQuantiles = cc_CV_quantiles,
    CrossValidationBootstrapQuantiles = cc_CVBoot_quantiles2,
    CrossValidationBootstrapPvalues = cc_CVBoot_pval
    
  ))
  
  
}

# --------------------------------------------------------------------------------------------------------------------------
# Internal Functions  ------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------

# Standard CCA function with additional functionality  ----------------------------------------------------------------------
          # The function can fit a CCA model in one dataset, and use the loadings (coefficients) generated from this to predict CCA variate (latent variable) scores in a new dataset.
          # In addition, the function can rotate the generated CCA loadings to a target matrix of X and Y loadings (ProcrustX and ProcrustY)
          # Note: assumes NO missing data

gb_CCA = function(X_FIT,Y_FIT,X_PRED=NULL,Y_PRED=NULL, 
                  ncomp=10, 
                  ProcrustX = NULL, ProcrustY = NULL,
                  SafetyChecks=FALSE){
  # browser()
  #Check some basic things 
  if (SafetyChecks){
    if (nrow(X_FIT)!=nrow(Y_FIT)) stop("nrow of X_FIT and Y_FIT do not match")
    if (!is.null(ProcrustX)){
      if (ncol(ProcrustX)!=ncomp) stop("ProcrustX should have same number of columns to ncomp")
    }
    if (!is.null(ProcrustY)){
      if (ncol(ProcrustY)!=ncomp) stop("ProcrustY should have same number of columns to ncomp")
    }
    if (!is.null(ncomp)){
      if ((ncomp>ncol(X_FIT)) | (ncomp>ncol(Y_FIT)) ) stop("ncomp should be equal to or less than smallest number of variables in X_FIT or X_FIT")
    }
  }
  
  #Set ncomp to the minimum of ncol
  ncomp = min(c(ncomp,ncol(X_FIT),ncol(Y_FIT)))
  
  #Scale FIT datasets
  X = as.matrix(dataPreparation::fastScale(X_FIT, verbose=FALSE))
  Y = as.matrix(dataPreparation::fastScale(Y_FIT, verbose=FALSE))
  
  #Scale the PRED matrix using the FIT m and sd  - X variables
  if (is.null(X_PRED)){
    X_PRED = X
  } else {
    # Scale the PRED matrix by the FIT matrix
    X_FIT_mat = as.matrix(X_FIT)
    X_Mean = matrixStats::colMeans2(X_FIT_mat)
    X_SD = matrixStats::colSds(X_FIT_mat)
    X_PRED = t((t(X_PRED)-X_Mean)/X_SD)
  }
  
  #Scale the PRED matrix using the FIT m and sd  - Y variables
  if (is.null(Y_PRED)){
    Y_PRED = Y
  } else {
    # Scale the PRED matrix by the FIT matrix
    Y_FIT_mat = as.matrix(Y_FIT)
    Y_Mean = matrixStats::colMeans2(Y_FIT_mat)
    Y_SD = matrixStats::colSds(Y_FIT_mat)
    Y_PRED = t((t(Y_PRED)-Y_Mean)/Y_SD)
  }
  
  # Estimate Eigenvalues and Eigenvectors... 
  Sxx  = crossprod(X,X)
  Syx  = crossprod(Y,X)
  Syy  = crossprod(Y,Y)
  Sxy  = crossprod(X,Y)
  
  #Slightly slower to use below code! 
  # Sxx  = t(X) %*% X
  # Syx  = t(Y) %*% X
  # Syy  = t(Y) %*% Y
  # Sxy  = t(X) %*% Y
  
  V = solve(Syy) %*% Syx %*% solve(Sxx) %*% Sxy
  W = solve(Sxx) %*% Sxy %*% solve(Syy) %*% Syx
  
  Eig_V = eigen(V)
  Eig_W = eigen(W)
  
  cc = sqrt(Eig_V$values)[1:ncomp] # Canonical correlations from FITTED MATRIX
  
  # Get loadings 
  Yloadings = apply(as.matrix(Eig_V$vectors[,1:ncomp]),2,as.numeric)
  Xloadings = apply(as.matrix(Eig_W$vectors[,1:ncomp]),2,as.numeric)
  
  # Optional Rotation of Loadings 
  
  if (!is.null(ProcrustX)){
    Xloadings = as.matrix(MCMCpack::procrustes(Xloadings  , ProcrustX)$X.new)
    
  }
  if (!is.null(ProcrustY)){
    Yloadings = as.matrix(MCMCpack::procrustes(Yloadings  , ProcrustY)$X.new)
  }
  
  # Estimate Canonical Variates
  
  Variates_Y =  Y_PRED %*% Yloadings
  colnames(Variates_Y) = paste0("Y", 1:ncomp)
  
  Variates_X =  X_PRED %*% Xloadings
  colnames(Variates_X) = paste0("X", 1:ncomp)
  
  Variates = cbind.data.frame(Variates_X, Variates_Y)
  
  cc_pred = Rfast::corpairs(Variates_X, Variates_Y)
  
  # Return Output 
  
  out = list(
    #Loadings 
    YLoadings = Yloadings,
    XLoadings = Xloadings,
    
    #Variates 
    Variates = Variates,
    Variates_X = Variates_X,
    Variates_Y = Variates_Y,
    
    #Canonical Correlations 
    cc_fit = cc,                            #Canonical correlations estimated from the input data (fitted canonical correlations)
    cc_pred = cc_pred                       #Canonical correlations estimated from the prediction data (predicted canonical correlations)
  )
  
  
  return(out)
  
}


# Function for VERY QUICK linear modelling  ---------------------------------------------------------------------------------------
          # returns an R-squared as output
          # WARNING - requires input to be scaled matrices!!
R2quickcalc = function(X,Y){
  # browser()
  out1= Rfast::lmfit(x=X,y=Y)
  Y = as.numeric(Y)
  return(1- sum(out1$residuals^2)/sum((Y-mean(Y))^2))
}

#Slower version, but does not require full rank matrices?
R2quickcalc_v2 = function(X,Y){
  # browser()
  out1= .lm.fit(x=X,y=Y)
  Y = as.numeric(Y)
  return(1- sum(out1$residuals^2)/sum((Y-mean(Y))^2))
}

# Estimate p=value from Bootstrap distribution ---------------------------------------------------------------------------------
boot_pval = function(x, null_val=0){
  x = na.omit(x)
  perc = length(which(x<null_val))/length(x)
  p_val = 1-abs(.50-perc)*2
  return(p_val)
}

# Function that is purely for simulation purposes ---------------------------------------------------------------------------------
# Automatically splis data in half - fits data in train half, and finds canonical corrlations in test half

gb_CCA_SplitHalfSim = function(X_FIT,Y_FIT, ncomp=10, ProcrustX = NULL, ProcrustY = NULL){
  
  df_folds = caret::createFolds(1:nrow(X_FIT), k=2)  #Split Data into halfs
  
  CCA_results = gb_CCA(X_FIT=X_FIT[df_folds[[1]],], Y_FIT=Y_FIT[df_folds[[1]],], X_PRED=X_FIT[df_folds[[2]],], Y_PRED=Y_FIT[df_folds[[2]],], ncomp=ncomp, ProcrustX = ProcrustX, ProcrustY = ProcrustY)
  
  CanonicalCorrelations = CCA_results$cc_pred # Estimated Canonical Correlations in hold-out dataset 
  
  SampleSize = nrow(X_FIT[df_folds[[2]],])
  
  CC_ConfidenceIntervals = t(sapply(CanonicalCorrelations, function(x) CorrelationCIEstimator(x, n=SampleSize, alpha = .05)$CI))
    rownames(CC_ConfidenceIntervals) = paste0("cc",1:nrow(CC_ConfidenceIntervals))
    colnames(CC_ConfidenceIntervals) = c("CI_LB","CI_UB")
  CC_pvalues =  t(sapply(CanonicalCorrelations, function(x) CorrelationCIEstimator(x, n=SampleSize, alpha = .05)$p))

  return(list(
    CanonicalCorrelations=CanonicalCorrelations,
    CC_ConfidenceIntervals=CC_ConfidenceIntervals,
    CC_pvalues=CC_pvalues
  ))
   
}


# Pearson Correlation Confidence Interval Calculator  ---------------------------------------------------------------------------------

CorrelationCIEstimator = function(r, n, alpha=0.05){
  Fr = atanh(r)                 # Fisher Z Transform
  SE = 1/((n-3)^.5)             # Standard Error
  CI = qnorm(c(alpha/2,1-alpha/2), mean=Fr, sd=SE)
  CI = tanh(CI)
  # p  = (1-pnorm(abs(Fr), mean=0, sd=SE))*2    # Fisher Z P value
  t = r*sqrt((n-2)/(1-r^2))       # P-value estimated from t-distribution
  p  = (1-pt(abs(t), df=n-2))*2 
  return(list(CI=CI,p=p))
}


