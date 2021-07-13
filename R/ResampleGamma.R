
###################################################################################
########## Function 4: Resampling function to get bootstrap variance for binary Y
##################################################################################

#' Resample for bootstrap variance for binary Y
#' 
#' Resampling function to get bootstrap variance for binary Y.
#' 
#' @param data synthetic data
#' @param indices indices ??
#' 
#' @return numeric vector of regression coefficients
#' 
#' @references Reference: Gu, T., Taylor, J.M.G. and Mukherjee, B. (2021) Regression 
#' inference for multiple populations by integrating summary-level data using stacked 
#' imputations <arXiv:http://arxiv.org/abs/2106.06835>. 
#' 
#' @export
Resample.gamma.binaryY <- function(data, indices){
  d = data[indices,]
  datan = d
  
  ################ Step 1: convert the external model information into the synthetic data
  data.combined = Create.Synthetic(nrep=nrep, 
                                   datan=datan, 
                                   Y='Y', 
                                   XB = c('X1','X2','X3','X4','B1','B2'), 
                                   Ytype='binary',
                                   parametric = c('Yes','No'),
                                   betaHatExt_list=betaHatExt_list, 
                                   sigmaHatExt_list=NULL)
  
  #### For the external risk calculator 2 (k=2), we need to manually create the synthetic data, and then row-bind with the dataset data.combined
  X.2 = c('X1','X2','X3','X4') #predictors used in the external risk calculator 2
  B.2 = c('B1','B2')  #unobserved variables in the external risk calculator 2
  n = dim(datan)[1]   #total size of the internal data
  m = n*nrep          #total number of synthetic data for S=2
  data.2 = data.frame(matrix(ncol = 8, nrow = m))
  colnames(data.2) = c('Y', 'X1', 'X2', 'X3', 'X4', 'B1', 'B2', 'S')
  data.2[,B.2] = NA  #unobserved variables remain missing 
  data.2[,X.2] = do.call("rbind", replicate(nrep, datan[,X.2], simplify = FALSE)) #replicate the observed X.2 nrep times
  data.2$Y = rbinom(m, 1, predict(fit.rf, newdata=data.2[,X.2])) #generate synthetic Y values via the external risk calculator 2 (random forest model)
  data.2$Y[is.na(data.2$Y)] = 1
  data.2$S = 2
  data.2 = cbind("Int" = 1, data.2)
  
  ##### combine all synthetic data
  data.combined = rbind(data.combined, data.2)
  
  ################ Step 2: Multiple imputation (impute missingness ignoring the outcome Y)
  pred_matrix = mice::make.predictorMatrix(data.combined)
  pred_matrix[c('Int',"Y",'X1','X2','X3','S'),] = 0
  pred_matrix[,c('Int','Y','S')] = 0
  imp_method = mice::make.method(data.combined)
  imp_method[c('X4','B1','B2')] = c('norm','norm','logreg') #choose your desired imputation method for each missingness
  data.combined$B2 = factor(data.combined$B2)
  imputes = mice::mice(data.combined,
                       m=M,
                       predictorMatrix=pred_matrix,
                       method=imp_method,
                       printFlag=F)
  
  ################ Step 3: Stack M imputed datasets
  stack = mice::complete(imputes, action="long", include = FALSE)
  stack$B2 = as.numeric(as.character(stack$B2))
  
  ################ Step 4: Calculate population-specific weights
  ##### Internal data only
  fit.gamma.I = glm(Y ~ X1 + X2 + X3 + X4 + B1 + B2, data = datan, family=binomial(link='logit'))
  gamma.I = coef(fit.gamma.I)
  
  ######## calculate the initial gamma for population S=1
  gamma.S1.origin= Initial.estimates(datan=datan, 
                                     gamma.I=gamma.I, 
                                     X=c('X1','X2','X3'), 
                                     B=c('X4','B1','B2'), 
                                     beta=betaHatExt_list[['Ext1']], 
                                     Btype=c('continuous','continuous','binary'))
  
  ###### ***Best fitted GLM for the external risk calculator 2:
  ###### ***When the external model is a non-parametric model with unknown form, this extra step is needed to calculate the self-generated beta estimates
  data.E2 = data.frame(matrix(ncol = 8, nrow = n*100))
  colnames(data.E2) = c('Y', 'X1', 'X2', 'X3', 'X4', 'B1', 'B2', 'S')
  
  data.E2$X1 = rep(datan$X1, 100)
  data.E2$X2 = rep(datan$X2, 100)
  data.E2$X3 = rep(datan$X3, 100)
  data.E2$X4 = rep(datan$X4, 100)
  data.E2$B1 = data.E2$B2 = NA
  data.E2$Y = rbinom(n*100, 1, predict(fit.rf, newdata=data.E2[,c('X1', 'X2', 'X3', 'X4')]))
  data.E2$Y[is.na(data.E2$Y)]=1
  
  fit.E2 = glm(Y ~ X1 + X2 + X3 + X4, data = data.E2, family=binomial(link='logit'))
  beta.E2 = coef(fit.E2)
  names(beta.E2) = c('int', 'X1', 'X2', 'X3', 'X4')
  betaHatExt_list = list(Ext1 = beta.E1, Ext2 = beta.E2)
  
  ###### calculate the initial gamma for population S=2
  gamma.S2.origin = Initial.estimates(datan=datan, 
                                      gamma.I=gamma.I, 
                                      beta=betaHatExt_list[['Ext2']], 
                                      X=c('X1','X2','X3','X4'), 
                                      B=c('B1','B2'), 
                                      Btype=c('continuous','binary'))
  
  ############ calculate weights for each observation by population
  stack$wt = 0
  stack[stack$S == 0, 'wt'] = dbinom(stack[stack$S == 0, 'Y'], 1, prob = expit(data.matrix(stack[stack$S==0, c('Int','X1','X2','X3','X4','B1','B2')])%*%matrix(gamma.I, q, 1)))
  stack[stack$S == 1, 'wt'] = dbinom(stack[stack$S == 1, 'Y'], 1, prob = expit(data.matrix(stack[stack$S==1, c('Int','X1','X2','X3','X4','B1','B2')])%*%matrix(gamma.S1.origin, q, 1)))
  stack[stack$S == 2, 'wt'] = dbinom(stack[stack$S == 2, 'Y'], 1, prob = expit(data.matrix(stack[stack$S==2, c('Int','X1','X2','X3','X4','B1','B2')])%*%matrix(gamma.S2.origin, q, 1)))
  stack = as.data.frame(stack %>% group_by(.id) %>% mutate(wt = wt / sum(wt)))   ## weights need to be re-scaled to 1 within individuals
  if(sum(is.na(stack$wt))>0){
    stack[is.na(stack$wt)==TRUE,]$wt = 0
  }
  
  ################ Step 5: Point Estimation
  fit = glm(Y~X1+X2+X3+X4+B1+B2+factor(S), data=stack, family=binomial(link='logit'), weights = stack$wt)
  
  return(coef(fit))
}

########################################################################################
########## Function 5: Resampling function to get bootstrap variance for continuous Y
#######################################################################################

#' Resample for bootstrap variance continuous Y
#' 
#' Resampling function to get bootstrap variance for continuous Y.
#' 
#' @param data synthetic data
#' @param indices indices ??
#' 
#' @return numeric vector of regression coefficients
#' 
#' @references Reference: Gu, T., Taylor, J.M.G. and Mukherjee, B. (2021) Regression 
#' inference for multiple populations by integrating summary-level data using stacked 
#' imputations <arXiv:http://arxiv.org/abs/2106.06835>. 
#' 
#' @export

Resample.gamma.continuousY <- function(data, indices){
  d = data[indices,]
  data.n = d
  
  ################ Step 1: convert the external model information into the synthetic data
  #### Function Create.Synthetic() can create synthetic data for both parametric model 1 and model 2
  data.combined = Create.Synthetic(datan=datan, 
                                   nrep=nrep, 
                                   Y='Y', 
                                   XB = c('X1','X2','B1','B2'), 
                                   Ytype='continuous',
                                   parametric = c('Yes','Yes'),
                                   betaHatExt_list=betaHatExt_list, 
                                   sigmaHatExt_list=sigmaHatExt_list)
  
  ################ Step 2: Multiple imputation (impute missingness ignoring the outcome Y)
  pred_matrix = mice::make.predictorMatrix(data.combined)
  pred_matrix[c('Int',"Y",'X1','S'),] = 0
  pred_matrix[,c('Int','Y','S')] = 0
  imp_method = mice::make.method(data.combined)
  imp_method[c('X2','B1','B2')] = c('norm','norm','logreg') #choose your desired imputation method for each missingness
  data.combined$B2 = factor(data.combined$B2)
  imputes = mice::mice(data.combined,
                       m=M,
                       predictorMatrix=pred_matrix,
                       method=imp_method,
                       printFlag=F)
  
  ################ Step 3: Stack M imputed datasets
  stack = mice::complete(imputes, action="long", include = FALSE)
  stack$B2 = as.numeric(as.character(stack$B2))
  
  ################ Step 4: Calculate population-specific weights
  ##### Internal data only
  fit.gamma.I = glm(Y ~ X1 + X2 + B1 + B2, data = datan, family=gaussian())
  gamma.I = coef(fit.gamma.I)
  
  ######## calculate the initial gamma for population S=1 and S=2
  gamma.S1.origin = Initial.estimates(datan=datan, 
                                      gamma.I=gamma.I, 
                                      beta=betaHatExt_list[['Ext1']], 
                                      X='X1', 
                                      B=c('X2','B1','B2'), 
                                      Btype=c('continuous','continuous','binary'))
  gamma.S2.origin = Initial.estimates(datan=datan, 
                                      gamma.I=gamma.I, 
                                      beta=betaHatExt_list[['Ext2']], 
                                      X=c('X1','X2'), 
                                      B=c('B1','B2'), 
                                      Btype=c('continuous','binary'))
  
  ############ calculate weights for each observation by population
  stack$wt = 0
  stack[stack$S == 0, 'wt'] = dnorm(stack[stack$S == 0, 'Y'], data.matrix(stack[stack$S==0, c('Int','X1','X2','B1','B2')])%*%matrix(gamma.I, q, 1))
  stack[stack$S == 1, 'wt'] = dnorm(stack[stack$S == 1, 'Y'], data.matrix(stack[stack$S==1, c('Int','X1','X2','B1','B2')])%*%matrix(gamma.S1.origin, q, 1))
  stack[stack$S == 2, 'wt'] = dnorm(stack[stack$S == 2, 'Y'], data.matrix(stack[stack$S==2, c('Int','X1','X2','B1','B2')])%*%matrix(gamma.S2.origin, q, 1))
  stack = as.data.frame(stack %>% group_by(.id) %>% mutate(wt = wt / sum(wt)))   ## weights need to be re-scaled to 1 within individuals
  if(sum(is.na(stack$wt))>0){
    stack[is.na(stack$wt)==TRUE,]$wt = 0
  }
  
  ################ Step 5: Point Estimation
  fit = glm(Y~X1+X2+B1+B2+factor(S), data=stack, family=gaussian(), weights = stack$wt)
  
  return(coef(fit))
}
