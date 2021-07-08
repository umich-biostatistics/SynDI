##### Table of Contents:
# Function 1: Expit function
# Function 2: Create the synthetic data
# Function 3: Calculate the initial estimates for external populations
# Function 4: Resampling function to get bootstrap variance for binary Y in Example 1
# Function 5: Resampling function to get bootstrap variance for continuous Y in Example 2

# Reference: Gu, T., Taylor, J.M.G. and Mukherjee, B. (2021) Regression inference for multiple populations by integrating summary-level data using stacked imputations
#            <arXiv:http://arxiv.org/abs/2106.06835>. 

####################################
########## Function 1: Expit function
###################################
expit <- function(x){
  y = exp(x)/(1+exp(x))
  return(y)
}

##################################################
########## Function 2: Create the synthetic data
#################################################
# datan: internal data only
# nrep: number of replication when creating the synthetic data
# Y: outcome name, e.g. Y='Y'
# XB: all covariate names for both X and B in the target model, e.g. XB=c('X1','X2','X3','X4','B1','B2')
# Ytype: the type of outcome Y, either 'binary' or 'continuous'.
# parametric: choice of "Yes" or "No" for each external model. Specify whether the external model is paramtric or not, e.g. parametric=c('Yes','No')
# betaHatExt_list: a list of parameter estimates of the external models. The order needs to be the same as listed in XB, and variable name is required. See example for details.
# sigmaHatExt_list: a list of sigma^2 for continuous outcome fitted from linear regression. If not available or the outcome type is binary, set sigmaHatExt_list=NULL
# reference: Gu, T., Taylor, J.M.G., and Mukherjee, B. (2019). Synthetic data method to incorporate external information into a current study. Canadian Journal of Statistics47, 580–603.

Create.Synthetic <- function(datan, nrep, Y='Y', XB=c('X1','X2','X3','X4','B1','B2'), Ytype='binary', parametric=c('Yes','No'), betaHatExt_list, sigmaHatExt_list=NULL){
  n = dim(datan)[1]
  m = n*nrep
  data = data.frame(matrix(ncol = length(XB)+2, nrow = 1))
  colnames(data) = c(Y, XB, 'S')
  for(k in 1:length(betaHatExt_list)){
    if(parametric[k] == 'Yes'){
      X.k = names(betaHatExt_list[[paste0('Ext',k)]])[-1]
      B.k = base::setdiff(XB, X.k)
      
      data.k = data.frame(matrix(ncol = length(XB)+1, nrow = m))
      colnames(data.k) = c('Y', XB)
      data.k[,B.k] = NA
      data.k[,X.k] = do.call("rbind", replicate(nrep, data.frame(datan[,X.k]), simplify = FALSE))
      
      #### If binary, generate Y~Bernoulli(p); if continuous, generate Y~Normal
      if(Ytype=='binary'){
        
        data.k[,Y] = rbinom(m, 1, prob = expit(data.matrix(cbind(1, data.k[,X.k])) %*% as.matrix(betaHatExt_list[[k]], length(betaHatExt_list[[k]]), 1)))
        
      }else if(Ytype=='continuous'){
        
        #### If sigma^2 is not available from the external model, estimate using the internal data
        if(is.null(sigmaHatExt_list)){
          
          sigmaHatExt = 1/n*sum((datan[,Y] - datan[,X.k]%*%betaHatExt_list[[k]])^2)
          
        }else{
          
          sigmaHatExt = sigmaHatExt_list[[k]]
          
        }
        data.k[,Y] = rnorm(m, data.matrix(cbind(1, data.k[,X.k])) %*% as.matrix(betaHatExt_list[[k]], length(betaHatExt_list[[k]]), 1), sigmaHatExt)
      }
      data.k$S = k
      data = rbind(data, data.k) 
    }else if(parametric[k] == 'No'){
      next
    }
  }
  data = data[-1,]
  
  datan$S = 0
  combined = rbind(datan, data)
  combined = cbind("Int" = 1, combined)
  
  return(combined)
}


###################################################################################
########## Function 3: Calculate the initial estimates for external populations
##################################################################################
# datan: internal data only
# gamma.I: regression estimates using internal data only (datan)
# X: a vector of predictor that were used in the external study, e.g. X = c('X1','X2','X3')
# B: a vector of covariates that were not used in the external study, e.g. B=c('X4','B1','B2')
# beta: a vector of external model estimates, the vector order should be the same as listed in X, e.g. names(beta) = c("int", "X1", "X2", "X3")
# Btype: a vector of type of B, either continuous or binary. If "continuous", linear regression will be used; if "binary", logistic regression will be used. More types can be implemented manually.
# reference: Neuhaus, J. and Jewell, N. (1993). A geometric approach to assess bias due to omitted covariates in generalized linear models. Biometrika 80,807–815.

Initial.estimates <- function(datan, gamma.I, X=c('X1','X2','X3'), B=c('X4','B1','B2'), beta=betaHatExt_list[['Ext1']], Btype=c('continuous','continuous','binary')){
  E.B = Var.B = NULL
  E.B.X = delta.B.binary = list()
  for(i in 1:length(B)){
    XB = c(X,B)
    RHS = XB[1:(length(X)+i-1)]
    RHS.eq = X[1]
    if(length(X)>1 | i>1){
      for(p in 2:(length(X)+i-1)){
        RHS.eq = paste0(RHS.eq, '+', XB[p])
      }
    }
    equation = as.formula(paste0(B[i],'~',RHS.eq))

    if(Btype[i]=='continuous'){
      fit = lm(equation, data=datan)
      E.Bi.X = coef(fit)
      Var.Bi = sigma(fit)
      E.Bi = sum(c(1, colMeans(data.frame(datan[,RHS])))*E.Bi.X)
    }else if(Btype[i]=='binary'){
      fit = glm(equation, data=datan, family = binomial(link='logit'))
      E.Bi.X = coef(fit)
      E.Bi = expit(sum(c(1,colMeans(datan[,RHS]))*E.Bi.X))
      Var.Bi = E.Bi*(1-E.Bi)
      
      #####calculate delta for the binary variable
      delta.B = NULL
      for(p.X in 1:length(X)){
        add = rep(0, length(RHS))
        add[p.X] = 1
        delta.Bi = expit(sum(c(1, colMeans(datan[RHS])+add)*E.Bi.X))-expit(sum(c(1,colMeans(datan[RHS]))*E.Bi.X))
        delta.B = c(delta.B, delta.Bi)
      }
      delta.B.binary[[i]] = delta.B
    }
    E.B = c(E.B, E.Bi)
    Var.B = c(Var.B, Var.Bi)
    E.B.X[[i]] = E.Bi.X
  }
  
  ## covariance matrix
  if(length(B) > 1){
    tmp = diag(gamma.I[B]) %*% cov(datan[,B]) %*% diag(gamma.I[B])
    cov.Bi.Bj = tmp[upper.tri(tmp)] #gamma.B1*gamma.B2*cov.B1.B2
  }else if(length(B) == 1){
    cov.Bi.Bj = 0
  }
  SolveEquation = function(gamma.int){
    tmp1 = gamma.int+ sum(gamma.I[B] * E.B)
    RHS = expit(tmp1)*(1 + 1/2*(1-exp(tmp1))/(1+exp(tmp1))^2*(sum(gamma.I[B]^2*Var.B)+2*sum(cov.Bi.Bj)))
    LHS = expit(beta[1])
    return(RHS-LHS)
  }
  
  ###### intercept
  gamma.int = stats::uniroot(SolveEquation, interval = c(-10,10))$root

  ###### slopes
  tmp1 = gamma.int + sum(gamma.I[B] * E.B)
  E.mu0     = expit(tmp1)*(1 + 1/2*(1-exp(tmp1))/(1+exp(tmp1))^2*(sum(gamma.I[B]^2*Var.B)+2*sum(cov.Bi.Bj)))
  E.mu0.sqr = (expit(tmp1))^2*(1 + (2-exp(tmp1))/(1+exp(tmp1))^2*(sum(gamma.I[B]^2*Var.B)+2*sum(cov.Bi.Bj)))
  V.mu0 = E.mu0.sqr - E.mu0^2
  
  gamma.X = NULL
  for(p in 1:length(X)){
    tmp = NULL
    for(q in 1:length(B)){
      if(Btype[q]=='binary'){
        tmp = c(tmp, delta.B.binary[[q]][p])
      }else{
        count = q
        sum = E.B.X[[q]][p+1]
        while(count > 1){
          sum = sum + E.B.X[[q]][1+length(X)+(count-1)]*E.B.X[[count-1]][p+1] 
          count = count - 1
          }
        tmp = c(tmp, sum) #e.g. delta.B1.X1+delta.B1.X4*delta.X4.X1
      }
    }
    gamma.p = beta[p+1]/(1-V.mu0/(E.mu0*(1-E.mu0)))-sum(gamma.I[B]*tmp)
    gamma.X = c(gamma.X, gamma.p)
  }
  gamma.origin = c(gamma.int, gamma.X, gamma.I[B])
  return(gamma.origin)
}



###################################################################################
########## Function 4: Resampling function to get bootstrap variance for binary Y
##################################################################################
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



