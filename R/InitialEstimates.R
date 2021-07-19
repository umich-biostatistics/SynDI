
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

#' Internal estimation
#' 
#' Calculate the initial estimates for external populations.
#' 
#' @param datan internal data only
#' @param gamma.I regression estimates using internal data only (datan)
#' @param X a vector of predictor that were used in the external study, e.g. X = 
#' c('X1','X2','X3')
#' @param B a vector of covariates that were not used in the external study, e.g. 
#' B=c('X4','B1','B2')
#' @param beta a vector of external model estimates, the vector order should be 
#' the same as listed in X, e.g. names(beta) = c("int", "X1", "X2", "X3")
#' @param Btype a vector of type of B, either continuous or binary. If "continuous", 
#' linear regression will be used; if "binary", logistic regression will be used. 
#' More types can be implemented manually.
#' 
#' @references 
#' Neuhaus, J. and Jewell, N. (1993). A geometric approach to assess bias due to 
#' omitted covariates in generalized linear models. Biometrika 80,807–815.
#' 
#' Gu, T., Taylor, J.M.G. and Mukherjee, B. (2021) Regression 
#' inference for multiple populations by integrating summary-level data using stacked 
#' imputations [arxiv.org/abs/2106.06835](http://arxiv.org/abs/2106.06835).
#' 
#' @examples 
#' \donttest{
#' library(mvtnorm)
#' library(mice)
#' library(arm)
#' library(dplyr)
#' library(StackImpute)
#' library(randomForest)
#' library(boot)
#' library(broom)
#' 
#' k = 2  #total number of external models
#' n = 1000  #internal study size
#' nrep = 3 #when generating the synthetic data, replicate the observed X 
#' #nrep times 
#' m1 = 10000  #size of external study 1
#' m2 = 10000  #size of external study 2
#' p1 = 3 #length of X for external calculator 1, excluding the intercept
#' p2 = 4 #length of X for external calculator 2, excluding the intercept
#' q = 7 #length of (X,B) in the full model, including the intercept
#' M = 100 #number of multiple imputation, denoted as M in the manuscript
#' 
#' gamma.S0.true = c(-1, rep(-1,6))   #true target model parameter for the internal 
#' #population
#' gamma.S1.true = c(2, rep(-1,6))    #true target model parameter for external 
#' #population 1
#' gamma.S2.true = c(3, rep(-1,6),-0.1,-0.1) ##true target model parameter for 
#' #external population 2, including non-linear terms -0.1*X1^2-0.1*X2*X3
#' set.seed(2333)
#' 
#' #########Population 1 has (Y,X1)
#' data.m1 = data.frame(matrix(ncol = 7, nrow = m1))
#' names(data.m1) = c('Y', 'X1', 'X2','X3', 'X4', 'B1','B2')
#' data.m1[,c('X1', 'X2', 'X3')] = MASS::mvrnorm(m1, rep(0,3), diag(0.7,3)+0.3)
#' data.m1$X4 = rnorm(m1, mean = 0.2*(data.m1$X1+data.m1$X2+data.m1$X3), sd=1)
#' data.m1$B1 = rnorm(m1, mean = 0.2*(data.m1$X1+data.m1$X2+data.m1$X3)+0.1*
#'                      data.m1$X4, sd=1)
#' data.m1[,c('X1', 'X2', 'X3', 'X4', 'B1')] = apply(data.m1[,c('X1', 'X2', 'X3', 
#'                                      'X4', 'B1')], 2, function(x) x-mean(x))
#' data.m1$B2 = rbinom(m1, 1, expit(0.2*(data.m1$X1+data.m1$X2+data.m1$X3)+0.1*
#'                                    (data.m1$X4+data.m1$B1)))
#' data.m1$Y = rbinom(m1, 1, prob = expit(data.matrix(cbind(1, 
#'                    data.m1[,c('X1', 'X2', 'X3', 'X4', 'B1','B2')])) %*% 
#'                                          matrix(gamma.S1.true, q, 1)))
#' #sum(data.m1$Y)/m1   ###prevalence=0.6454
#' 
#' fit.E1 = glm(Y~X1+X2+X3, data = data.m1, family = binomial(link='logit'))
#' 
#' ####### Beta estimates of external model 1
#' beta.E1 = coef(fit.E1)
#' names(beta.E1) = c('int', 'X1','X2','X3')
#' betaHatExt_list = list(Ext1 = beta.E1)
#' 
#' #######Population 2 has (Y,X1,X2,X3,X4,X5)
#' data.m2 = data.frame(matrix(ncol = 7, nrow = m2))
#' names(data.m2) = c('Y', 'X1', 'X2','X3', 'X4', 'B1','B2')
#' data.m2[,c('X1', 'X2', 'X3')] = MASS::mvrnorm(m2, rep(0,3), diag(0.7,3)+0.3)
#' data.m2$X4 = rnorm(m2, mean = 0.2*(data.m2$X1+data.m2$X2+data.m2$X3), sd=1)
#' data.m2$B1 = rnorm(m2, mean = 0.2*(data.m2$X1+data.m2$X2+data.m2$X3)+0.1*
#'                      data.m2$X4, sd=1)
#' data.m2[,c('X1', 'X2', 'X3', 'X4', 'B1')] = apply(data.m2[,c('X1', 'X2', 
#'                            'X3', 'X4', 'B1')], 2, function(x) x-mean(x))
#' data.m2$B2 = rbinom(m2, 1, expit(0.2*(data.m2$X1+data.m2$X2+data.m2$X3)+
#'                                    0.1*(data.m2$X4+data.m2$B1)))
#' data.m2$X1.sqr = data.m2$X1*data.m2$X1
#' data.m2$X2.X3 = data.m2$X2*data.m2$X3
#' data.m2$Y = rbinom(m2, 1, prob = expit(data.matrix(cbind(1, data.m2[,c('X1', 
#'                        'X2', 'X3', 'X4', 'B1','B2','X1.sqr','X2.X3')])) %*% 
#'                                          matrix(gamma.S2.true, q+2, 1)))
#' 
#' #######Risk calculator 2
#' fit.rf = randomForest(Y~X1+X2+X3+X4, data=data.m2)
#' 
#' datan = data.frame(matrix(ncol = 8, nrow = n))
#' names(datan) = c('Y', 'X1', 'X2', 'X3', 'X4', 'B1', 'B2','S')
#' datan[,c('X1', 'X2', 'X3')] = MASS::mvrnorm(n, rep(0,3), diag(0.7,3)+0.3)
#' datan$X4 = rnorm(n, mean = 0.2*(datan$X1+datan$X2+datan$X3), sd=1)
#' datan$B1 = rnorm(n, mean = 0.2*(datan$X1+datan$X2+datan$X3)+0.1*datan$X4, sd=1)
#' datan[,c('X1', 'X2', 'X3', 'X4', 'B1')] = apply(datan[,c('X1', 'X2', 'X3', 
#'                                    'X4', 'B1')], 2, function(x) x-mean(x))
#' datan$B2 = rbinom(n, 1, prob = expit(0.2*(datan$X1+datan$X2+datan$X3)+
#'                                        0.1*(datan$X4+datan$B1)))
#' datan$Y = rbinom(n, 1, prob = expit(data.matrix(cbind(1, datan[,c('X1', 'X2', 
#'                   'X3', 'X4', 'B1','B2')])) %*% matrix(gamma.S0.true, q, 1)))
#' 
#' #### Function Create.Synthetic() can only create synthetic data for parametric 
#' # model 1
#' data.combined = Create.Synthetic(nrep=nrep, 
#'                                  datan=datan, 
#'                                  Y='Y', 
#'                                  XB = c('X1','X2','X3','X4','B1','B2'), 
#'                                  Ytype='binary',
#'                                  parametric = c('Yes','No'),
#'                                  betaHatExt_list=betaHatExt_list, 
#'                                  sigmaHatExt_list=NULL)
#'                                  
#' #### For the external risk calculator 2 (k=2), we need to manually create the 
#' # synthetic data, and then row-bind with the dataset data.combined
#' X.2 = c('X1','X2','X3','X4') #predictors used in the external risk calculator 2
#' B.2 = c('B1','B2')  #unobserved variables in the external risk calculator 2
#' n = dim(datan)[1]   #total size of the internal data
#' m = n*nrep          #total number of synthetic data for S=2
#' data.2 = data.frame(matrix(ncol = 8, nrow = m))
#' colnames(data.2) = c('Y', 'X1', 'X2', 'X3', 'X4', 'B1', 'B2', 'S')
#' data.2[,B.2] = NA  #unobserved variables remain missing 
#' #replicate the observed X.2 nrep times
#' data.2[,X.2] = do.call("rbind", replicate(nrep, datan[,X.2], simplify = FALSE)) 
#' #generate synthetic Y values via the external risk calculator 2 
#' #(random forest model)
#' data.2$Y = rbinom(m, 1, predict(fit.rf, newdata=data.2[,X.2])) 
#' data.2$Y[is.na(data.2$Y)] = 1
#' data.2$S = 2
#' data.2 = cbind("Int" = 1, data.2)
#' 
#' ##### combine all synthetic data
#' data.combined = rbind(data.combined, data.2)
#' 
#' ### Impute missingness ignoring the outcome Y
#' pred_matrix = mice::make.predictorMatrix(data.combined)
#' pred_matrix[c('Int',"Y",'X1','X2','X3','S'),] = 0
#' pred_matrix[,c('Int','Y','S')] = 0
#' imp_method = mice::make.method(data.combined)
#' imp_method[c('X4','B1','B2')] = c('norm','norm','logreg') #choose your desired 
#' # imputation method for each missingness
#' data.combined$B2 = factor(data.combined$B2)
#' imputes = mice::mice(data.combined,
#'                      m=M,
#'                      predictorMatrix=pred_matrix,
#'                      method=imp_method,
#'                      printFlag=F)
#' 
#' stack = mice::complete(imputes, action="long", include = FALSE)
#' stack$B2 = as.numeric(as.character(stack$B2))
#' 
#' ##### Internal data only
#' fit.gamma.I = glm(Y ~ X1 + X2 + X3 + X4 + B1 + B2, data = datan, 
#'                   family=binomial(link='logit'))
#' gamma.I = coef(fit.gamma.I)
#' 
#' ######## calculate the initial gamma for population S=1
#' gamma.S1.origin= Initial.estimates(datan=datan, 
#'                                    gamma.I=gamma.I, 
#'                                    X=c('X1','X2','X3'), 
#'                                    B=c('X4','B1','B2'), 
#'                                    beta=betaHatExt_list[['Ext1']], 
#'                                    Btype=c('continuous','continuous','binary'))
#' 
#' }
#' 
#' @export
Initial.estimates <- function(datan, gamma.I, X, B, beta, Btype){
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
    tmp1 = gamma.int + sum(gamma.I[B] * E.B)
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
