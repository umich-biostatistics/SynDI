
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

#' Create the synthetic data
#'
#' Creates a synthetic data set from internal data and external models.
#' 
#' @param datan internal data only
#' @param nrep number of replication when creating the synthetic data
#' @param Y outcome name, e.g. Y='Y'
#' @param XB all covariate names for both X and B in the target model, e.g. XB=c('X1','X2','X3','X4','B1','B2')
#' @param Ytype the type of outcome Y, either 'binary' or 'continuous'.
#' @param parametric choice of "Yes" or "No" for each external model. Specify whether the external model is paramtric or not, e.g. parametric=c('Yes','No')
#' @param betaHatExt_list a list of parameter estimates of the external models. The order needs to be the same as listed in XB, and variable name is required. See example for details.
#' @param sigmaHatExt_list a list of sigma^2 for continuous outcome fitted from linear regression. If not available or the outcome type is binary, set sigmaHatExt_list=NULL
#'
#' @references Reference: Gu, T., Taylor, J.M.G. and Mukherjee, B. (2021) Regression 
#' inference for multiple populations by integrating summary-level data using stacked 
#' imputations [arxiv.org/abs/2106.06835](http://arxiv.org/abs/2106.06835). 
#' 
#' @examples 
#' \donttest{
#' library(MASS)
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
#'                        'X3', 'X4', 'B1')], 2, function(x) x-mean(x))
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
#'                                                          'X4', 'B1')], 2, function(x) x-mean(x))
#' datan$B2 = rbinom(n, 1, prob = expit(0.2*(datan$X1+datan$X2+datan$X3)+
#'                                        0.1*(datan$X4+datan$B1)))
#' datan$Y = rbinom(n, 1, prob = expit(data.matrix(cbind(1, datan[,c('X1', 'X2', 
#'                                                                   'X3', 'X4', 'B1','B2')])) %*% matrix(gamma.S0.true, q, 1)))
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
#' }
#'
#' @export
Create.Synthetic <- function(datan, nrep, Y, XB, Ytype = 'binary', parametric, betaHatExt_list, sigmaHatExt_list = NULL){
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
