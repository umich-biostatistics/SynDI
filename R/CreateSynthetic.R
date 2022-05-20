
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
#' imputations [https://arxiv.org/abs/2106.06835](https://arxiv.org/abs/2106.06835). 
#' 
#' @examples 
#' data(create_synthetic_example)
#' 
#' nrep = create_synthetic_example$nrep
#' datan = create_synthetic_example$datan
#' betaHatExt_list = create_synthetic_example$betaHatExt_list
#' 
#' data.combined = Create.Synthetic(nrep = nrep, datan = datan, Y = 'Y', 
#'     XB = c('X1', 'X2', 'X3', 'X4', 'B1', 'B2'), Ytype = 'binary', 
#'     parametric = c('Yes', 'No'), betaHatExt_list = betaHatExt_list, 
#'     sigmaHatExt_list = NULL)
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
