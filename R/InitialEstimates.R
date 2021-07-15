
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
