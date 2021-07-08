rm(list = ls())
library(mvtnorm)
library(mice)
library(arm)
library(dplyr)
library(StackImpute)
library(randomForest)
library(boot)
library(broom)

source("SynDI.R")

##########################################################
##########################################################
# Example 1: Binary Y, simulation II in the manuscript 
##########################################################
##########################################################

################# Settings ###########################
k = 2  #total number of external models
n = 1000  #internal study size
nrep = 3 #when generating the synthetic data, replicate the observed X nrep times 
m1 = 10000  #size of external study 1
m2 = 10000  #size of external study 2
p1 = 3 #length of X for external calculator 1, excluding the intercept
p2 = 4 #length of X for external calculator 2, excluding the intercept
q = 7 #length of (X,B) in the full model, including the intercept
M = 100 #number of multiple imputation, denoted as M in the manuscript

gamma.S0.true = c(-1, rep(-1,6))   #true target model parameter for the internal population
gamma.S1.true = c(2, rep(-1,6))    #true target model parameter for external population 1
gamma.S2.true = c(3, rep(-1,6),-0.1,-0.1) ##true target model parameter for external population 2, including non-linear terms -0.1*X1^2-0.1*X2*X3

###############################################################################
###Create external models:
###obtain beta estimates from external model 1: Y~X1+X2+X3
###fit random forest model for external model 2: Y~X1+X2+X3+X4
###############################################
set.seed(2333)

#########Population 1 has (Y,X1)
data.m1 = data.frame(matrix(ncol = 7, nrow = m1))
names(data.m1) = c('Y', 'X1', 'X2','X3', 'X4', 'B1','B2')
data.m1[,c('X1', 'X2', 'X3')] = MASS::mvrnorm(m1, rep(0,3), diag(0.7,3)+0.3)
data.m1$X4 = rnorm(m1, mean = 0.2*(data.m1$X1+data.m1$X2+data.m1$X3), sd=1)
data.m1$B1 = rnorm(m1, mean = 0.2*(data.m1$X1+data.m1$X2+data.m1$X3)+0.1*data.m1$X4, sd=1)
data.m1[,c('X1', 'X2', 'X3', 'X4', 'B1')] = apply(data.m1[,c('X1', 'X2', 'X3', 'X4', 'B1')], 2, function(x) x-mean(x))
data.m1$B2 = rbinom(m1, 1, expit(0.2*(data.m1$X1+data.m1$X2+data.m1$X3)+0.1*(data.m1$X4+data.m1$B1)))
data.m1$Y = rbinom(m1, 1, prob = expit(data.matrix(cbind(1, data.m1[,c('X1', 'X2', 'X3', 'X4', 'B1','B2')])) %*% matrix(gamma.S1.true, q, 1)))
#sum(data.m1$Y)/m1   ###prevalence=0.6454

fit.E1 = glm(Y~X1+X2+X3, data = data.m1, family = binomial(link='logit'))

####### Beta estimates of external model 1
beta.E1 = coef(fit.E1)
names(beta.E1) = c('int', 'X1','X2','X3')
betaHatExt_list = list(Ext1 = beta.E1)

#######Population 2 has (Y,X1,X2,X3,X4,X5)
data.m2 = data.frame(matrix(ncol = 7, nrow = m2))
names(data.m2) = c('Y', 'X1', 'X2','X3', 'X4', 'B1','B2')
data.m2[,c('X1', 'X2', 'X3')] = MASS::mvrnorm(m2, rep(0,3), diag(0.7,3)+0.3)
data.m2$X4 = rnorm(m2, mean = 0.2*(data.m2$X1+data.m2$X2+data.m2$X3), sd=1)
data.m2$B1 = rnorm(m2, mean = 0.2*(data.m2$X1+data.m2$X2+data.m2$X3)+0.1*data.m2$X4, sd=1)
data.m2[,c('X1', 'X2', 'X3', 'X4', 'B1')] = apply(data.m2[,c('X1', 'X2', 'X3', 'X4', 'B1')], 2, function(x) x-mean(x))
data.m2$B2 = rbinom(m2, 1, expit(0.2*(data.m2$X1+data.m2$X2+data.m2$X3)+0.1*(data.m2$X4+data.m2$B1)))
data.m2$X1.sqr = data.m2$X1*data.m2$X1
data.m2$X2.X3 = data.m2$X2*data.m2$X3
data.m2$Y = rbinom(m2, 1, prob = expit(data.matrix(cbind(1, data.m2[,c('X1', 'X2', 'X3', 'X4', 'B1','B2','X1.sqr','X2.X3')])) %*% matrix(gamma.S2.true, q+2, 1)))

#######Risk calculator 2
fit.rf = randomForest(Y~X1+X2+X3+X4, data=data.m2)

#########################################
#Create an internal dataset of size n
#########################################
datan = data.frame(matrix(ncol = 8, nrow = n))
names(datan) = c('Y', 'X1', 'X2', 'X3', 'X4', 'B1', 'B2','S')
datan[,c('X1', 'X2', 'X3')] = MASS::mvrnorm(n, rep(0,3), diag(0.7,3)+0.3)
datan$X4 = rnorm(n, mean = 0.2*(datan$X1+datan$X2+datan$X3), sd=1)
datan$B1 = rnorm(n, mean = 0.2*(datan$X1+datan$X2+datan$X3)+0.1*datan$X4, sd=1)
datan[,c('X1', 'X2', 'X3', 'X4', 'B1')] = apply(datan[,c('X1', 'X2', 'X3', 'X4', 'B1')], 2, function(x) x-mean(x))
datan$B2 = rbinom(n, 1, prob = expit(0.2*(datan$X1+datan$X2+datan$X3)+0.1*(datan$X4+datan$B1)))
datan$Y = rbinom(n, 1, prob = expit(data.matrix(cbind(1, datan[,c('X1', 'X2', 'X3', 'X4', 'B1','B2')])) %*% matrix(gamma.S0.true, q, 1)))

##################################################################################
### Step 1: convert the external model information into the synthetic data
##################################################################################
#### Function Create.Synthetic() can only create synthetic data for parametric model 1
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

##################################################################################
### Step 2: Multiple imputation
##################################################################################
### Impute missingness ignoring the outcome Y
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

##################################################################################
### Step 3: Stack M imputed datasets
##################################################################################
stack = mice::complete(imputes, action="long", include = FALSE)
stack$B2 = as.numeric(as.character(stack$B2))

##################################################################################
### Step 4: Calculate population-specific weights
##################################################################################
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

##################################################################################
### Step 5: Estimation
##################################################################################

#### point estimation
fit = glm(Y~X1+X2+X3+X4+B1+B2+factor(S), data=stack, family=binomial(link='logit'), weights = stack$wt)
gamma.proposed = coef(fit)

##### Lauren's variance estimation
Louiscovar = StackImpute::Louis_Information(fit, stack, M = M)
var.Louis = diag(solve(Louiscovar))

##### Bootstrap variance 
#***readers need to modify the existing function Resample.gamma.binaryY() to match their own Steps 1-5
results.boot = boot(data=datan, 
                    statistic=Resample.gamma.binaryY,
                    R=2) ##R=2 is for illustration purpose to save running time. Typically a larger R number, e.g.R=500 is used for reliable estimates.***it may take hours to finish running 
broom::tidy(results.boot)$std.error^2

