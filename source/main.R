#-------------------------------------
# simulation for regression
# source function 
source("~/Desktop/Lei/AIMER2/github/sim.R")
source("~/Desktop/Lei/AIMER2/github/find_t.R")
source("~/Desktop/Lei/AIMER2/github/sparsePC.R")
source("~/Desktop/Lei/AIMER2/github/sparsePC_reg.R")
source("~/Desktop/Lei/AIMER2/github/compare_reg_sim.R")

# set parameter
n = 100 #
p = 1000 #
s = 10 #
d = 5 #
SNRx = 10 #
SNRy = 10 #
d.est = 5 #
nsol = 10 #
niter = 50

# output holder
SEC = matrix(NA, nrow = niter, ncol = 7)
MSE = matrix(NA, nrow = niter, ncol = 7)
NCOV = matrix(NA, nrow = niter, ncol = 7)
COV.PRECISION = matrix(NA, nrow = niter, ncol = 7)
COV.RECALL = matrix(NA, nrow = niter, ncol = 7)

for (i in 1:niter) {
  print(i)
  # generate simulation data
  set.seed(492 + i)
  data = sim_reg(n = n, p = p, s = s, d = d, SNRx = SNRx, SNRy = SNRy)
  
  # compare results
  fit = compare_reg_sim(Xtrain = data$Xtrain, Ytrain = data$Ytrain, 
                        Xvalidate = data$Xvalidate, Yvalidate = data$Yvalidate,
                        Xtest = data$Xtest, Ytest = data$Ytest,
                        d = d.est, nsol = nsol, s = s)
  
  # store output
  SEC[i, ] = fit$sec
  MSE[i, ] = fit$mse
  NCOV[i, ] = fit$ncov
  COV.PRECISION[i, ] = fit$cov.precision
  COV.RECALL[i, ] = fit$cov.recall
}

# plot
par(mar=c(3,4,3,1), oma=c(0,0.5,0,0))
par(mfrow=c(1,4))
boxplot(log(MSE[,c(3, 4, 1, 2)]), main = "MSE", las = 1, horizontal = T, 
        cex.axis = 0.8, cex.main=0.9
        , names = c("Ridge", "Lasso", "sparsePC", "True"))
boxplot(log(NCOV[,c(3, 4, 1, 2)]), main = "# of Selected Features", las = 1, horizontal = T, 
        cex.axis = 0.8, cex.main=0.9
        , names = c("Ridge", "Lasso", "sparsePC", "True"))
boxplot(COV.PRECISION[,c(3, 4, 1, 2)], main = "Precision", las = 1, horizontal = T, 
        cex.axis = 0.8, cex.main=0.9
        , names = c("Ridge", "Lasso", "sparsePC", "True"))
boxplot(COV.RECALL[,c(3, 4, 1, 2)], main = "Recall", las = 1, horizontal = T, 
        cex.axis = 0.8, cex.main=0.9
        , names = c("Ridge", "Lasso", "sparsePC", "True"))
par(mfrow=c(1,1))




#-------------------------------------
# simulation for classification
# source comparison function 
source("~/Desktop/Lei/AIMER2/github/sim.R")
source("~/Desktop/Lei/AIMER2/github/find_t.R")
source("~/Desktop/Lei/AIMER2/github/sparsePC.R")
source("~/Desktop/Lei/AIMER2/github/sparsePC_class.R")
source("~/Desktop/Lei/AIMER2/github/compare_class_sim.R")

# set parameter
n = 100 #
p = 1000 #
s = 10 #
d = 5 #
SNRx = 10 #
SNRy = 10 #
d.est = 5 #
nsol = 10 #
niter = 50

# output holder
SEC = matrix(NA, nrow = niter, ncol = 5)
ACCURACY = matrix(NA, nrow = niter, ncol = 5)
NCOV = matrix(NA, nrow = niter, ncol = 5)
COV.PRECISION = matrix(NA, nrow = niter, ncol = 5)
COV.RECALL = matrix(NA, nrow = niter, ncol = 5)

for (i in 1:niter) {
  print(i)
  # generate simulation data
  set.seed(492 + i)
  data = sim_class(n = n, p = p, s = s, d = d, SNRx = SNRx, SNRy = SNRy)
  
  # compare results
  fit = compare_class_sim(Xtrain=data$Xtrain, Ytrain=data$Ytrain,
                          Xvalidate=data$Xvalidate, Yvalidate=data$Yvalidate,
                          Xtest=data$Xtest, Ytest=data$Ytest, 
                          d = d.est, nsol = nsol, s = s)
  
  # store output
  SEC[i, ] = fit$sec
  ACCURACY[i, ] = fit$accuracy
  NCOV[i, ] = fit$ncov
  COV.PRECISION[i, ] = fit$cov.precision
  COV.RECALL[i, ] = fit$cov.recall
}

par(mar=c(3,4,3,1), oma=c(0,0.5,0,0))
par(mfrow=c(1,4))
boxplot(log(MSE[,c(3, 4, 1, 2)]), main = "MSE", las = 1, horizontal = T, 
        cex.axis = 0.8, cex.main=0.9
        , names = c("Ridge", "Lasso", "sparsePC", "True"))
boxplot(log(ACCURACY[,c(3, 4, 1, 2)]), main = "# of Selected Features", las = 1, horizontal = T, 
        cex.axis = 0.8, cex.main=0.9
        , names = c("Ridge", "Lasso", "sparsePC", "True"))
boxplot(COV.PRECISION[,c(3, 4, 1, 2)], main = "Precision", las = 1, horizontal = T, 
        cex.axis = 0.8, cex.main=0.9
        , names = c("Ridge", "Lasso", "sparsePC", "True"))
boxplot(COV.RECALL[,c(3, 4, 1, 2)], main = "Recall", las = 1, horizontal = T, 
        cex.axis = 0.8, cex.main=0.9
        , names = c("Ridge", "Lasso", "sparsePC", "True"))
par(mfrow=c(1,1))