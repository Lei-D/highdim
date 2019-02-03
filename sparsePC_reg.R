# install packages from github for first time user
# library(devtools)
# install_github("Lei-D/irlba", ref="lei")
# install_github("Lei-D/fps", ref="fps_irlba")
library(irlba)
library(fps)

# this function train the model on training data
# then predict on validation and test data
sparsePC_reg = function(Xtrain, Ytrain, Xvalidate, Yvalidate, Xtest, Ytest,
                        d, nsol){
  n = nrow(Xtrain)
  p = ncol(Xtrain)
  
  # standardize data matrix
  xmeans = colMeans(Xtrain)
  xsd = apply(Xtrain, 2, sd)
  Xtrain = scale(Xtrain, center = xmeans, scale = xsd)
  Xvalidate = scale(Xvalidate, center = xmeans, scale = xsd)
  Xtest = scale(Xtest, center = xmeans, scale = xsd)
  
  # standardize response
  ymean = mean(Ytrain)
  ysd = sd(Ytrain)
  Ytrain = scale(Ytrain, center = ymean, scale = ysd)
  
  # generate sample covariance matrix
  S = crossprod(Xtrain) / n
  
  # compute sparse projection matrix
  approximate = fps(S, ndim = d, ncomp = d, 
                    lambda = seq(from = 1, to = 0, length.out = nsol),
                    maxnvar = p, maxiter = 100)
  
  # output holder
  mse.validate = rep(NA, nsol)
  mse.test = rep(NA, nsol)
  Betahat = matrix(NA, nrow = nsol, ncol = p)
  
  for (i in 1:nsol) {
    l = sort(diag(approximate$projection[[i]]))
    if (i == 1) {
      t = max(l)
    } else {
      t = find_t(l)
    }
    
    # set rows and columns that has small diagonal values to be 0
    small = which(diag(approximate$projection[[i]]) <= t)
    approximate$projection[[i]][small, ] = 0
    approximate$projection[[i]][ ,small] = 0
    
    if (all(approximate$projection[[i]] == 0) == F) {
      # decompose projection matrix
      decomp = svd(approximate$projection[[i]], nu = d, nv = 0)
      vhat = decomp$u
      vhat[which(abs(vhat) < 1e-10)] = 0
      pchat = Xtrain %*% vhat
      
      # fit linear model
      model = lm(Ytrain~pchat)
      gamhat = as.matrix(model$coefficients[-1], ncol = 1)
      betahat = vhat %*% gamhat
      Betahat[i, ] = betahat
      
    } else {
      # use ymean to predict
      betahat = rep(0, p)
      Betahat[i, ] = betahat
    }
    
    # predict on validation set
    Yvalidate.hat = ymean + Xvalidate %*% betahat * ysd
    mse.validate[i] = mean((Yvalidate.hat - Yvalidate)^2)
    
    # predict on test set
    Ytest.hat = ymean + Xtest %*% betahat * ysd
    mse.test[i] = mean((Ytest.hat - Ytest)^2)
  }
  
  # return output
  out = list(betahat = Betahat,
             mse.validate= mse.validate,
             mse.test = mse.test)
  return(out)
}