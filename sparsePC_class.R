# install packages from github for first time user
# library(devtools)
# install_github("Lei-D/irlba", ref="lei")
# install_github("Lei-D/fps", ref="fps_irlba")
library(irlba)
library(fps)
library(psych) # for logistic function

# this function compute sparsePC for classification
sparsePC_class = function(Xtrain, Ytrain, Xvalidate, Yvalidate, Xtest, Ytest,
                          d, nsol, screening = TRUE){
  n = nrow(Xtrain)
  p = ncol(Xtrain)
  
  # standardize data matrix
  xmeans = colMeans(Xtrain)
  xsd = apply(Xtrain, 2, sd)
  xsd = ifelse(xsd == 0, 1, xsd)
  Xtrain = scale(Xtrain, center = xmeans, scale = xsd)
  Xvalidate = scale(Xvalidate, center = xmeans, scale = xsd)
  Xtest = scale(Xtest, center = xmeans, scale = xsd)
  
  # generate sample covariance matrix
  S = crossprod(Xtrain) / n
  
  # compute projection matrix
  approximate = fps(S, ndim = d, ncomp = d, 
                    lambda = seq(from = 1, to = 0, length.out = nsol),
                    maxnvar = p, maxiter = 100)
  
  # output holder
  mse.validate = rep(NA, nsol)
  accuracy.validate = rep(NA, nsol)
  mse.test = rep(NA, nsol)
  accuracy.test = rep(NA, nsol)
  Betahat = matrix(NA, nrow = nsol, ncol = p)
  
  for (i in 1:nsol) {
    if (screening == TRUE) {
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
    }
    
    if (all(approximate$projection[[i]] == 0) == F) {
      # decompose projection matrix
      decomp = svd(approximate$projection[[i]], nu = d, nv = 0)
      vhat = decomp$u
      vhat[which(abs(vhat) < 1e-10)] = 0
      pchat = Xtrain %*% vhat
      
      # fit logistic
      model = glm(Ytrain~pchat, family = binomial(link = 'logit'))
      gamhat = as.vector(model$coefficients[-1])
      intercept = model$coefficients[1]
      betahat = vhat %*% t(t(gamhat))
      Betahat[i, ] = betahat
      
      # find the best cutoff
      Ytrain.prob = logistic(intercept + Xtrain %*% betahat)
      cutoffs = sort(Ytrain.prob)
      accuracy = rep(NA, n)
      for (j in 1:n){
        prediction = ifelse(Ytrain.prob >= cutoffs[j], 1, 0)
        accuracy[j] = sum(Ytrain == prediction) / n * 100
      }
      cutoffs = cutoffs[which(accuracy == max(accuracy))]
      cutoff.best = cutoffs[which.min(abs(cutoffs - 0.5))]
      # cutoff.best = 0.5
      
      # predict on validation set
      Yvalidate.prob = logistic(intercept + Xvalidate %*% betahat)
      mse.validate[i] = mean((Yvalidate.prob - Yvalidate)^2)
      Yvalidate.hat = ifelse(Yvalidate.prob >= cutoff.best, 1, 0)
      accuracy.validate[i] = sum(Yvalidate == Yvalidate.hat) / length(Yvalidate) * 100
      
      # predict on test set
      Ytest.prob = logistic(intercept + Xtest %*% betahat)
      mse.test[i] = mean((Ytest.prob - Ytest)^2)
      Ytest.hat = ifelse(Ytest.prob >= cutoff.best, 1, 0)
      accuracy.test[i] = sum(Ytest == Ytest.hat) / length(Ytest) * 100
    } else {
      # predict by majority
      betahat = rep(0, p)
      Betahat[i, ] = betahat
      
      # predict on validation set
      Yvalidate.hat = rep(NA, length(Yvalidate))
      Yvalidate.hat[1:length(Yvalidate)] = ifelse(sum(Ytrain == 1) >= sum(Ytrain == 0), 
                                                  rep(1, length(Yvalidate)), 
                                                  rep(0, length(Yvalidate)))
      mse.validate[i] = mean((Yvalidate.hat - Yvalidate)^2)
      accuracy.validate[i] = sum(Yvalidate == Yvalidate.hat) / length(Yvalidate) * 100
      
      # predict on test set
      Ytest.hat = rep(NA, length(Ytest))
      Ytest.hat[1:length(Ytest)] = ifelse(sum(Ytrain==1)>=sum(Ytrain==0), 
                                          rep(1, length(Ytest)), 
                                          rep(0, length(Ytest)))
      mse.test[i] = mean((Ytest.hat - Ytest)^2)
      accuracy.test[i] = sum(Ytest == Ytest.hat) / length(Ytest) * 100
    }
  }
  
  # return output
  out = list(betahat = Betahat,
             mse.validate = mse.validate,
             accuracy.validate = accuracy.validate,
             mse.test = mse.test,
             accuracy.test = accuracy.test)
  return(out)
}