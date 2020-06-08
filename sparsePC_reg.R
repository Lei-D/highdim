# install packages from github for first time user
# library(devtools)
# install_github("Lei-D/irlba", ref="lei")
# install_github("Lei-D/fps", ref="fps_irlba")
library(irlba)
library(fps)


sparsePC_reg = function(Xtrain, Ytrain, Xvalidate, Yvalidate, Xtest, Ytest,
                       d, nsol, maxnvar, screening = TRUE){
  n = nrow(Xtrain)
  p = ncol(Xtrain)
  
  # center input data X
  xmeans = colMeans(Xtrain)
  xsd = apply(Xtrain, 2, sd)
  Xtrain = scale(Xtrain, center = xmeans, scale = xsd)
  Xvalidate = scale(Xvalidate, center = xmeans, scale = xsd)
  Xtest = scale(Xtest, center = xmeans, scale = xsd)
  
  # center input data Y
  ymean = mean(Ytrain)
  ysd = sd(Ytrain)
  Ytrain = scale(Ytrain, center = ymean, scale = ysd)
  
  # output holder
  mse.train = matrix(NA, nrow=length(d), ncol = nsol)
  mse.validate = matrix(NA, nrow=length(d), ncol = nsol)
  mse.test = matrix(NA, nrow=length(d), ncol = nsol)
  Betahat = matrix(NA, nrow = nsol*length(d), ncol = p+2)
  Vhat.norm = matrix(NA, nrow = nsol*length(d), ncol = p+2)
  
  # generate sample covariance matrix
  S = crossprod(Xtrain) / n
  
  for(k in 1:length(d)){
    # run fps
    approximate = fps(S, ndim = d[k], ncomp = d[k], 
                      lambda = seq(from = 1, to = 0, length.out = nsol),
                      maxnvar = maxnvar, maxiter = 100)
    
    # check on each solution
    for (i in 1:nsol) {
      # compute row-wise l2 norm of vhat
      vhat.norm = diag(approximate$projection[[i]])
      Vhat.norm[(k-1) * nsol + i, 1] = k
      Vhat.norm[(k-1) * nsol + i, 2] = i
      Vhat.norm[(k-1) * nsol + i, 3:(p+2)] = vhat.norm
      
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
        decomp = svd(approximate$projection[[i]], nu = d[k], nv = 0)
        vhat = decomp$u
        vhat[which(abs(vhat) < 1e-10)] = 0
        pchat = Xtrain %*% vhat
        
        # fit linear model
        model = lm(Ytrain~pchat)
        gamhat = as.matrix(model$coefficients[-1], ncol = 1)
        betahat = vhat %*% gamhat
      } else {
        # use ymean to predict
        betahat = rep(0, p)
      }
      
      # store outputs
      Betahat[(k-1) * nsol + i, 1] = k
      Betahat[(k-1) * nsol + i, 2] = i
      Betahat[(k-1) * nsol + i, 3:(p+2)] = betahat
      
      # predict on training set
      Ytrain.hat = ymean + Xtrain %*% betahat * ysd
      mse.train[k, i] = mean((Ytrain.hat - (ymean + Ytrain * ysd))^2)
      
      # predict on validation set
      Yvalidate.hat = ymean + Xvalidate %*% betahat * ysd
      mse.validate[k, i] = mean((Yvalidate.hat - Yvalidate)^2)
      
      # predict on test set
      Ytest.hat = ymean + Xtest %*% betahat * ysd
      mse.test[k, i] = mean((Ytest.hat - Ytest)^2)
    }
  }
  
  # return output
  out = list(betahat = Betahat,
             mse.train = mse.train,
             mse.validate= mse.validate,
             mse.test = mse.test,
             vhat.norm = Vhat.norm)
  return(out)
}