# install packages from github for first time user
# library(devtools)
# install_github("Lei-D/irlba", ref="lei")
# install_github("Lei-D/fps", ref="fps_irlba")
library(irlba)
library(fps)
source("~/Desktop/Lei/AIMER2/github/find_t.R")

sparsePC.reg.real = function(Xtrain, Ytrain, 
                             Xvalidate = NULL, Yvalidate = NULL, 
                             d, lambda, maxnvar, screening = TRUE){
  n = nrow(Xtrain)
  p = ncol(Xtrain)
  nsol = length(lambda)
  
  # center input data X
  xmeans = colMeans(Xtrain)
  xsd = apply(Xtrain, 2, sd)
  Xtrain = scale(Xtrain, center = xmeans, scale = xsd)
  
  # center input data Y
  ymean = mean(Ytrain)
  ysd = sd(Ytrain)
  Ytrain = scale(Ytrain, center = ymean, scale = ysd)
  
  # output holder
  mse.train = matrix(NA, nrow=length(d), ncol = nsol)
  Betahat = matrix(NA, nrow = nsol*length(d), ncol = p+2)
  Vhat.norm = matrix(NA, nrow = nsol*length(d), ncol = p+2)
  
  if (is.null(Xvalidate) == F) {
    Xvalidate = scale(Xvalidate, center = xmeans, scale = xsd)
    mse.validate = matrix(NA, nrow=length(d), ncol = nsol)
  } else {
    mse.validate = NULL
  }
  
  # generate sample covariance matrix
  S = crossprod(Xtrain) / n
  
  for(k in 1:length(d)){
    # run fps
    approximate = fps(S, ndim = d[k], ncomp = d[k], 
                      lambda = lambda,
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
      if (is.null(Yvalidate) == F) {
        Yvalidate.hat = ymean + Xvalidate %*% betahat * ysd
        mse.validate[k, i] = mean((Yvalidate.hat - Yvalidate)^2)
      }
    }
  }
  
  # return output
  out = list(betahat = Betahat,
             mse.train = mse.train,
             mse.validate = mse.validate,
             vhat.norm = Vhat.norm)
  return(out)
}


sparsePC.reg.real.CV = function(X, Y, d, lambda, Kfold, 
                                maxnvar, screening = TRUE){
  n = nrow(X)
  p = ncol(X)
  nsol = length(lambda)
  
  # center input data X
  xmeans = colMeans(X)
  xsd = apply(X, 2, sd)
  X = scale(X, center = xmeans, scale = xsd)
  
  # center input data Y
  ymean = mean(Y)
  ysd = sd(Y)
  Y = scale(Y, center = ymean, scale = ysd)
  
  # split data for CV
  dataCV = MykfoldCV(x = X, y = Y, k = Kfold)  
  
  # output holder
  mse.validate = array(NA, dim = c(length(d), nsol, Kfold))

  for(k in 1:Kfold){
    # prepare training data and testing data
    Xvalidate = dataCV$dataset[[k]]$testing.x
    Yvalidate = dataCV$dataset[[k]]$testing.y
    Xtrain = dataCV$dataset[[k]]$training.x
    Ytrain = dataCV$dataset[[k]]$training.y
    
    fit = sparsePC.reg.real(Xtrain = Xtrain, Ytrain = Ytrain,
                            Xvalidate = Xvalidate, Yvalidate = Yvalidate,
                            d = d, 
                            lambda = lambda, 
                            maxnvar = p, 
                            screening = TRUE)
    mse.validate[,,k] = fit$mse.validate
  }
  mse = apply(mse.validate, c(1,2), mean)
  best = which(mse == min(mse), arr.ind = TRUE)[1, ]
  
  # fit the model with whole data
  final = sparsePC.reg.real(Xtrain = X, Ytrain = Y,
                          d = d[best[1]], 
                          lambda = lambda, 
                          maxnvar = p, 
                          screening = TRUE)
  
  # return output
  out = list(mse.validate = mse.validate,
             mse = mse,
             best = best,
             d = d[best[1]],
             lambda = lambda[best[2]],
             betahat = final$betahat[best[2], 3:(p+2)])
  return(out)
}


MykfoldCV = function(x, y, k){
  ## Purpose: prepare training set and testing set for kfold cross-validation
  ## Inputs: an n*p matrix x as data matrix;
  ##         an n*1 vector y as response;
  ##         k is the number of folds
  ## Outputs: a list of 2 components, first is the list of dataset including k lists, 
  ##          each includes sets of training.x, training.y, testing.x, testing.y;
  ##          second is the folds list which indicates the index for each fold.
  
  # randomly select index for partition of rows of x
  n = nrow(x)
  folds = vector("list", k)
  breaks = round(seq(from = 1, to = (n + 1), length = (k + 1)))
  cv.order = sample(1:n)
  dataset = vector("list", k)
  for(i in 1:k){
    # prepare index for the ith fold
    folds[[i]] = cv.order[(breaks[i]):(breaks[i + 1] - 1)]
    # generate training set and testing set for ith fold
    testing.x = x[folds[[i]], ]
    testing.y = y[folds[[i]]]
    training.x = x[-folds[[i]], ]
    training.y = y[-folds[[i]]]
    dataset[[i]] = list(testing.x = testing.x, testing.y = testing.y, 
                         training.x = training.x, training.y=training.y)
  }
  out = list(dataset=dataset, folds=folds)
  return(out)
}
