library(glmnet)

# this function compare sparsePC with relative methods include lasso, ridge, etc.
compare_reg = function(Xtrain, Ytrain, Xvalidate, Yvalidate, Xtest, Ytest, 
                       d, nsol){
  p = ncol(Xtrain)
  
  # output holder
  sec = rep(NA, 6)
  mse = rep(NA, 6)
  ncov = rep(NA, 6)
  
  # sparsePC
  t1 = proc.time()
  sparsePC.mod = sparsePC_reg(Xtrain = Xtrain, Ytrain = Ytrain,
                              Xvalidate = Xvalidate, Yvalidate = Yvalidate, 
                              Xtest = Xtest, Ytest = Ytest,
                              d = d, nsol = nsol)
  t2 = proc.time()
  sec[1] = t2[3]-t1[3]
  best = which.min(sparsePC.mod$mse.validate)
  mse[1] = sparsePC.mod$mse.test[best]
  ncov[1] = sum(sparsePC.mod$betahat[best, ] != 0)
  
  # ridge
  t1 = proc.time()
  ridge.mod = glmnet(x = Xtrain, y = Ytrain, alpha = 0, nlambda = nsol)
  ridge.pre = predict(ridge.mod, newx = Xvalidate, type = 'response')
  ridge.mse = apply(ridge.pre, 2, function(yhat){mean((yhat - Yvalidate)^2)})
  t2 = proc.time()
  sec[2] = t2[3]-t1[3]
  mse[2] = mean((predict(ridge.mod, newx = Xtest, 
                         s = ridge.mod$lambda[which.min(ridge.mse)]) - Ytest)^2)
  ncov[2] = ncol(Xtrain)
  
  # lasso
  t1 = proc.time()
  lasso.mod = glmnet(x = Xtrain, y = Ytrain, alpha = 1, nlambda = nsol)
  lasso.pre = predict(lasso.mod, newx = Xvalidate, type = 'response')
  lasso.mse = apply(lasso.pre, 2, function(yhat){mean((yhat - Yvalidate)^2)})
  t2 = proc.time()
  sec[3] = t2[3]-t1[3]
  mse[3] = mean((predict(lasso.mod, newx = Xtest, 
                         s = lasso.mod$lambda[which.min(lasso.mse)]) - Ytest)^2)
  ncov[3] = sum(coef(lasso.mod)[-1, which.min(lasso.mse)] != 0)
  
  # Adaptive Lasso using Ridge coefficients as weights
  t1 = proc.time()
  w.ridge = 1/abs(as.vector(coef(ridge.mod)[-1, which.min(ridge.mse)]))
  adp.lasso1.mod = glmnet(x = Xtrain, y = Ytrain, alpha = 1, 
                          nlambda = nsol, penalty.factor = w.ridge)
  adp.lasso1.pre = predict(adp.lasso1.mod, newx = Xvalidate, type = 'response')
  adp.lasso1.mse = apply(adp.lasso1.pre, 2, function(yhat){mean((yhat - Yvalidate)^2)})
  t2 = proc.time()
  sec[4] = t2[3]-t1[3]
  mse[4] = mean((predict(adp.lasso1.mod, newx = Xtest,
                         s = adp.lasso1.mod$lambda[which.min(adp.lasso1.mse)]) - Ytest)^2)
  ncov[4] = sum(coef(adp.lasso1.mod)[-1, which.min(adp.lasso1.mse)] != 0)
  
  # Adaptive Lasso using Lasso coefficients as weights
  t1 = proc.time()
  w.lasso = 1/abs(as.vector(coef(lasso.mod)[-1, which.min(lasso.mse)]))
  for (i in 1:p) {
    if (w.lasso[i] == Inf) {
      w.lasso[i] = 1e+9
    }
  }
  adp.lasso2.mod = glmnet(x = Xtrain, y = Ytrain, alpha = 1, 
                          nlambda = nsol, penalty.factor = w.lasso)
  adp.lasso2.pre = predict(adp.lasso2.mod, newx = Xvalidate, type = 'response')
  adp.lasso2.mse = apply(adp.lasso2.pre, 2, function(yhat){mean((yhat - Yvalidate)^2)})
  t2 = proc.time()
  sec[5] = t2[3]-t1[3]
  mse[5] = mean((predict(adp.lasso2.mod, newx = Xtest,
                         s = adp.lasso2.mod$lambda[which.min(adp.lasso2.mse)]) - Ytest)^2)
  ncov[5] = sum(coef(adp.lasso2.mod)[-1, which.min(adp.lasso2.mse)] != 0)
  
  # use mean of Ytrain to predict
  sec[6] = 0
  mse[6] = mean((mean(Ytrain) - Ytest)^2)
  ncov[6] = 0
  
  output = list(sec = sec, mse = mse, ncov = ncov,
                sparsePC.mod = sparsePC.mod, best = best)
  return(output)
}