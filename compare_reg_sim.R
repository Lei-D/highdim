library(glmnet)

# compare various models on simulation data
compare_reg_sim = function(Xtrain, Ytrain, Xvalidate, Yvalidate, Xtest, Ytest,
                           d, nsol, s){
  p = ncol(Xtrain)
  
  # output holder
  sec = rep(NA, 7)
  mse = rep(NA, 7)
  ncov = rep(NA, 7)
  cov.precision = rep(NA, 7)
  cov.recall = rep(NA, 7)
  
  # fps
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
  cov.precision[1] = sum(sparsePC.mod$betahat[best,1:s] != 0) / ncov[1]
  cov.recall[1] = sum(sparsePC.mod$betahat[best,1:s] != 0) / s
  
  # OLS for right features
  t1 = proc.time()
  OLS.mod = lm(Ytrain ~ .
               , data = data.frame(Xtrain = Xtrain[,1:s], Ytrain = Ytrain))
  t2 = proc.time()
  sec[2] = t2[3]-t1[3]
  mse[2] = mean((predict(OLS.mod, 
                         newdata = data.frame(Xtrain = Xtest[ ,1:s])) - Ytest)^2)
  ncov[2] = s
  cov.precision[2] = 1
  cov.recall[2] = 1
  
  # ridge
  t1 = proc.time()
  ridge.mod = glmnet(x = Xtrain, y = Ytrain, alpha = 0, nlambda = nsol)
  ridge.pre = predict(ridge.mod, newx = Xvalidate, type = 'response')
  ridge.mse = apply(ridge.pre, 2, function(yhat){mean((yhat - Yvalidate)^2)})
  t2 = proc.time()
  sec[3] = t2[3]-t1[3]
  mse[3] = mean((predict(ridge.mod, newx = Xtest, 
                         s = ridge.mod$lambda[which.min(ridge.mse)]) - Ytest)^2)
  ncov[3] = p
  cov.precision[3] = s / p
  cov.recall[3] = 1
  
  # lasso
  t1 = proc.time()
  lasso.mod = glmnet(x = Xtrain, y = Ytrain, alpha = 1, nlambda = nsol)
  lasso.pre = predict(lasso.mod, newx = Xvalidate, type = 'response')
  lasso.mse = apply(lasso.pre, 2, function(yhat){mean((yhat - Yvalidate)^2)})
  t2 = proc.time()
  sec[4] = t2[3]-t1[3]
  mse[4] = mean((predict(lasso.mod, newx = Xtest, 
                         s = lasso.mod$lambda[which.min(lasso.mse)]) - Ytest)^2)
  ncov[4] = sum(coef(lasso.mod)[-1, which.min(lasso.mse)] != 0)
  cov.precision[4] =  sum(coef(lasso.mod)[-1, which.min(lasso.mse)][1:s] != 0)/ ncov[4]
  cov.recall[4] = sum(coef(lasso.mod)[-1, which.min(lasso.mse)][1:s] != 0) / s
  
  # Adaptive Lasso using Ridge coefficients as weights
  t1 = proc.time()
  w.ridge = 1/abs(as.vector(coef(ridge.mod)[-1, which.min(ridge.mse)]))
  adp.lasso1.mod = glmnet(x = Xtrain, y = Ytrain, alpha = 1, 
                          nlambda = nsol, penalty.factor = w.ridge)
  adp.lasso1.pre = predict(adp.lasso1.mod, newx = Xvalidate, type = 'response')
  adp.lasso1.mse = apply(adp.lasso1.pre, 2, function(yhat){mean((yhat - Yvalidate)^2)})
  t2 = proc.time()
  sec[5] = t2[3]-t1[3]
  mse[5] = mean((predict(adp.lasso1.mod, newx = Xtest,
                         s = adp.lasso1.mod$lambda[which.min(adp.lasso1.mse)]) - Ytest)^2)
  ncov[5] = sum(coef(adp.lasso1.mod)[-1, which.min(adp.lasso1.mse)] != 0)
  cov.precision[5] =  sum(coef(adp.lasso1.mod)[-1, which.min(adp.lasso1.mse)][1:s] != 0)/ ncov[5]
  cov.recall[5] = sum(coef(adp.lasso1.mod)[-1, which.min(adp.lasso1.mse)][1:s] != 0) / s
  
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
  sec[6] = t2[3]-t1[3]
  mse[6] = mean((predict(adp.lasso2.mod, newx = Xtest,
                         s = adp.lasso2.mod$lambda[which.min(adp.lasso2.mse)]) - Ytest)^2)
  ncov[6] = sum(coef(adp.lasso2.mod)[-1, which.min(adp.lasso2.mse)] != 0)
  cov.precision[6] =  sum(coef(adp.lasso2.mod)[-1, which.min(adp.lasso2.mse)][1:s] != 0)/ ncov[6]
  cov.recall[6] = sum(coef(adp.lasso2.mod)[-1, which.min(adp.lasso2.mse)][1:s] != 0) / s
  
  # use mean of Ytrain to predict
  sec[7] = 0
  mse[7] = mean((mean(Ytrain) - Ytest)^2)
  ncov[7] = 0
  cov.precision[7] = 0
  cov.recall[7] = 0
  
  output = list(sec = sec, mse = mse, ncov = ncov, 
                cov.precision = cov.precision, 
                cov.recall = cov.recall,
                sparsePC.mod = sparsePC.mod, best = best)
  return(output)
}