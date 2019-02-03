library(glmnet)

# this function compare sparsePC for classification with reletive methods
compare_class = function(Xtrain, Ytrain, Xvalidate, Yvalidate, Xtest, Ytest,
                         d, nsol){
  p = ncol(Xtrain)
  
  # output holder
  sec = rep(NA, 4)
  accuracy = rep(NA, 4)
  ncov = rep(NA, 4)
  
  # fps
  t1 = proc.time()
  sparsePC.mod = sparsePC_class(Xtrain = Xtrain, Ytrain = Ytrain,
                                Xvalidate = Xvalidate, Yvalidate = Yvalidate,
                                Xtest = Xtest, Ytest = Ytest,
                                d = d, nsol = nsol)
  t2 = proc.time()
  sec[1] = t2[3]-t1[3]
  best = which(sparsePC.mod$mse.validate==min(sparsePC.mod$mse.validate[sparsePC.mod$accuracy.validate==max(sparsePC.mod$accuracy.validate)]))
  accuracy[1] = sparsePC.mod$accuracy.test[best]
  ncov[1] = sum(sparsePC.mod$betahat[best,] != 0)
  
  # logistic ridge
  t1 = proc.time()
  lr.mod = glmnet(x = Xtrain, y = Ytrain, alpha = 0, 
                  nlambda = nsol, family = "binomial")
  lr.pre = predict(lr.mod, newx = Xvalidate, type = 'response')
  lr.mse = apply(lr.pre, 2, 
                 function(yhat){mean((yhat - Yvalidate)^2)})
  lr.hat = predict(lr.mod, newx = Xtest, type = 'class',
                   s = lr.mod$lambda[which.min(lr.mse)])
  t2 = proc.time()
  sec[2] = t2[3]-t1[3]
  accuracy[2] = sum(Ytest == lr.hat) / length(Ytest) * 100
  ncov[2] = sum(coef(lr.mod)[-1, which.min(lr.mse)] != 0)
  
  # logistic lasso
  t1 = proc.time()
  ll.mod = glmnet(x = Xtrain, y = Ytrain, alpha = 1, 
                  nlambda = nsol, family = "binomial")
  ll.pre = predict(ll.mod, newx = Xvalidate, type = 'response')
  ll.mse = apply(ll.pre, 2, 
                 function(yhat){mean((yhat - Yvalidate)^2)})
  ll.hat = predict(ll.mod, newx = Xtest, type = 'class',
                   s = ll.mod$lambda[which.min(ll.mse)])
  t2 = proc.time()
  sec[3] = t2[3]-t1[3]
  accuracy[3] = sum(Ytest == ll.hat) / length(Ytest) * 100
  ncov[3] = sum(coef(ll.mod)[-1, which.min(ll.mse)] != 0)
  
  # use majority class
  t1 = proc.time()
  Yhat = rep(NA, length(Ytest))
  Yhat[1:length(Ytest)] = ifelse(sum(Ytrain==1)>=sum(Ytrain==0), 
                                 rep(1, length(Ytest)), rep(0, length(Ytest)))
  t2 = proc.time()
  sec[4] = t2[3]-t1[3]
  accuracy[4] = sum(Ytest == Yhat) / length(Ytest) * 100
  ncov[4] = 0
  
  # return output
  output = list(sec = sec, accuracy = accuracy, ncov = ncov, 
                sparsePC.mod = sparsePC.mod, best = best)
  return(output)
}