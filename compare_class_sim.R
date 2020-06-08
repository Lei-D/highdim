library(glmnet)
source("~/Desktop/Lei/AIMER2/github/sparsePC_class.R")

# compare various models on simulation data
compare_class_sim = function(Xtrain, Ytrain, Xvalidate, Yvalidate, Xtest, Ytest,
                             d, nsol, s){
  p = ncol(Xtrain)
  
  # output holder
  sec = rep(NA, 6)
  accuracy = rep(NA, 6)
  ncov = rep(NA, 6)
  cov.precision = rep(NA, 6)
  cov.recall = rep(NA, 6)
  
  # sparsePC
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
  cov.precision[1] = sum(sparsePC.mod$betahat[best,1:s] != 0) / ncov[1]
  cov.recall[1] = sum(sparsePC.mod$betahat[best,1:s] != 0) / s
  
  # logistic for right features
  t1 = proc.time()
  log.mod = glm(Ytrain~., 
                data = data.frame(Xtrain = Xtrain[,1:s], Ytrain = Ytrain), 
                family = binomial(link = 'logit'))
  t2 = proc.time()
  sec[2] = t2[3]-t1[3]
  log.mod.prd = predict(log.mod, newdata = data.frame(Xtrain = Xtest[ ,1:s]))
  log.mod.hat = ifelse(log.mod.prd >= 0, 1, 0)
  accuracy[2] = sum(Ytest == log.mod.hat) / length(Ytest) * 100
  ncov[2] = s
  cov.precision[2] = 1
  cov.recall[2] = 1
  
  # logistic ridge
  t1 = proc.time()
  lr.mod = glmnet(x = Xtrain, y = Ytrain, alpha = 0, 
                  nlambda = nsol, family = "binomial")
  lr.pre = predict(lr.mod, newx = Xvalidate, type = 'response')
  lr.mse = apply(lr.pre, 2, 
                 function(yhat){mean((yhat - Yvalidate)^2)})
  lr.prd.test = predict(lr.mod, newx = Xtest, type = 'response',
                        s = lr.mod$lambda[which.min(lr.mse)])
  lr.hat = ifelse(lr.prd.test >= 0.5, 1, 0)
  t2 = proc.time()
  sec[3] = t2[3]-t1[3]
  accuracy[3] = sum(Ytest == lr.hat) / length(Ytest) * 100
  ncov[3] = sum(coef(lr.mod)[-1, which.min(lr.mse)] != 0)
  cov.precision[3] = sum(coef(lr.mod)[-1, which.min(lr.mse)][1:s] != 0)/ncov[3]
  cov.recall[3] = 1
  
  # logistic lasso
  t1 = proc.time()
  ll.mod = glmnet(x = Xtrain, y = Ytrain, alpha = 1, 
                  nlambda = nsol, family = "binomial")
  ll.pre = predict(ll.mod, newx = Xvalidate, type = 'response')
  ll.mse = apply(ll.pre, 2, 
                 function(yhat){mean((yhat - Yvalidate)^2)})
  ll.prd.test = predict(ll.mod, newx = Xtest, type = 'response',
                        s = ll.mod$lambda[which.min(ll.mse)])
  ll.hat = ifelse(ll.prd.test >= 0.5, 1, 0)
  t2 = proc.time()
  sec[4] = t2[3]-t1[3]
  accuracy[4] = sum(Ytest == ll.hat) / length(Ytest) * 100
  ncov[4] = sum(coef(ll.mod)[-1, which.min(ll.mse)] != 0)
  cov.precision[4] = sum(coef(ll.mod)[-1, which.min(ll.mse)][1:s] != 0)/ ncov[4]
  cov.recall[4] = sum(coef(ll.mod)[-1, which.min(ll.mse)][1:s] != 0)/ s
  
  # use majority class
  t1 = proc.time()
  Yhat = rep(NA, length(Ytest))
  Yhat[1:length(Ytest)] = ifelse(sum(Ytrain==1)>=sum(Ytrain==0), 
                                 rep(1, length(Ytest)), rep(0, length(Ytest)))
  t2 = proc.time()
  sec[5] = t2[3]-t1[3]
  accuracy[5] = sum(Ytest == Yhat) / length(Ytest) * 100
  ncov[5] = 0
  cov.precision[5] = 0
  cov.recall[5] = 0
  
  # sparsePC without screening
  t1 = proc.time()
  sparsePC.mod = sparsePC_class(Xtrain = Xtrain, Ytrain = Ytrain,
                                Xvalidate = Xvalidate, Yvalidate = Yvalidate,
                                Xtest = Xtest, Ytest = Ytest,
                                d = d, nsol = nsol, screening = FALSE)
  t2 = proc.time()
  sec[6] = t2[3]-t1[3]
  best = which(sparsePC.mod$mse.validate==min(sparsePC.mod$mse.validate[sparsePC.mod$accuracy.validate==max(sparsePC.mod$accuracy.validate)]))
  accuracy[6] = sparsePC.mod$accuracy.test[best]
  ncov[6] = sum(sparsePC.mod$betahat[best,] != 0)
  cov.precision[6] = sum(sparsePC.mod$betahat[best,1:s] != 0) / ncov[6]
  cov.recall[6] = sum(sparsePC.mod$betahat[best,1:s] != 0) / s
  
  # store output
  output = list(sec = sec, accuracy = accuracy, ncov = ncov, 
                cov.precision = cov.precision, 
                cov.recall = cov.recall,
                sparsePC.mod = sparsePC.mod, best = best)
  return(output)
}