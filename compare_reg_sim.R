library(glmnet)
library(pcLasso)
library(dimreduce)
library(tidyverse)
source("~/Desktop/Lei/AIMER2/github/AIMER.R")
source("~/Desktop/Lei/AIMER2/github/sparsePC_reg.R")
source("~/Desktop/Lei/AIMER2/github/find_t.R")

# compare various models on simulation data
compare_reg_sim = function(Xtrain, Ytrain, Xvalidate, Yvalidate, Xtest, Ytest,
                           d, nsol, s, d.aimer){
  n = nrow(Xtrain)
  p = ncol(Xtrain)
  
  # output holder
  sec = rep(NA, 13)
  mse = rep(NA, 13)
  ncov = rep(NA, 13)
  cov.precision = rep(NA, 13)
  cov.recall = rep(NA, 13)
  betahat = data.frame(Method = character(),
                       Index = integer(),
                       Value = numeric(),
                       stringsAsFactors = FALSE)
  # sparsePC
  t1 = proc.time()
  sparsePC.mod = sparsePC_reg(Xtrain = Xtrain, Ytrain = Ytrain,
                              Xvalidate = Xvalidate, Yvalidate = Yvalidate,
                              Xtest = Xtest, Ytest = Ytest,
                              d = d, nsol = nsol, maxnvar = p)
  t2 = proc.time()
  sec[1] = t2[3]-t1[3]
  sparsePC.best = which.min(sparsePC.mod$mse.validate)
  mse[1] = sparsePC.mod$mse.test[sparsePC.best]
  ncov[1] = sum(sparsePC.mod$betahat[sparsePC.best, 3:(p+2)] != 0)
  cov.precision[1] = ifelse(ncov[1] == 0, 0,
                            sum(sparsePC.mod$betahat[sparsePC.best,3:(s+2)] != 0) / ncov[1])
  cov.recall[1] = sum(sparsePC.mod$betahat[sparsePC.best,3:(s+2)] != 0) / s
  # compute ROC
  sparsePC.TPR = rep(NA, p)
  sparsePC.FPR = rep(NA, p)
  l = sparsePC.mod$vhat.norm[sparsePC.best,3:(p+2)]
  threshold = sort(l, decreasing = TRUE)
  for (j in 1:p) {
    sparsePC.TPR[j] = sum(l[1:s] >= threshold[j]) / s
    sparsePC.FPR[j] = sum(l[(s+1):p] >= threshold[j]) / (p-s)
  }
  ROC = data.frame(Index = 1:(length(sparsePC.TPR) + 1),
                   TPR = c(0, sparsePC.TPR), 
                   FPR = c(0,sparsePC.FPR), 
                   Method = "SuffPCR",
                   stringsAsFactors = FALSE)
  betahat = bind_rows(betahat, 
                      data.frame(Method = "SuffPCR",
                                 Index = 1:s,
                                 Value = sparsePC.mod$betahat[sparsePC.best, 3:(s+2)],
                                 stringsAsFactors = FALSE))
  
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
  betahat = bind_rows(betahat, 
                      data.frame(Method = "Oracle",
                                 Index = 1:s,
                                 Value = coef(OLS.mod)[-1],
                                 stringsAsFactors = FALSE))
  
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
  betahat = bind_rows(betahat, 
                      data.frame(Method = "Ridge",
                                 Index = 1:s,
                                 Value = ridge.mod$beta[1:s,which.min(ridge.mse)],
                                 stringsAsFactors = FALSE))
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
  cov.precision[4] = ifelse(ncov[4] == 0, 0,
                            sum(coef(lasso.mod)[-1, which.min(lasso.mse)][1:s] != 0)/ ncov[4])
  cov.recall[4] = sum(coef(lasso.mod)[-1, which.min(lasso.mse)][1:s] != 0) / s
  # compute ROC
  nlam.lasso = ncol(coef(lasso.mod))
  lasso.TPR = rep(NA, nlam.lasso + 1)
  lasso.FPR = rep(NA, nlam.lasso + 1)
  for (i in 1:nlam.lasso) {
    lasso.TPR[i] = sum(coef(lasso.mod)[2:(s+1), i] != 0)/s
    lasso.FPR[i] = sum(coef(lasso.mod)[(s+2):(p+1), i] != 0)/(p-s)
  }
  lasso.TPR[nlam.lasso + 1] = 1
  lasso.FPR[nlam.lasso + 1] = 1
  ROC = bind_rows(ROC, data.frame(Index = 1:(length(lasso.TPR) + 1),
                                  TPR = c(0, lasso.TPR), 
                                  FPR = c(0, lasso.FPR), 
                                  Method = 'Lasso',
                                  stringsAsFactors = FALSE))
  betahat = bind_rows(betahat, 
                      data.frame(Method = "Lasso",
                                 Index = 1:s,
                                 Value = lasso.mod$beta[1:s,which.min(lasso.mse)],
                                 stringsAsFactors = FALSE))
  
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
  cov.precision[5] = ifelse(ncov[5] == 0, 0,
                            sum(coef(adp.lasso1.mod)[-1, which.min(adp.lasso1.mse)][1:s] != 0)/ ncov[5])
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
  cov.precision[6] =  ifelse(ncov[6] == 0, 0,
                             sum(coef(adp.lasso2.mod)[-1, which.min(adp.lasso2.mse)][1:s] != 0)/ ncov[6])
  cov.recall[6] = sum(coef(adp.lasso2.mod)[-1, which.min(adp.lasso2.mse)][1:s] != 0) / s
  
  # use mean of Ytrain to predict
  sec[7] = 0
  mse[7] = mean((mean(Ytrain) - Ytest)^2)
  ncov[7] = 0
  cov.precision[7] = 0
  cov.recall[7] = 0
  
  # Elastic net
  t1 = proc.time()
  mse.min = rep(NA, nsol)
  for (j in 1:nsol) {
    ElasticNet.mod = glmnet(x = Xtrain, y = Ytrain
                            , alpha = seq(0.1, 0.9, length.out = nsol)[j]
                            , nlambda = nsol)
    ElasticNet.pre = predict(ElasticNet.mod, newx = Xvalidate
                             , type = 'response')
    ElasticNet.mse = apply(ElasticNet.pre, 2
                           , function(yhat){mean((yhat - Yvalidate)^2)})
    mse.min[j] = min(ElasticNet.mse)
  }
  ElasticNet.mod = glmnet(x = Xtrain, y = Ytrain
                          , alpha = seq(0.1, 0.9, length.out = nsol)[which.min(mse.min)]
                          , nlambda = nsol)
  ElasticNet.pre = predict(ElasticNet.mod, newx = Xvalidate
                           , type = 'response')
  ElasticNet.mse = apply(ElasticNet.pre, 2
                         , function(yhat){mean((yhat - Yvalidate)^2)})
  t2 = proc.time()
  sec[8] = t2[3]-t1[3]
  mse[8] = mean((predict(ElasticNet.mod, newx = Xtest,
                         s = ElasticNet.mod$lambda[which.min(ElasticNet.mse)]) - Ytest)^2)
  ncov[8] = sum(coef(ElasticNet.mod)[-1, which.min(ElasticNet.mse)] != 0)
  cov.precision[8] =  ifelse(ncov[8] == 0, 0,
                             sum(coef(ElasticNet.mod)[-1, which.min(ElasticNet.mse)][1:s] != 0)/ ncov[8])
  cov.recall[8] = sum(coef(ElasticNet.mod)[-1, which.min(ElasticNet.mse)][1:s] != 0) / s
  # compute ROC
  nlam.elnet = ncol(coef(ElasticNet.mod))
  ElasticNet.TPR = rep(NA, nlam.elnet + 1)
  ElasticNet.FPR = rep(NA, nlam.elnet + 1)
  for (i in 1:nlam.elnet) {
    ElasticNet.TPR[i] = sum(coef(ElasticNet.mod)[2:(s+1), i] != 0)/s
    ElasticNet.FPR[i] = sum(coef(ElasticNet.mod)[(s+2):(p+1), i] != 0)/(p-s)
  }
  
  ElasticNet.FPR[nlam.elnet + 1] = 1
  ElasticNet.TPR[nlam.elnet + 1] = 1
  ROC = bind_rows(ROC, data.frame(Index = 1:(length(ElasticNet.TPR)+1),
                                  TPR = c(0, ElasticNet.TPR), 
                                  FPR = c(0, ElasticNet.FPR), 
                                  Method = "ElasticNet",
                                  stringsAsFactors = FALSE))
  betahat = bind_rows(betahat, 
                      data.frame(Method = "ElasticNet",
                                 Index = 1:s,
                                 Value = coef(ElasticNet.mod)[-1, which.min(ElasticNet.mse)][1:s],
                                 stringsAsFactors = FALSE))
  
  # SPC
  t1 = proc.time()
  nCovs = floor(seq(from = d+2, to = p/20, length.out = nsol))
  mse.validate = rep(NA, nsol)
  SPC.TPR = rep(NA, nsol + 1)
  SPC.FPR = rep(NA, nsol + 1)
  for (i in 1:nsol) {
    fit = approximatPCR(x = Xtrain, y = Ytrain, ncomp = d, 
                        nCov = nCovs[i], PCAmethod = "SPC")
    Yvalidate.pre = predict(fit, newx = Xvalidate, PCAmethod = "SPC")
    mse.validate[i] = mean((Yvalidate - Yvalidate.pre)^2)
    SPC.TPR[i] = sum(fit$keep[1:s])/s
    SPC.FPR[i] = sum(fit$keep[(s+1):p])/(p-s)
  }
  SPC.best = which.min(mse.validate)
  fit.best = approximatPCR(x = Xtrain, y = Ytrain, ncomp = d, 
                           nCov = nCovs[SPC.best], PCAmethod = "SPC")
  Ytest.pre = predict(fit.best, newx = Xtest, PCAmethod = "SPC")
  t2 = proc.time()
  sec[9] = t2[3]-t1[3]
  mse[9] = mean((Ytest - Ytest.pre)^2)
  ncov[9] = sum(fit.best$keep)
  cov.precision[9] = ifelse(ncov[9] == 0, 0,
                            sum(fit.best$keep[1:s]) / ncov[9])
  cov.recall[9] = sum(fit.best$keep[1:s])/s 
  SPC.TPR[nsol + 1] = 1
  SPC.FPR[nsol + 1] = 1
  ROC = bind_rows(ROC, data.frame(Index = 1:(length(SPC.TPR) + 1),
                                  TPR = c(0, SPC.TPR), 
                                  FPR = c(0, SPC.FPR), 
                                  Method = "SPC",
                                  stringsAsFactors = FALSE))
  SPC.est = rep(0, p)
  SPC.est[fit.best$keep] = fit.best$betahat[,1]
  betahat = bind_rows(betahat, 
                      data.frame(Method = "SPC",
                                 Index = 1:s,
                                 Value = SPC.est[1:s],
                                 stringsAsFactors = FALSE))
  
  # AIMER
  t1 = proc.time()
  # mse.min = matrix(NA, nrow = nsol, ncol = nsol)
  # for (i in 1:nsol) {
  #   for (j in 1:nsol) {
  #     AIMER.mod = approximatPCR(x = Xtrain, y = Ytrain, 
  #                               ncomp = d.aimer, nCov = nCovs[i], 
  #                               nCov.select = nCovs[j],
  #                               PCAmethod = "cs.select")
  #     AIMER.pre = predict(AIMER.mod, newx = Xvalidate, PCAmethod = "cs.select")
  #     mse.min[i, j] = mean((AIMER.pre - Yvalidate)^2)
  #   }
  # }
  # best = which(mse.min == min(mse.min), arr.ind = TRUE)
  # AIMER.mod = approximatPCR(x = Xtrain, y = Ytrain, 
  #                           ncomp = d.aimer, nCov = nCovs[best[1]],
  #                           nCov.select = nCovs[best[2]],
  #                           PCAmethod = "cs.select")
  mse.min = rep(NA, nsol)
  AIMER.TPR = rep(NA, nsol + 1)
  AIMER.FPR = rep(NA, nsol + 1)
  for (i in 1:nsol) {
    AIMER.mod = approximatPCR(x = Xtrain, y = Ytrain, 
                              ncomp = d.aimer, nCov = 10, 
                              nCov.select = nCovs[i],
                              PCAmethod = "cs.select")
    AIMER.pre = predict(AIMER.mod, newx = Xvalidate, PCAmethod = "cs.select")
    mse.min[i] = mean((AIMER.pre - Yvalidate)^2)
    AIMER.TPR[i] = sum(AIMER.mod$b.keep[1:s])/s
    AIMER.FPR[i] = sum(AIMER.mod$b.keep[(s+1):p])/(p-s)
  }
  best = which.min(mse.min)
  AIMER.mod = approximatPCR(x = Xtrain, y = Ytrain, 
                            ncomp = d.aimer, nCov = 10,
                            nCov.select = nCovs[best],
                            PCAmethod = "cs.select")
  AIMER.pre = predict(AIMER.mod, newx = Xtest, PCAmethod = "cs.select")
  t2 = proc.time()
  sec[10] = t2[3]-t1[3]
  mse[10] = mean((AIMER.pre - Ytest)^2)
  ncov[10] = sum(AIMER.mod$b.keep)
  cov.precision[10] =  ifelse(ncov[10] == 0, 0,
                              sum(AIMER.mod$b.keep[1:s]) / ncov[10])
  cov.recall[10] = sum(AIMER.mod$b.keep[1:s]) / s
  AIMER.TPR[nsol + 1] = 1
  AIMER.FPR[nsol + 1] = 1
  ROC = bind_rows(ROC, data.frame(Index = 1:(length(AIMER.TPR) + 1),
                                  TPR = c(0, AIMER.TPR), 
                                  FPR = c(0, AIMER.FPR), 
                                  Method = "AIMER",
                                  stringsAsFactors = FALSE))  
  AIMER.est = rep(0, p)
  AIMER.est[AIMER.mod$b.keep] = AIMER.mod$betahat
  betahat = bind_rows(betahat, 
                      data.frame(Method = "AIMER",
                                 Index = 1:s,
                                 Value = AIMER.est[1:s],
                                 stringsAsFactors = FALSE))
  
  # pcLasso
  t1 = proc.time()
  mse.min = rep(NA, 10)
  for (j in 1:nsol) {
    pcLasso.mod = pcLasso(x = Xtrain, y = Ytrain
                          , family = "gaussian"
                          , ratio = seq(0.5, 0.9, length.out = nsol)[j]
                          , nlam = 10)
    pcLasso.pre = predict(pcLasso.mod, xnew = Xvalidate)
    pcLasso.mse = apply(pcLasso.pre, 2
                        , function(yhat){mean((yhat - Yvalidate)^2)})
    mse.min[j] = min(pcLasso.mse)
  }
  pcLasso.mod = pcLasso(x = Xtrain, y = Ytrain
                        , family = "gaussian"
                        , ratio = seq(0.5, 0.9, length.out = nsol)[which.min(mse.min)]
                        , nlam = 10)
  pcLasso.pre = predict(pcLasso.mod, xnew = Xvalidate)
  pcLasso.mse = apply(pcLasso.pre, 2
                      , function(yhat){mean((yhat - Yvalidate)^2)})
  t2 = proc.time()
  sec[11] = t2[3]-t1[3]
  mse[11] = mean((predict(pcLasso.mod, xnew = Xtest)[,which.min(pcLasso.mse)] - Ytest)^2)
  ncov[11] = sum(pcLasso.mod$beta[, which.min(pcLasso.mse)] != 0)
  cov.precision[11] = ifelse(ncov[11] == 0, 0,
                             sum(pcLasso.mod$beta[1:s, which.min(pcLasso.mse)] != 0) / ncov[11])
  cov.recall[11] = sum(pcLasso.mod$beta[, which.min(pcLasso.mse)] != 0)/s 
  
  # ISPCA
  t1 = proc.time()
  # find estimated pc
  ISPCA.mod = try(ispca(Xtrain, Ytrain, nctot = d, verbose = F, alpha = 0.1), silent = T) 
  t2 = proc.time()
  sec[12] = t2[3]-t1[3]
  if (class(ISPCA.mod) != "try-error") {
    # fit linear model on pc
    lm.mod = lm(Ytrain~ISPCA.mod$z)
    pc.pred = predict(ISPCA.mod, Xtest)
    lm.pred = predict(lm.mod, as.data.frame(pc.pred))
    mse[12] = mean((lm.pred - Ytest)^2)
    V.row.norm = apply(ISPCA.mod$w, 1, function(x)sum(x^2))
    ncov[12] = sum(V.row.norm != 0)
    cov.precision[12] = ifelse(ncov[12] == 0, 0,
                               sum(V.row.norm[1:s] != 0) / ncov[12])
    cov.recall[12] = sum(V.row.norm[1:s] != 0) / s
  }
  betahat = bind_rows(betahat, 
                      data.frame(Method = "ISPCA",
                                 Index = 1:s,
                                 Value = rep(NA, s),
                                 stringsAsFactors = FALSE))
  
  # sparsePC without screening
  t1 = proc.time()
  sparsePC.mod = sparsePC_reg(Xtrain = Xtrain, Ytrain = Ytrain,
                              Xvalidate = Xvalidate, Yvalidate = Yvalidate,
                              Xtest = Xtest, Ytest = Ytest,
                              d = d, nsol = nsol, maxnvar = p,
                              screening = FALSE)
  t2 = proc.time()
  sec[13] = t2[3]-t1[3]
  best = which.min(sparsePC.mod$mse.validate)
  mse[13] = sparsePC.mod$mse.test[best]
  ncov[13] = sum(sparsePC.mod$betahat[best, ] != 0)
  cov.precision[13] = ifelse(ncov[13] == 0, 0,
                             sum(sparsePC.mod$betahat[best,1:s] != 0) / ncov[13])
  cov.recall[13] = sum(sparsePC.mod$betahat[best,1:s] != 0) / s
  betahat = bind_rows(betahat, 
                      data.frame(Method = "FPS",
                                 Index = 1:s,
                                 Value = sparsePC.mod$betahat[best, 1:s],
                                 stringsAsFactors = FALSE))
  # store output
  output = list(sec = sec, mse = mse, ncov = ncov, 
                cov.precision = cov.precision, 
                cov.recall = cov.recall,
                sparsePC.mod = sparsePC.mod,
                sparsePC.best = sparsePC.best,
                ROC = ROC, 
                betahat = betahat)
  return(output)
}