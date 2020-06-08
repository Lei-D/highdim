# install packages from github for first time user
# library(devtools)
# install_github("Lei-D/irlba", ref="lei")
# install_github("Lei-D/fps", ref="fps_irlba")
library(irlba)
library(fps)

# data generating function for K=M=3 case (cs has advantage)
factorModelSim2 = function(n, p, p1, lambdas, theta1, theta2, SNRx, SNRy){
    V = matrix(c(rep(1, 3*p1), rep(c(1, -1/4, -3/4), times=c(p1, p1, p1)), 
                 rep(c(1/10, -7/20, 1/4), times=c(p1, p1, p1))), ncol = 3)
    V = rbind(V,matrix(0, p-3*p1, 3)) # still sparse, the remaining rows have norm 0
    V = V %*% solve(diag(c(sqrt(sum(V[, 1]^2))
                           , sqrt(sum(V[, 2]^2))
                           , sqrt(sum(V[, 3]^2))))) # normalize
    sigma0 = sqrt(sum(lambdas))/sqrt(p)/SNRx
    Lambda = diag(lambdas) 
    D = diag(lambdas + sigma0^2)
    W = V %*% sqrt(Lambda)
    w = W[2*p1+1, ]
    theta3 = (- w[1]*theta1 - w[2]*theta2) / w[3]
    theta = c(theta1, theta2, theta3)
    SigXY = W %*% theta
    beta = W %*% solve(D) %*% theta
    U = matrix(rnorm(3 * n), ncol = 3)
    sigma1 = as.numeric(sqrt((crossprod(sqrt(Lambda)%*%t(V)%*%beta) + (sigma0^2) * crossprod(beta))/n)/SNRy)
    
    # add noise and generate traing set
    Xtrain = U %*% sqrt(Lambda) %*% t(V) + sigma0 * matrix(rnorm(n*p), n, p)
    Ytrain = U %*% theta + sigma1 * rnorm(n)
    
    # add noise and generate validation set
    Xvalidate = U %*% sqrt(Lambda) %*% t(V) + sigma0 * matrix(rnorm(n*p), n, p)
    Yvalidate = U %*% theta + sigma1 * rnorm(n)
    
    # add noise and generate test set
    Xtest = U %*% sqrt(Lambda) %*% t(V) + sigma0 * matrix(rnorm(n*p), n, p)
    Ytest = U %*% theta + sigma1 * rnorm(n)
    
    ret = list(Xtrain = Xtrain,Ytrain = Ytrain
               , Xvalidate = Xvalidate, Yvalidate = Yvalidate
               , Xtest = Xtest, Ytest = Ytest
               , beta = beta, theta = theta, SigXY = SigXY
               , U = U, V = V, Lambda = Lambda)
    return(ret)
}


ROC = function(Xtrain, d, nsol){
    n = nrow(Xtrain)
    p = ncol(Xtrain)
    
    # center input data X
    xmeans = colMeans(Xtrain)
    xsd = apply(Xtrain, 2, sd)
    Xtrain = scale(Xtrain, center = xmeans, scale = xsd)
    
    # output holder
    TPR = matrix(NA, nrow = 1000, ncol = nsol)
    FPR = matrix(NA, nrow = 1000, ncol = nsol)
    
    # generate sample covariance matrix
    S = crossprod(Xtrain) / n
    
    # run fps
    approximate = fps(S, ndim = d, ncomp = d, 
                      lambda = seq(from = 1, to = 0, length.out = nsol),
                      maxnvar = p, maxiter = 100)
    
    # check on each solution
    for (i in 1:nsol) {
        # decompose projection matrix
        decomp = svd(approximate$projection[[i]], nu = d, nv = 0)
        vhat = decomp$u
        vhat[which(abs(vhat) < 1e-10)] = 0
        l = apply(vhat, 1, function(x)sum(x^2))
        threshold = sort(l, decreasing = TRUE)
        for (j in 1:1000) {
            TPR[j, i] = sum(l[1:15] >= threshold[j]) / 15
            FPR[j, i] = sum(l[16:1000] >= threshold[j]) / 985
        }
    }
    out = list(TPR = TPR, FPR = FPR)
    return(out)
}

# generate simulation data
set.seed(19812)
data1 = factorModelSim2(n = 100, p = 1000, p1 = 5, lambdas = c(10, 5, 1)
                        , theta1 = 1, theta2 = 1
                        , SNRx = 5, SNRy = 5)
sim1.ROC = ROC(Xtrain = data1$Xtrain, d = 3, nsol = 10)
df.SuffPCR = data.frame(TPR = sim1.ROC$TPR[, 9], 
                        FPR = sim1.ROC$FPR[, 9], 
                        method = "SuffPCR")

# Lasso
library(glmnet)
lasso.mod = glmnet(x = data1$Xtrain, y = data1$Ytrain, 
                   alpha = 1, 
                   nlambda = 200,
                   thresh = 1e-12,
                   lambda.min.ratio = 1e-8)
# output holder
nlam.lasso = ncol(coef(lasso.mod))
lasso.TPR = rep(NA, nlam.lasso+1)
lasso.FPR = rep(NA, nlam.lasso+1)

for (i in 1:nlam.lasso) {
    lasso.TPR[i] = sum(coef(lasso.mod)[2:16, i] != 0)/15
    lasso.FPR[i] = sum(coef(lasso.mod)[17:1001, i] != 0)/985
}
lasso.TPR[nlam.lasso+1] = 1
lasso.FPR[nlam.lasso+1] = 1
df.lasso = data.frame(TPR = lasso.TPR, FPR = lasso.FPR, method = 'Lasso')

# ElasticNet
ElasticNet.mod = glmnet(x = data1$Xtrain, y = data1$Ytrain, 
                   alpha = 0.5, 
                   nlambda = 200,
                   thresh = 1e-12,
                   lambda.min.ratio = 1e-8)
# output holder
nlam.elnet = ncol(coef(ElasticNet.mod))
ElasticNet.TPR = rep(NA, nlam.elnet+1)
ElasticNet.FPR = rep(NA, nlam.elnet+1)

for (i in 1:nlam.elnet) {
    ElasticNet.TPR[i] = sum(coef(ElasticNet.mod)[2:16, i] != 0)/15
    ElasticNet.FPR[i] = sum(coef(ElasticNet.mod)[17:1001, i] != 0)/985
}

ElasticNet.FPR[nlam.elnet+1] = 1
ElasticNet.TPR[nlam.elnet+1] = 1
df.elnet = data.frame(TPR = ElasticNet.TPR, FPR = ElasticNet.FPR, method = "ElasticNet")

# SPC
marginalRegressionT <- function (x, y) {
    ## Purpose: gets t-statistics via marginal regression
    ## Inputs: an n x p matrix x
    ##             an n vector y
    ## Outputs: a p vector of t-statistics for the univariate (marginal regressions)
    ##                of each column of x on y
    y <- as.vector(y)
    n <- length(y)
    cx <- scale(x, scale=FALSE) # centered only
    cy <- y - mean(y)
    sxx <- colSums(cx^2)
    sxy <- crossprod(cx , cy) # t(x) %*% y but faster
    syy <- sum(cy^2)
    numer <- sxy/sxx
    std <- sqrt((syy/sxx - numer^2)/(n - 2))
    return(numer / std)
}

tvalue = rep(NA, 1000)
for(i in 1:1000){
    tvalue[i] = marginalRegressionT(x = data1$Xtrain[, i], y = data1$Ytrain)
}
tvalue.sort = sort(tvalue, decreasing = TRUE)

# output holder
SPC.TPR = rep(NA, 1000)
SPC.FPR = rep(NA, 1000)

for(j in 1:1000){
    SPC.TPR[j] = sum(tvalue[1:15] >= tvalue.sort[j]) / 15
    SPC.FPR[j] = sum(tvalue[16:1000] >= tvalue.sort[j]) / 985
}
df.SPC = data.frame(TPR = SPC.TPR, FPR = SPC.FPR, method = "SPC")

# combine output
library(tidyverse)
dat = bind_rows(df.SuffPCR, df.lasso, df.elnet, df.SPC)
dat$method = factor(dat$method, 
                    levels = c("SuffPCR", "Lasso", "ElasticNet", "SPC"))
save(dat, file = "~/Desktop/Lei/AIMER2/rebuttal/ROC/ROC.RData")



load("~/Desktop/Lei/AIMER2/rebuttal/ROC/ROC.RData")
ggplot(data = dat, aes(x = FPR, y = TPR, colour = method)) + 
  geom_line(size = 2) + 
  coord_cartesian(xlim=c(0,0.05)) +
  scale_color_manual(values = viridis::viridis_pal()(10)[c(10, 8, 6, 5)]) + 
  theme(text = element_text(size=18,family='Palatino'))








