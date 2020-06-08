## To be done when building a package:
## 1. MySvd: prompt error when min(dim(x)) < ncomponent (in cross-validation, 
##           we need to make sure nrow of training x is larger than ncomponent)
## 2. estimateSuperPCA, findThreshold: check if ncomp < min(nCovs)
## 3. approximatPCR: check if ncomp < nCov
## 4. might only create functions for PCAmethod "cs.select"
## 5. may need to change the order of betahat into the original order 
##    in the output of approximatPCR and change the predict function accordingly.


## Main changes comparing to v8:
## 1. keep approximatPCR as in v7, which only work for fixed nCov and ncomp; 
##    change findThreshold to do cross-validation over nCov, ncomp, and nCov.select
## 2. in predict SPC.lasso, specify s="lambda.min"


## Main changes comparing to v7:
## 1. add cross validation for number of components in superPCA and findThreshold


## Main changes comparing to v6:
## 1. add argument nCovs(meaning number of Covariates) in estimateSuperPCA, 
##    and make it to be the cross-validated instead of values of threshold.
## 2. delete argument thresholds, UQ.tStats, LQ.tStats, UQ.betahat, LQ.betahat
## 3. delete argument svd="MySvd" 
## 4. delete argument type=c("GCV", "looCV", "kfoldCV"). Delete the result about "GCV" and "loocV".
##    Right now we use kfold cross-validation by default.
## 5. delete some unuseful functions
## 6. delete code for PCAmethod = "cs.keep"
## 7. change the name of the methods: change from "regular" to "SPC", 
##    from "regular.lasso" to "SPC.lasso", from "cs.all" to "cs", 
##    from "cs.all.select" to "cs.select"
## 8. In MySvd function, if ncomponent <= min(nrow(x), ncol(x)) * .5, use irlba, 
##    otherwise use svd.


## Main changes comparing to v5dm:
## 1. Create new function createThresholds
## 2. Use recursive algorithm in function superPCA for method regular.lasso, cs.all.select, 
##    and cs.keep
## 3. Include intercept in regular.lasso


## Main changes comparing to v5:
## 1. removed the "power inflation" from the estimation of d in CS
## 2. Added irlba on the rest of MySVD
## 3. Commented out spatial stuff.
## 4. Add PCAmethod regular.lasso
## 5. Add PCAmethod cs.all.select


## Main changes comparing to v4: 
## 1. Adjust column-sampling method in superPCA and predict.supervisedPCA function.
## 2. In predict.supervisedPCA function, use the means we derived on the training data to 
##    center new X.
## 3. Change the generating process for potential thresholds: define the maximum threshold as 
##    the 95% quantile of abs(tStats) and the minimum as the 5% quantile of abs(tStats).
## 4. Add a new covariance matrix type: BandDiagonal.
## 5. Use irlba() instead of svd() in MySvd().
## 6. Add cs.keep method.

require(glmnet)
require(irlba)

estimateSuperPCA <- function(x, y, ncomps, 
             nCovs = NULL, 
             nCovs.min = ifelse(is.null(nCovs), max(ncomps)+2, min(nCovs)), 
             nCovs.max = ifelse(is.null(nCovs), nrow(x), max(nCovs)),
             nthresh = ifelse(is.null(nCovs), 25, length(nCovs)),
             nCovs.select = NULL, 
             nCovs.min.select = ifelse(is.null(nCovs.select), max(ncomps)+2, min(nCovs.select)), 
             nCovs.max.select = ifelse(is.null(nCovs.select), nrow(x), max(nCovs.select)),
             nthresh.select = ifelse(is.null(nCovs.select), 25, length(nCovs.select)),
             PCAmethod = c("SPC", "SPC.lasso", "cs", "cs.select"),
             kfold = 10,
             progress = FALSE){
    ## Purpose: the main function, automatically chooses best number of covariates,
    ##          best number of component for SPC/SPC.lasso/cs, 
    ##          and best number of covariates in selecting step for cs.select
    ##          by minimizing the CV error, returns the estimated 
    ##          supervised PCA object with the optimal values for tuning parameters.
    ## Inputs: data matrices: x (dim(n,p), and y (an n-vector) required
    ##         ncomps: required, denotes number of components, 
    ##                can be a number or a vector of integers
    ##         nCovs: optional, a vector of possible numbers of covariates
    ##         nCovs.min: optional, the smallest number of covariates
    ##                    default as max(ncomps)+2
    ##         nCovs.max: optional, the largest number of covariates
    ##                    default as number of rows of x
    ##         nthresh: optional, how many nCov to be tested
    ##                  default as 25
    ##         nCovs.select: optional, a vector of possible numbers of covariates 
    ##                       in final selecting process
    ##         nCovs.min.select: optional, the smallest number of covariates
    ##                    in final selecting process, default as max(ncomps)+2
    ##         nCovs.max.select: optional, the largest number of covariates
    ##                    in final selecting process, default as number of rows of x
    ##         nthresh.select: optional, how many nCov.select to be tested 
    ##                  in final selecting process, default as 25    
    ##         kfold: required, indicates the number of k in kfold cross-validation,
    ##                default as 10
    ##         progress: logical, print how many folds have been calculated and the 
    ##                   the current number of covariates
    ## Outputs: two lists the first is 'unfixTuning' of class 'supervisedPCACV' 
    ##                            which contains tuning parameter information
    ##                the second is 'fixTuning' of class 'supervisedPCA' 
    ##                            which is the estimated model
    if(PCAmethod == "SPC" | PCAmethod == "SPC.lasso"){
        unfixTuning <- findThreshold.SPC(x=x, y=y, ncomps=ncomps, nCovs=nCovs, 
                                         nCovs.min=nCovs.min, nCovs.max=nCovs.max,
                                         nthresh=nthresh,
                                         PCAmethod=PCAmethod,
                                         kfold=kfold, progress=progress)
        fixTuning <- approximatPCR(x=x, y=y, ncomp=unfixTuning$ncomp.best, 
                                   nCov=unfixTuning$nCov.best, PCAmethod=PCAmethod)
    } else if (PCAmethod == "cs") {
        unfixTuning <- findThreshold.cs(x=x, y=y, ncomps=ncomps, nCovs=nCovs,
                                        nCovs.min=nCovs.min, nCovs.max=nCovs.max,
                                        nthresh=nthresh,
                                        kfold=kfold, progress=progress)
        fixTuning <- approximatPCR(x=x, y=y, ncomp=unfixTuning$ncomp.best, 
                                   nCov=unfixTuning$nCov.best, PCAmethod=PCAmethod)
    } else if (PCAmethod == "cs.select") {
        unfixTuning <- findThreshold.select(x=x, y=y, ncomps=ncomps, nCovs=nCovs,
                                            nCovs.min=nCovs.min, nCovs.max=nCovs.max,
                                            nthresh=nthresh,
                                            nCovs.select=nCovs.select,
                                            nCovs.min.select=nCovs.min.select,
                                            nCovs.max.select=nCovs.max.select,
                                            nthresh.select=nthresh.select, 
                                            kfold=kfold, progress=progress)
        fixTuning <- approximatPCR(x=x, y=y, ncomp=unfixTuning$ncomp.best, 
                                   nCov=unfixTuning$nCov.best, 
                                   nCov.select=unfixTuning$nCov.select.best,
                                   PCAmethod=PCAmethod)
    }
    out <- append(unfixTuning, fixTuning)
    class(out) <- c('supervisedPCACV', 'supervisedPCA')
    return(out)
}


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



approximatPCR <- function (x, y, ncomp, nCov, nCov.select=NULL,
                      PCAmethod=c("SPC", "SPC.lasso", "cs", "cs.select"),
                      tStats=NULL){
    ## Purpose: performs SPC/SPC.lasso/cs/cs.select for fixed ncomp, nCov and nCov.select
    ## Inputs: data matrices: x (dim(n,p), and y (an n-vector) required
    ##         ncomp: required, number of components
    ##         nCov: required, number of covariates to compute SPC/SPC.lasso/cs
    ##         nCov.select: required when PCAmethod="cs.select", can be < or = or > than nCov
    ##         PCAmethod: required
    ##         tStates: optional
    ## Outputs: a list with many components, includes much of the input, of class 'supervisedPCA'
    ##          keep: logical, which columns of x are used to compute SPC/cs
    ##          tStats: see above
    ##          gamhat: regression coefficients on UD
    ##          betahat: regression coefficients on the kept x columns(order is changed)
    ##          ymean: the mean of y
    ##          xmeans: the column means of x, for out-of-sample prediction
    ##          u: approximated left singular vectors
    ##          d: approximated sigular values
    ##          v: approximated right sigular vectors
    if(is.null(tStats)){
        tStats <- marginalRegressionT(x=x, y=y) 
    }
    threshold <- quantile(abs(tStats), 1-nCov/ncol(x))
    keep <- abs(tStats) > threshold
    switch(PCAmethod,
           SPC = {
               smallx <- x[ , keep, drop=FALSE]
               smallx.svd <- MySvd(x=smallx, ncomponent=ncomp)
               u <- smallx.svd$u
               d <- smallx.svd$d
               v <- smallx.svd$v
               xmeans <- smallx.svd$xmeans
               ymean <- mean(y)
               ## the model is y-ymean = ud %*% gam
               ## so estimated gam is inv(d) %*% t(u) %*% (y-ymean)
               gamhat <- crossprod(scale(u, center=FALSE, scale=d), y-ymean) #dim(ncomp x 1)
               ## The betas are v %*% gamhat
               betahat <- v %*% gamhat #dim(keep x 1)
               out <- list(keep = keep, tStats = tStats, gamhat = gamhat, betahat = betahat, 
                           ymean = ymean, xmeans = xmeans, u = u, v = v, d = d, x = x, y = y, 
                           nCov = nCov, ncomp = ncomp)
           },
           SPC.lasso = {
               supervisedPCA <- approximatPCR(x=x, y=y, ncomp=ncomp, nCov=nCov, 
                                              PCAmethod="SPC", tStats=tStats)
               yhat <- predict(obj=supervisedPCA, newx = NULL, PCAmethod="SPC")
               cvfit <- cv.glmnet(x=x[ , keep, drop=FALSE], y=yhat, alpha=1, 
                                  intercept=TRUE, lambda.min.ratio=0.0001)  
               betahat <- coef(cvfit, s = "lambda.min")
               # betahat includes intercept
               out <- list(keep = keep, tStats = tStats, gamhat = supervisedPCA$gamhat,
                           betahat = betahat, ymean = supervisedPCA$ymean, 
                           xmeans = supervisedPCA$xmeans,
                           u = supervisedPCA$u, v = supervisedPCA$v, d = supervisedPCA$d,
                           x = x, y = y, nCov = nCov, ncomp = ncomp, 
                           b.keep = (abs(betahat)>0), cvfit = cvfit) # b.keep includes intercept
           },
           cs = {
               x1 <- x[ , keep, drop=FALSE]
               x.new <- cbind(x1, x[ , !keep, drop=FALSE])
               L.S <- crossprod(x.new, x1)  ## t(x.new) %*% x1
               svd <- MySvd(x=L.S, ncomponent=ncomp)
               v <- svd$u[ , 1:ncomp, drop=FALSE]
               d <- svd$d[1:ncomp]^(1/2) #* (p/l)^(1/4)
               u <- scale(x.new %*% v, center=FALSE, scale=d)
               xmeans <- colMeans(x.new)
               ymean <- mean(y)
               ## the model is y-ymean = ud %*% gam
               ## so estimated gam is inv(d) %*% t(u) %*% (y-ymean)
               gamhat <- crossprod(scale(u, center=FALSE, scale=d), y-ymean) #dim(ncomp x 1)
               ## The betas are v %*% gamhat
               betahat <- v %*% gamhat #dim(p x 1)
               out <- list(keep = keep, tStats = tStats, gamhat = gamhat, 
                           betahat = betahat, ymean = ymean, xmeans = xmeans, 
                           u = u, v = v, d = d, x = x, y = y, nCov = nCov, 
                           ncomp = ncomp, b.keep = (abs(betahat)>0))
           },
           cs.select = {
               cs <- approximatPCR(x=x, y=y, ncomp=ncomp, nCov=nCov, 
                                   PCAmethod="cs", tStats=tStats)
               b.threshold <- quantile(abs(cs$betahat), 1-nCov.select/ncol(x))
               b.keep <- abs(cs$betahat) > b.threshold
               out <- list(keep = keep, tStats = tStats, gamhat = cs$gamhat,
                           betahat = cs$betahat[b.keep], 
                           ymean = cs$ymean, xmeans = cs$xmeans, 
                           u = cs$u, v = cs$v, d = cs$d, 
                           x = x, y = y, nCov = nCov, ncomp = ncomp, 
                           b.keep = b.keep, nCov.select = nCov.select)
           }
    )
    class(out) <- 'supervisedPCA'
    return(out)
}




findThreshold.SPC <- function (x, y, ncomps, nCovs = NULL, 
                               nCovs.min = ifelse(is.null(nCovs), max(ncomps)+2, min(nCovs)), 
                               nCovs.max = ifelse(is.null(nCovs), nrow(x), max(nCovs)),
                               nthresh = ifelse(is.null(nCovs), 25, length(nCovs)),
                               PCAmethod = c("SPC", "SPC.lasso"),
                               kfold = 10, progress = FALSE) {
    ## Purpose: find the optimal number of covariates and number of components
    ##          using kfold cross-validation for SPC and SPC.lasso
    ## Inputs: data matrices: x (dim(n,p), and y (an n-vector) required
    ##         ncomps: required, denotes number of components, 
    ##                can be a number or a vector
    ##         nCovs: optional, a vector of possible numbers of covariates
    ##         nCovs.min: optional, the smallest number of covariates
    ##                    default as max(ncomps)+2
    ##         nCovs.max: optional, the largest number of covariates
    ##                    default as number of rows of x
    ##         nthresh: optional, how many nCov to be tested
    ##                  default as 25
    ##         PCAmethod: required, either SPC or SPC.lasso
    ##         kfold: required, indicates the number of k in kfold cross-validation,
    ##                default as 10
    ##         progress: logical, print how many nCov/folds have been calculated
    ## Outputs: an object of class 'supervisedPCACV', a list with the following components
    ##          nCov.best：the best number of covariates which gives smallest mse in CV
    ##          ncomp.best: the best number of components which gives smallest mse in CV
    ##          ncomps: all the tested ncomps
    ##          nCovs: all the tested numbers of covariates
    ##          CVmse: mse in CV with respect to each value of nCovs, ncomps, kfold
    ##          mse: average mse in CV over k folds
    if(is.null(nCovs)){
        nCovs <- round(seq(from=nCovs.min, to=nCovs.max, length.out=nthresh))
    } 
    dataCV <- MykfoldCV(x=x, y=y, k=kfold)
    # CVmse store all the mse data over various fold, ncomp, nCov
    CVmse <- array(NA, dim=c(nthresh, length(ncomps), kfold))
    for(k in 1:kfold){
        # prepare training data and testing data
        testingX <- dataCV$dataset[[k]]$testing.x
        testingY <- dataCV$dataset[[k]]$testing.y
        trainingX <- dataCV$dataset[[k]]$training.x
        trainingY <- dataCV$dataset[[k]]$training.y
        # compute tStats for training data
        tStats <- marginalRegressionT(x=trainingX, y=trainingY) 
        for(i in 1:nthresh) {
            threshold <- quantile(abs(tStats), 1-nCovs[i]/ncol(trainingX))
            keep <- abs(tStats) > threshold
            smallx <- trainingX[ , keep, drop=FALSE]
            smallx.svd <- MySvd(x=smallx, ncomponent=max(ncomps))
            xmeans <- smallx.svd$xmeans
            ymean <- mean(trainingY)
            for(j in 1:length(ncomps)){
                u <- smallx.svd$u[, 1:ncomps[j]]
                d <- smallx.svd$d[1:ncomps[j]]
                v <- smallx.svd$v[, 1:ncomps[j]]
                ## the model is y-ymean = ud %*% gam
                ## so estimated gam is inv(d) %*% t(u) %*% (y-ymean)
                gamhat <- crossprod(scale(u, center=FALSE, scale=d), trainingY-ymean) #dim(ncomp x 1)
                ## The betas are v %*% gamhat
                betahat <- v %*% gamhat #dim(keep x 1)
                if(PCAmethod == "SPC"){
                    smallx.test <- testingX[ , keep, drop=FALSE]
                    smallx.test <- scale(smallx.test, center=xmeans, scale=FALSE)
                    testingY.hat <- ymean + smallx.test %*% betahat
                    CVmse[i, j, k] <- mean((testingY - testingY.hat)^2)
                } else if (PCAmethod == "SPC.lasso") {
                    smallx.scale <- scale(smallx, center=xmeans, scale=FALSE)
                    trainingY.hat <- ymean + smallx.scale %*% betahat
                    cvfit <- cv.glmnet(x=smallx, y=trainingY.hat, alpha=1, 
                                       intercept=TRUE, lambda.min.ratio=0.0001)  
                    smallx.test <- testingX[ , keep, drop=FALSE]
                    testingY.hat <- predict(cvfit, newx=smallx.test, s="lambda.min")
                    CVmse[i, j, k] <- mean((testingY - testingY.hat)^2)
                }
            }
            if(progress==TRUE) cat('nCov', nCovs[i],'\n')
        }
        if(progress==TRUE) cat('Fold', k, ' / ', kfold,'\n')  
    }
    mse <- findMeanForArray3(CVmse)
    out <- list(nCov.best = nCovs[which(mse == min(mse), arr.ind = TRUE)[1]], 
                ncomp.best = ncomps[which(mse == min(mse), arr.ind = TRUE)[2]],
                ncomps = ncomps, nCovs = nCovs, CVmse = CVmse, mse = mse)
    class(out) <- 'supervisedPCACV'
    return(out)
}




findThreshold.cs <- function (x, y, ncomps, nCovs = NULL, 
                    nCovs.min = ifelse(is.null(nCovs), max(ncomps)+2, min(nCovs)), 
                    nCovs.max = ifelse(is.null(nCovs), nrow(x), max(nCovs)),
                    nthresh = ifelse(is.null(nCovs), 25, length(nCovs)),
                    kfold = 10, progress = FALSE) {
    ## Purpose: find the optimal number of covariates and number of components
    ##          using kfold cross-validation for cs
    ## Inputs: data matrices: x (dim(n,p), and y (an n-vector) required
    ##         ncomps: required, denotes number of components, 
    ##                can be a number or a vector
    ##         nCovs: optional, a vector of possible numbers of covariates
    ##         nCovs.min: optional, the smallest number of covariates
    ##                    default as max(ncomps)+2
    ##         nCovs.max: optional, the largest number of covariates
    ##                    default as number of rows of x
    ##         nthresh: optional, how many nCov to be tested
    ##                  default as 25
    ##         kfold: required, indicates the number of k in kfold cross-validation,
    ##                default as 10
    ##         progress: logical, print how many nCov/folds have been calculated
    ## Outputs: an object of class 'supervisedPCACV', a list with the following components
    ##          nCov.best：the best number of covariates which gives smallest mse in CV
    ##          ncomp.best: the best number of components which gives smallest mse in CV
    ##          ncomps: all the tested ncomps
    ##          nCovs: all the tested numbers of covariates
    ##          CVmse: mse in CV with respect to each value of  nCovs, ncomps, kfold
    ##          mse: average mse in CV over k folds
    if(is.null(nCovs)){
        nCovs <- round(seq(from=nCovs.min, to=nCovs.max, length.out=nthresh))
    } 
    dataCV <- MykfoldCV(x=x, y=y, k=kfold)
    # CVmse store all the mse data over various fold, ncomp, nCov
    CVmse <- array(NA, dim=c(nthresh, length(ncomps), kfold))
    for(k in 1:kfold){
        # prepare training data and testing data
        testingX <- dataCV$dataset[[k]]$testing.x
        testingY <- dataCV$dataset[[k]]$testing.y
        trainingX <- dataCV$dataset[[k]]$training.x
        trainingY <- dataCV$dataset[[k]]$training.y
        # compute tStats for training data
        tStats <- marginalRegressionT(x=trainingX, y=trainingY) 
        for(i in 1:nthresh) {
            threshold <- quantile(abs(tStats), 1-nCovs[i]/ncol(trainingX))
            keep <- abs(tStats) > threshold
            x1 <- trainingX[ , keep, drop=FALSE]
            x.new <- cbind(x1, trainingX[ , !keep, drop=FALSE])
            L.S <- crossprod(x.new, x1)  ## t(x.new) %*% x1
            svd <- MySvd(x=L.S, ncomponent=max(ncomps))
            xmeans <- colMeans(x.new)
            ymean <- mean(trainingY)
            for(j in 1:length(ncomps)){
                v <- svd$u[ , 1:ncomps[j], drop=FALSE]
                d <- svd$d[1:ncomps[j]]^(1/2) #* (p/l)^(1/4)
                u <- scale(x.new %*% v, center=FALSE, scale=d)
                ## the model is y-ymean = ud %*% gam
                ## so estimated gam is inv(d) %*% t(u) %*% (y-ymean)
                gamhat <- crossprod(scale(u, center=FALSE, scale=d), trainingY-ymean) #dim(ncomp x 1)
                ## The betas are v %*% gamhat
                betahat <- v %*% gamhat #dim(p x 1)
                testingX.new <- cbind(testingX[ , keep, drop=FALSE], testingX[ , !keep, drop=FALSE])
                testingX.new <- scale(testingX.new, center=xmeans, scale=FALSE)
                testingY.hat <- ymean + testingX.new %*% betahat
                CVmse[i, j, k] <- mean((testingY - testingY.hat)^2)
            }
            if(progress==TRUE) cat('nCov', nCovs[i],'\n')
        }
        if(progress==TRUE) cat('Fold', k, ' / ', kfold,'\n')  
    }
    mse <- findMeanForArray3(CVmse)
    out <- list(nCov.best = nCovs[which(mse == min(mse), arr.ind = TRUE)[1]], 
                ncomp.best = ncomps[which(mse == min(mse), arr.ind = TRUE)[2]],
                ncomps = ncomps, nCovs = nCovs, CVmse = CVmse, mse = mse)
    class(out) <- 'supervisedPCACV'
    return(out)
}



findThreshold.select <- function (x, y, ncomps, nCovs = NULL, 
          nCovs.min = ifelse(is.null(nCovs), max(ncomps)+2, min(nCovs)), 
          nCovs.max = ifelse(is.null(nCovs), nrow(x), max(nCovs)),
          nthresh = ifelse(is.null(nCovs), 25, length(nCovs)),
          nCovs.select = NULL, 
          nCovs.min.select = ifelse(is.null(nCovs.select), max(ncomps)+2, min(nCovs.select)), 
          nCovs.max.select = ifelse(is.null(nCovs.select), nrow(x), max(nCovs.select)),
          nthresh.select = ifelse(is.null(nCovs.select), 25, length(nCovs.select)),
          kfold = 10, progress = FALSE) {
    ## Purpose: find the optimal number of covariates, number of components, 
    ##          and number of covariates in selecting process for cs.select method
    ##          using kfold cross-validation
    ## Inputs: data matrices: x (dim(n,p), and y (an n-vector) required
    ##         ncomps: required, denotes number of components, 
    ##                can be a number or a vector
    ##         nCovs: optional, a vector of possible numbers of covariates
    ##         nCovs.min: optional, the smallest number of covariates
    ##                    default as max(ncomps)+2
    ##         nCovs.max: optional, the largest number of covariates
    ##                    default as number of rows of x
    ##         nthresh: optional, how many nCov to be tested
    ##                  default as 25
    ##         nCovs.select: optional, a vector of possible numbers of covariates 
    ##                       in final selecting process
    ##         nCovs.min.select: optional, the smallest number of covariates
    ##                    in final selecting process, default as max(ncomps)+2
    ##         nCovs.max.select: optional, the largest number of covariates
    ##                    in final selecting process, default as number of rows of x
    ##         nthresh.select: optional, how many nCov.select to be tested 
    ##                  in final selecting process, default as 25    
    ##         kfold: required, indicates the number of k in kfold cross-validation,
    ##                default as 10
    ##         progress: logical, print how many nCov/folds have been calculated
    ## Outputs: an object of class 'supervisedPCACV', a list with the following components
    ##          nCov.select.best：the best number of covariates in selecting process
    ##                            which gives smallest mse in CV
    ##          ncomp.best: the best number of components which gives smallest mse in CV
    ##          nCov.best：the best number of covariates which gives smallest mse in CV  
    ##          ncomps: all the tested ncomps
    ##          nCovs: all the tested numbers of covariates
    ##          nCovs.select: all the tested numbers of covariates in selecting process
    ##          CVmse: mse in CV with respect to each value of  nCovs, ncomps, kfold
    ##          mse: average mse in CV over k folds
    if(is.null(nCovs)){
        nCovs <- round(seq(from=nCovs.min, to=nCovs.max, length.out=nthresh))
    } 
    if(is.null(nCovs.select)){
        nCovs.select <- round(seq(from=nCovs.min.select, to=nCovs.max.select, 
                                  length.out=nthresh.select))
    } 
    dataCV <- MykfoldCV(x=x, y=y, k=kfold)
    # CVmse store all the mse data over various fold, ncomp, nCov, nCov.select
    CVmse <- array(NA, dim=c(nthresh.select, length(ncomps), nthresh, kfold))
    for(k in 1:kfold){
        # prepare training data and testing data
        testingX <- dataCV$dataset[[k]]$testing.x
        testingY <- dataCV$dataset[[k]]$testing.y
        trainingX <- dataCV$dataset[[k]]$training.x
        trainingY <- dataCV$dataset[[k]]$training.y
        # compute tStats for training data
        tStats <- marginalRegressionT(x=trainingX, y=trainingY) 
        for(i in 1:nthresh) {
            threshold <- quantile(abs(tStats), 1-nCovs[i]/ncol(trainingX))
            keep <- abs(tStats) > threshold
            x1 <- trainingX[ , keep, drop=FALSE]
            x.new <- cbind(x1, trainingX[ , !keep, drop=FALSE])
            L.S <- crossprod(x.new, x1)  ## t(x.new) %*% x1
            svd <- MySvd(x=L.S, ncomponent=max(ncomps))
            xmeans <- colMeans(x.new)
            ymean <- mean(trainingY)
            for(j in 1:length(ncomps)){
                v <- svd$u[ , 1:ncomps[j], drop=FALSE]
                d <- svd$d[1:ncomps[j]]^(1/2) #* (p/l)^(1/4)
                u <- scale(x.new %*% v, center=FALSE, scale=d)
                ## the model is y-ymean = ud %*% gam
                ## so estimated gam is inv(d) %*% t(u) %*% (y-ymean)
                gamhat <- crossprod(scale(u, center=FALSE, scale=d), trainingY-ymean) #dim(ncomp x 1)
                ## The betas are v %*% gamhat
                betahat <- v %*% gamhat #dim(p x 1)
                # select process
                for(l in 1:nthresh.select){
                    b.threshold <- quantile(abs(betahat), 1-nCovs.select[l]/ncol(trainingX))
                    b.keep <- abs(betahat) > b.threshold
                    testingX.new <- cbind(testingX[, keep, drop=FALSE], testingX[, !keep, drop=FALSE])
                    testingX.new <- scale(testingX.new, center=xmeans, scale=FALSE)
                    testingY.hat <- ymean + testingX.new[, b.keep, drop=FALSE] %*% betahat[b.keep]
                    CVmse[l, j, i, k] <- mean((testingY - testingY.hat)^2)
                }
            }
            if(progress==TRUE) cat('nCov', nCovs[i],'\n')
        }
        if(progress==TRUE) cat('Fold', k, ' / ', kfold,'\n')  
    }
    mse <- findMeanForArray4(CVmse)
    out <- list(nCov.select.best = nCovs.select[which(mse == min(mse), arr.ind = TRUE)[1]], 
                ncomp.best = ncomps[which(mse == min(mse), arr.ind = TRUE)[2]],
                nCov.best = nCovs[which(mse == min(mse), arr.ind = TRUE)[3]],
                ncomps = ncomps, nCovs = nCovs, nCovs.select = nCovs.select, 
                CVmse = CVmse, mse = mse)
    class(out) <- 'supervisedPCACV'
    return(out)
}





## The following are methods for the output of the main functions. They work just like the generic functions in base R
## i.e. you can say predict(superPCAout) to get the predicted values

predict.supervisedPCA <- function(obj, newx = NULL, 
                                  PCAmethod=c("SPC", "SPC.lasso", "cs", "cs.select")){
    ## Inputs: an object of class 'supervisedPCA' as output by approximatPCR
    ##         an optional matrix of points for predictions, if missing, then uses original data
    ## Outputs: predictions for newx design matrix
    if(is.null(newx)) newx <- obj$x  # can be used to compute fitted values
    switch (PCAmethod,
            SPC = {
                smallx <- newx[ , obj$keep, drop=FALSE]
                smallx <- scale(smallx, center=obj$xmeans, scale=FALSE)
                ## Note: at this point, their code does some screwy scalings (superpc.predict ~30).
                ## Check results later.
                yhat <- obj$ymean + smallx %*% obj$betahat
                return(yhat)
            },
            SPC.lasso = {
                smallx <- newx[ , obj$keep, drop=FALSE]
                yhat <- predict(obj$cvfit, newx=smallx, s="lambda.min")
                return(yhat)
            },
            cs = {
                x.new <- cbind(newx[ , obj$keep, drop=FALSE], newx[ , !obj$keep, drop=FALSE])
                x.new <- scale(x.new, center=obj$xmeans, scale=FALSE)
                yhat <- obj$ymean + x.new %*% obj$betahat
                return(yhat)
            },
            cs.select = {
                x.new <- cbind(newx[ , obj$keep, drop=FALSE], newx[ , !obj$keep, drop=FALSE])
                x.new <- scale(x.new, center=obj$xmeans, scale=FALSE)
                yhat <- obj$ymean + x.new[, obj$b.keep, drop=FALSE] %*% obj$betahat
                return(yhat)
            }
    )
}


## Some helper functions which are not meant to be called directly

MySvd <- function(x, ncomponent=min(nrow(x), ncol(x))) {
    ## Purpose: performs svd quickly, if ncomponent <= min(nrow(x), ncol(x)) * .5, use irlba,
    ##          otherwise use svd directly
    ## Inputs: an n x p matrix x
    ##             an (optional) number of desired components, defaults to "all"
    ## Outputs: a list with components u, d, v, and xmeans
    ##             u is a matrix of left singular vectors of dimension n x ncomponent
    ##             v is a matrix of right singular vectors (not t(v)) of dimension p x ncomponent
    ##             d is a vector of singular values of length ncomponent
    ##             xmeans is a p-vector of the column means of x
    xmeans <- colMeans(x)
    x <- scale(x, center = xmeans, scale=FALSE)
# I think we  should require the input of ncomp to be larger than 0
    
#     if(ncomponent==0){
#         s <- svd(x, nu=0, nv=0) 
#         return(list(d = s$d, xmeans=xmeans))
#     }
    
    if (ncomponent <= min(nrow(x), ncol(x)) * .5) {
        s <- irlba(A=x, nv=ncomponent, nu=ncomponent)
        u <- s$u[,1:ncomponent,drop=FALSE] 
        d <- s$d[1:ncomponent]
        v <- s$v[,1:ncomponent,drop=FALSE]
    } else {
        s <- svd(x=x, nv=ncomponent, nu=ncomponent)
        u <- s$u[,1:ncomponent,drop=FALSE] 
        d <- s$d[1:ncomponent]
        v <- s$v[,1:ncomponent,drop=FALSE]
    }
    return(list(u=u, d=d, v=v, xmeans=xmeans)) # v not t(v)
}


MykfoldCV <- function(x, y, k=10){
    ## Purpose: prepare training set and testing set for kfold cross-validation
    ## Inputs: an n*p matrix x as data matrix;
    ##         an n*1 vector y as response;
    ##         k is the number of folds, by default is 10.
    ## Outputs: a list of 2 components, first is the list of dataset including k lists, 
    ##          each includes sets of training.x, training.y, testing.x, testing.y;
    ##          second is the folds list which indicates the index for each fold.
    
    # randomly select index for partition of rows of x
    n <- nrow(x)
    folds <- vector("list", k)
    breaks <- round(seq(from = 1, to = (n + 1), length = (k + 1)))
    cv.order <- sample(1:n)
    dataset <- vector("list", k)
    for(i in 1:k){
        # prepare index for the ith fold
        folds[[i]] <- cv.order[(breaks[i]):(breaks[i + 1] - 1)]
        # generate training set and testing set for ith fold
        testing.x <- x[folds[[i]], ]
        testing.y <- y[folds[[i]]]
        training.x <- x[-folds[[i]], ]
        training.y <- y[-folds[[i]]]
        dataset[[i]] <- list(testing.x = testing.x, testing.y = testing.y, 
                         training.x = training.x, training.y=training.y)
    }
    out <- list(dataset=dataset, folds=folds)
    return(out)
}


findMeanForArray3 <- function(x){
    ## Purpose: find the mean over the 3rd dimention of a 3-dimentional array
    ## Inputs: x is a 3-dimentional array
    ## Outputs: a matrix of the mean over the 3rd dimention
    n <- dim(x)[1]
    p <- dim(x)[2]
    q <- dim(x)[3]
    out <- matrix(NA, n, p)
    for (i in 1:n){
        for(j in 1:p){
            tot <- 0
            for(k in 1:q){
                tot <- tot + x[i, j, k]
            }
            out[i, j] <- tot/q
        }
    }
    return(out)
}



findMeanForArray4 <- function(x){
    ## Purpose: find the mean over the 4th dimention of a 4-dimentional array
    ## Inputs: x is a 4-dimentional array
    ## Outputs: a 3-dimentional array
    n <- dim(x)[1]
    p <- dim(x)[2]
    q <- dim(x)[3]
    m <- dim(x)[4]
    out <- array(NA, c(n, p, q))
    for (i in 1:n){
        for(j in 1:p){
            for(k in 1:q){
                tot <- 0
                for(l in 1:m){
                    tot <- tot + x[i, j, k, l]
                }
                out[i, j, k] <- tot/m 
            }
        }
    }
    return(out)
}
