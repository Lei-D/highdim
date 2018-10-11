library(fps)
library(irlba)

PCR.fps.logistic <- function(X, Y, Xtest, Ytest, ncomp, ncov, covPercent){
    n = nrow(X)
    
    # center X, Y
    xmeans = colMeans(X)
    X = scale(X, center=xmeans, scale=FALSE)
    
    # fps
    approximate = fps(crossprod(X)/n, ndim=ncomp, nsol=2, maxnvar=ncov, covPercent=covPercent)
    decomp = irlba(approximate$projection[[2]], nv=0, nu=ncomp)
    vhat = decomp$u
    pchat = X %*% vhat
    
    # fit logistic
    model = glm(Y~pchat, family=binomial(link='logit'))
    gamhat = as.vector(model$coefficients[-1])
    betahat = vhat %*% t(t(gamhat))
    fitted = X %*% betahat
    
    # find the best cutoff
    cutoffs = sort(fitted)
    accuracy = rep(NA, n)
    for (i in 1:n){
        prediction = ifelse(fitted >= cutoffs[i], 1, 0) 
        accuracy[i] = sum(Y == prediction) / n * 100
    }
    cutoff.best = cutoffs[which.max(accuracy)]
    
    # predict
    Xtest = scale(Xtest, center=xmeans, scale=FALSE)
    yhat.odd = Xtest %*% betahat
    yhat = ifelse(yhat.odd >= cutoff.best, 1, 0)
    accuracy = sum(Ytest == yhat) / length(Ytest) * 100
    out = list(accuracy = accuracy, ncov = sum(betahat!=0))
    return(out)
}



PCR.fps.logistic.CV = function(X, Y, Xtest, Ytest, ncomp, ncov, kfold, covPercent) {
    dataCV = MykfoldCV(x=X, y=Y, k=kfold)
    # CVaccuracy store all the mse data over various fold, ncomp, NCov
    CVaccuracy = array(NA, dim=c(length(ncov), length(ncomp), kfold))
    for(k in 1:kfold){
        cat('fold', k,'\n')
        
        # prepare training data and testing data
        testingX = dataCV$dataset[[k]]$testing.x
        testingY = dataCV$dataset[[k]]$testing.y
        trainingX = dataCV$dataset[[k]]$training.x
        trainingY = dataCV$dataset[[k]]$training.y
        
        for(i in 1:length(ncov)) {
            cat('ncov', ncov[i],'\n')
            for (j in 1:length(ncomp)) {
                cat('ncomp', ncomp[j],'\n')
                m = PCR.fps.logistic(X=trainingX, Y=trainingY, Xtest=testingX, Ytest=testingY, 
                                     ncomp=ncomp[j], ncov=ncov[i], covPercent=covPercent)
                CVaccuracy[i,j,k] = m$accuracy
            }
        }
    }
    Accuracy = apply(CVaccuracy, c(1, 2), mean)
    ncov.best = ncov[which(Accuracy == max(Accuracy), arr.ind = TRUE)[1]]
    ncomp.best = ncomp[which(Accuracy == max(Accuracy), arr.ind = TRUE)[2]]
    out = list(ncov.best = ncov.best, ncomp.best = ncomp.best, Accuracy.best = max(Accuracy))
}





PCR.fps.logistic.CV2 = function(X, Y, Xtest, Ytest
                                , ncomp, ncov.max
                                , nsol, kfold
                                # , covPercent
                                ) {
    dataCV = MykfoldCV(x=X, y=Y, k=kfold)
    # CVaccuracy store all the prediction accuracy data over various fold, ncomp, ncov
    CVaccuracy = array(NA, dim=c(length(ncomp), nsol, kfold))
    ncov = array(NA, dim=c(length(ncomp), nsol, kfold))
    for(k in 1:kfold){
        cat('......fold', k, '/', kfold, '\n')
        
        # prepare training data and testing data
        testingX = dataCV$dataset[[k]]$testing.x
        testingY = dataCV$dataset[[k]]$testing.y
        trainingX = dataCV$dataset[[k]]$training.x
        trainingY = dataCV$dataset[[k]]$training.y
        
        for (i in 1:length(ncomp)) {
            n = nrow(trainingX)
            
            # center trainingX, trainingY
            xmeans = colMeans(trainingX)
            trainingX = scale(trainingX, center=xmeans, scale=FALSE)
            
            # fps
            t0 = Sys.time()
            approximate = fps(crossprod(trainingX)/n
                              , ndim=ncomp[i], nsol=nsol
                              # , maxnvar=ncov.max
                              , covPercent=covPercent
            )
            t1 = Sys.time()
            print(difftime(t1, t0, units="mins"))
            
            # compute accuracy for each solution
            for (j in 1:nsol) {
                decomp = irlba(approximate$projection[[j]], nv=0, nu=ncomp[i])
                vhat = decomp$u
                pchat = trainingX %*% vhat
                
                # fit logistic
                model = glm(trainingY ~ pchat, family=binomial(link='logit'))
                gamhat = as.vector(model$coefficients[-1])
                betahat = vhat %*% t(t(gamhat))
                fitted = trainingX %*% betahat
                
                # find the best cutoff
                # cutoffs = sort(fitted)
                # accuracy = rep(NA, n)
                # for (i in 1:n){
                #   prediction = ifelse(fitted >= cutoffs[i], 1, 0) 
                #   accuracy[i] = sum(trainingY == prediction) / n * 100
                # }
                # cutoff.best = cutoffs[which.max(accuracy)]
                cutoff.best = 0.5
                
                # predict
                testingX = scale(testingX, center=xmeans, scale=FALSE)
                yhat.odd = testingX %*% betahat
                yhat = ifelse(yhat.odd >= cutoff.best, 1, 0)
                CVaccuracy[i,j,k] = sum(testingY == yhat) / length(testingY) * 100
                ncov[i,j,k] = sum(betahat!=0)
            }
        }
    }
    Accuracy = apply(CVaccuracy, c(1, 2), mean)
    Ncov = apply(ncov, c(1, 2), mean)
    ncomp.best = ncomp[which(Accuracy == max(Accuracy), arr.ind = TRUE)[1]]
    ncov.best = ncov[which(Accuracy == max(Accuracy), arr.ind = TRUE)[1], ncov.best = ncov[which(Accuracy == max(Accuracy), arr.ind = TRUE)[2]]]
    out = list(ncov.best = ncov.best, ncomp.best = ncomp.best, Accuracy.best = max(Accuracy))
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

