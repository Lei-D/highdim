# source functions
source("/Users/lei/Desktop/AIMER2/empirical/PCR_fps_logistic.R")

# indicate which dataset you wanna run
dataset =  c(2)
# c(1:4, 6, 8, 10:12)
niter = 10 # the number of training/test spliting
for (i in dataset) {
    cat("dataset",i,"running...", "\n")
    # read in data
    Xtot = as.matrix(read.csv(paste("/Users/lei/Desktop/AIMER2/empirical/cleandata/", i, "X.csv", sep=""), header=F))
    Ytot = as.matrix(read.csv(paste("/Users/lei/Desktop/AIMER2/empirical/cleandata/", i, "Y.csv", sep=""), header=F))
    Ytot = factor(Ytot, levels=sort(unique(Ytot)), labels=c(0, 1))
    Accuracy = matrix(NA, nrow=length(dataset), ncol=niter)
    Ncomp = matrix(NA, nrow=length(dataset), ncol=niter)
    Ncov = matrix(NA, nrow=length(dataset), ncol=niter)
    for (j in 1:niter) {
        cat("...iteration", j, "\n")
        # choose half randomly to be training set, the other half is test set
        set.seed(j+276) 
        sample_index = sort(sample.int(nrow(Xtot), as.integer(nrow(Xtot)*.5), replace = FALSE))
        X = Xtot[sample_index,]
        Xtest = Xtot[-sample_index,]
        Y = Ytot[sample_index]
        Ytest = Ytot[-sample_index]
        fit = PCR.fps.logistic.CV(X=X, Y=Y, Xtest=Xtest, Ytest=Ytest
                                  , ncomp=c(3, 5, 10, 20, 50)
                                  , ncov=c(seq(20, 100, 5), seq(200, 500, 100))
                                  , kfold=5
                                  , covPercent=1
                                  )
        print(fit$Accuracy.best)
        print(fit$ncomp.best)
        print(fit$ncov.best)
        
        # store output
        Accuracy[i,j] = fit$Accuracy.best
        Ncomp[i,j] = fit$ncomp.best
        Ncov[i,j] = fit$ncov.best
        
    }
    save(Accuracy, Ncomp, Ncov, file = paste("/Users/lei/Desktop/AIMER2/empirical/output/50comp", i, ".RData", sep=""))
}
