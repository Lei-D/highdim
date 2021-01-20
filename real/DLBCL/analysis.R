source("~/Desktop/Lei/AIMER2/github/sparsePC_reg_real.R")

# read in data
staudt.x = matrix(scan("~/Desktop/Lei/AIMER/DLBCL/rawData/staudt.x"), ncol=240, byrow=TRUE) 
staudt.tim = scan("~/Desktop/Lei/AIMER/DLBCL/rawData/staudt.tim") 

# transform covariates and response
Xtot = t(staudt.x)
Xtot = scale(Xtot, center = T, scale = T)
Ytot = log(staudt.tim + 1)

# set lambda values
lambda.min = 0
lambda.max = 1
nlambda = 10
lambda = lambda.max * log10(seq(1, 10, length.out = nlambda)) + 
  lambda.min * (1 - log10(seq(1, 10, length.out = nlambda)))
lambda = sort(lambda, decreasing = T)

# suffPCR
set.seed(8372)
suffPCR = sparsePC.reg.real.CV(X = Xtot, Y = Ytot, d = 3,
                               lambda = lambda, Kfold = 5,
                               maxnvar = ncol(Xtot),
                               screening = TRUE)

# save output
save(suffPCR, file = "~/Desktop/Lei/AIMER2/ISMB/real/DLBCL/suffPCR.RData")


