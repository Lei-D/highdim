source("~/Desktop/Lei/AIMER2/github/sparsePC_reg_real.R")

# read in data
load("~/Desktop/Lei/AIMER/NSCLC/rawdata/Chemores-miRNA-LVS-normalized.RData")
load("~/Desktop/Lei/AIMER/NSCLC/rawdata/Chemores-clindat-public-Apr2014.RData")

# condiser relapse time for patients who relapse
Xtot = t(mirdat)
Xtot = scale(Xtot, center = T, scale = T)
Ytot = log(chemores.clindat$Relapse.time + 1)

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
save(suffPCR, file = "~/Desktop/Lei/AIMER2/ISMB/real/NSCLC/suffPCR.RData")


