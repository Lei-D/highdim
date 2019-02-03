# function to simulate U from Gaussian
sim_U = function(n, d){
  U = matrix(rnorm(n*d), nrow = n, ncol = d)
  return(U)
}

# function to simulate V
sim_V = function(p, s, d){
  # p is the number of columns in X
  # s is the number of non-zero rows in V
  # d is the population dimension of space spaned by Sigma
  
  # generate non-zero rows of V from Gaussian
  V = matrix(rnorm(s*d), nrow = s, ncol = d)
  # orthonormalize non-zero entries of V
  V = far::orthonormalization(V, basis = FALSE, norm = TRUE)
  # add zero rows into V
  V = rbind(V, matrix(0, nrow = p-s, ncol = d))
  return(V)
}

# function to simulate lambda
sim_lambda = function(d){
  # generate lambda from uniform
  lambda = sort(runif(n = d, min = 1, max = d), decreasing = TRUE)
  return(lambda)
}

# function to simulate theta
sim_theta = function(d){
  # generate theta from uniform
  theta = sort(runif(n = d, min = 1, max = d), decreasing = TRUE)
  return(theta)
}

# generate training, validation, test set for continuous response
sim_reg = function(n, p, s, d, SNRx, SNRy){
  # n is the number of rows in X
  # p is the number of columns in X
  # s is the number of non-zero rows in V
  # d is the population dimension of space spaned by Sigma
  
  # generate fixed variables
  U = sim_U(n, d)
  V = sim_V(p, s, d)
  lambda = sim_lambda(d)
  theta = sim_theta(d)
  sigma0 = sqrt(sum(lambda^2))/sqrt(p)/SNRx
  beta = V %*% diag(lambda/(lambda^2 + sigma0^2)) %*% theta
  # sigma1 = as.numeric(sqrt(crossprod(diag(lambda)%*%t(V)%*%beta)/SNRy - sigma0^2 * crossprod(beta)))
  sigma1 = as.numeric(sqrt((crossprod(diag(lambda)%*%t(V)%*%beta) + sigma0^2 * crossprod(beta))/SNRy))
  
  # add noise and generate traing set
  Xtrain = U %*% diag(lambda) %*% t(V) + sigma0 * matrix(rnorm(n * p), n, p)
  Ytrain = U %*% theta + sigma1 * rnorm(n)
  
  # add noise and generate validation set
  Xvalidate = U %*% diag(lambda) %*% t(V) + sigma0 * matrix(rnorm(n * p), n, p)
  Yvalidate = U %*% theta + sigma1 * rnorm(n)
  
  # add noise and generate test set
  Xtest = U %*% diag(lambda) %*% t(V) + sigma0 * matrix(rnorm(n * p), n, p)
  Ytest = U %*% theta + sigma1 * rnorm(n)
  
  output = list(Xtrain = Xtrain, Ytrain = Ytrain
                , Xvalidate = Xvalidate, Yvalidate = Yvalidate
                , Xtest = Xtest, Ytest = Ytest
                , U = U, V = V
                , lambda = lambda, theta = theta
                , beta = beta
                , sigma0 = sigma0, sigma1 = sigma1)
  return(output)
}

# generate training, validation, test set for binary response
sim_class = function(n, p, s, d, SNRx, SNRy){
  # n is the number of rows in X
  # p is the number of columns in X
  # s is the number of non-zero rows in V
  # d is the population dimension of space spaned by Sigma
  
  reg = sim_reg(n, p, s, d, SNRx, SNRy)
  Ytrain = ifelse(reg$Ytrain >= 0, 1, 0)
  Yvalidate = ifelse(reg$Yvalidate >= 0, 1, 0)
  Ytest = ifelse(reg$Ytest >= 0, 1, 0)
  
  output = list(Xtrain = reg$Xtrain, Ytrain = Ytrain
                , Xvalidate = reg$Xvalidate, Yvalidate = Yvalidate
                , Xtest = reg$Xtest, Ytest = Ytest
                , Ytrain_reg = reg$Ytrain
                , Yvalidate_reg = reg$Yvalidate
                , Ytest_reg = reg$Ytest
                , U = reg$U, V = reg$V
                , lambda = reg$lambda, theta = reg$theta
                , beta = reg$beta
                , sigma0 = reg$sigma0, sigma1 = reg$sigma1)
  return(output)
}
