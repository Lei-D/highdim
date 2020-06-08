# this function returns the best threshold t for a sorted list
find_t = function(list){
  # input list is a sorted list in ascending order
  
  if (max(list) == 0) {
    t = 0
  } else {
    p = length(list)
    
    # compute total variance when divide the list in position i
    total_var = rep(NA, p-1)
    for(i in 1:(p-1)){
      var1 = ifelse(i == 1, 0, var(list[1:i]))
      var2 = ifelse(i == p-1, 0, var(list[(i+1):p]))
      total_var[i] = i * var1 + (p-i) * var2
    }
    
    # compute first derivative of total variance
    first_Derivative = total_var[2:(p-1)] - total_var[1:(p-2)]
    
    # find the first point that moving it from the 2nd group to 1 group
    # increase the first derivative
    i = 2
    index = 1
    while (i < p - 2) {
      if (first_Derivative[i] - first_Derivative[i-1] > mean(abs(first_Derivative[1:(i-1)]))) {
        index = i
        break
      } else {
        i = i + 1
      }
    }
    t = list[index]
  }
  return(t)
}


# return values for lambda
compute_lambda = function(S, maxnvar, nsol, type){
  maxoffdiag = sort(compute_maxoffdiag(S), decreasing = T)
  lambda_max = maxoffdiag[1]
  lambda_min = maxoffdiag[maxnvar]
  # lambda_min = 0
  if (type == "linear") {
    lambda = seq(from = lambda_max, to = lambda_min, length.out = nsol)
  } else if (type == "loglinear") {
    lambda = lambda_min * log10(seq(1, 10, length.out = nsol)) + 
      lambda_max * (1 - log10(seq(1, 10, length.out = nsol)))
  }
  return(lambda)
}



# return the maximum off diagnal absolute value in each row 
compute_maxoffdiag = function(S){
  # S is the input matrix which is symmetric
  S = abs(S)
  n = nrow(S)
  out = rep(NA, n)
  out[1] = max(S[1, 2:n])
  for (i in 2:(n-1)) {
    out[i] = max(max(S[i, 1:(i-1)]), max(S[i, (i+1):n]))
  }
  out[n] = max(S[n, 1:(n-1)])
  return(out)
}
