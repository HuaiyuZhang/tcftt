#################################################################
# The simulation code for TCFU and TTi tests
# Created at 9/18/2017 by Huaiyu
# The code borrows old code from folder /Zhou/2nd order/, 
# the notation is consistent with the writing now.
#
# Lastest change on 1/25/2020
# Add bootstrap testing functions.


################ CDF and Cornish-Fisher expansion ###############
# The CDF under Ha and 2nd-order Cornish-Fisher expansion
# Input: the population moments.
# Return: CDF under Ha and Cornish-Fisher under H0. list(cdf_H11, CFqH0) 
CDF_and_CFq <-  function(n1=20, n2=30, mu1=0, mu2=0, sigma1=1, sigma2=0.5, 
                         gamma1=1, gamma2=0, tau1=60, tau2=0)
{
  # Lower-level parameters
  n = n1+n2
  lambda1 = n1/n
  lambda2 = n2/n
  d1 = sqrt(lambda2*sigma1^2/(lambda2*sigma1^2 + lambda1*sigma2^2))
  d2 = sqrt(lambda1*sigma2^2/(lambda2*sigma1^2 + lambda1*sigma2^2))
  d3 = d1^2/sqrt(lambda1)
  d4 = d2^2/sqrt(lambda2)
  delta = mu1 - mu2
  w = delta/sqrt(sigma1^2/n1 +sigma2^2/n2)
  
  # High-level parameters
  a1 = -gamma1*d1*d3/2 + gamma2*d2*d4/2
  a2 = w*a1   
  a3 = a1*2/3
  b1 = 3*w/8*( d3^2*(tau1+2) + d4^2*(tau2+2) )
  b2 = ( gamma1*d1*d3 - gamma2*d2*d4 )^2 + (d3^2+d4^2) + w^2/8*( d3^2*(tau1+2) + d4^2*(tau2+2) )
  b3 = w*( 7/8*( gamma1*d1*d3 - gamma2*d2*d4 )^2 + (d3^2+d4^2)/2 )
  b4 = -(tau1*d3^2+tau2*d4^2)/12 + 2/3*(gamma1*d1*d3 - gamma2*d2*d4)^2 + (d3^2+d4^2)/2 + a2^2/2
  b5 = a2*a3 
  b6 = a3^2/2
  
  # CDF under Ha. Second order Edgeworth expansion
  cdf_Ha <-  function(x, type = 2)
  {
    H1 = x-w
    H2 = (x-w)^2 - 1
    H3 = (x-w)^3 - 3*(x-w)
    H4 = (x-w)^4 - 6*(x-w)^2 + 3;   
    H5 = (x-w)^5-10*(x-w)^3+15*(x-w)
    p1 = -( a1 + a2*H1 + a3*H2 )
    p2 = -( b1 + b2*H1 + b3*H2 + b4*H3 + b5*H4 + b6*H5 )
    if (type == 1){
      return(pnorm(x-w) + p1*dnorm(x-w)/sqrt(n))
    }
    if (type == 2){
      return(pnorm(x-w) + p1*dnorm(x-w)/sqrt(n) + p2*dnorm(x-w)/n)
    }
    
  }
  
  # Theoretical Cornish-Fisher expansion 
  # ONLY VALID UNDER H0: mu1=mu2, delta=w=0
  # if cdf is pnorm(x) + p1(x)dnorm(x)/sqrt(n) + p2(x)dnorm(x)/n,
  # p11(x) = -p1(x); p21(x) = p1(x)p1'(x) - 1/2*x*p1(x)^2 - p2(x)
  # Input: alpha
  # Return: Cornish-Fisher expansion at alpha
  CFqH0  <-  function(alpha, type = 2)
  {
    z = qnorm(alpha)
    H1z = z
    H2z = z^2 - 1
    H3z = z^3 - 3*z
    H4z = z^4 - 6*z^2 + 3   
    H5z = z^5 - 10*z^3 + 15*z
    p1z = -( a1 + a2*H1z + a3*H2z )
    p1_d = -( a2 + a3*z*2 )
    p2z = -( b1 + b2*H1z + b3*H2z + b4*H3z + b5*H4z + b6*H5z )
    p11 = - p1z                        
    p21 = p1z*p1_d - 1/2*z*p1z^2 - p2z
    if (type == 1)
    {
      z + p11/sqrt(n) 
    }
    else if (type == 2)
    {
      z + p11/sqrt(n) + p21/n   
    }
    else
    {
      stop("Undefined type.")
    }
    
  }
  # Return the cdf under Ha, and CF under H0
  list(cdf_Ha = cdf_Ha,  CFqH0 = CFqH0)
}

# Test the function
# tr1=CDF_and_CFq( n1=30, n2=30, mu1=0, mu2=0, sigma1=1, sigma2=0.5, gamma1=5, gamma2=0, tau1=3, tau2=3)
# curve(tr1$cdf_Ha(x),from = -5, to =5)
# tr1$cdf_Ha(1.96)
# tr1$cdf_Ha(-1.96)
# tr1$CFqH0(0.05)
# tr1$CFqH0(0.95)

################ TCFU test ########################
# The TCFU test function
# type = c(1, 2) specifies the TCFU1 or TCFU2
# alternative = c('less', 'greater', 'two.sided')
# Return: The testing result for one set of data
TCFU <- function(x1, x2, effectSize = 0, alternative = 'greater', 
                 alpha = 0.05, type = 2)
{
  n1 = length(x1)
  n2 = length(x2)  
  n = n1+n2
  # Estimate the moments and compute the test statistic
  mean1 = drop(matrix(rep(1, n1), nrow=1) %*% matrix(x1)/n1)
  mean2 = drop(matrix(rep(1, n2), nrow=1) %*% matrix(x2)/n2)
  var1 = drop(matrix(rep(1, n1), nrow=1) %*% matrix((x1-mean1)^2)/(n1-1))
  var2 = drop(matrix(rep(1, n2), nrow=1) %*% matrix((x2-mean2)^2)/(n2-1))
  Sp = sqrt( var1/n1+var2/n2 ) 
  y1 = (x1-mean1)/sqrt(var1)
  y2 = (x2-mean2)/sqrt(var2)
  gamma1 = (n1/((n1-1)*(n1-2))) * sum(y1^3)
  gamma2 = (n2/((n2-1)*(n2-2))) * sum(y2^3)
  tau1 = n1* rep(1, n1) %*% ((x1-mean1)^4)/ ((n1-1)*var1)^2 +6/(n1+1)
  tau2 = n2* rep(1, n2) %*% ((x2-mean2)^4)/ ((n2-1)*var2)^2 +6/(n2+1)
  t.stat = (mean1 - mean2 - effectSize)/Sp
  CFq = CDF_and_CFq(n1=n1, n2=n2, mu1=0, mu2=0, sigma1=sqrt(var1), sigma2=sqrt(var2), 
                    gamma1=gamma1, gamma2= gamma2, tau1=tau1, tau2=tau2)$CFqH0
  
  cdf = CDF_and_CFq(n1=n1, n2=n2, mu1=0, mu2=0, sigma1=sqrt(var1), sigma2=sqrt(var2), 
                    gamma1=gamma1, gamma2= gamma2, tau1=tau1, tau2=tau2)$cdf_Ha
  # pvalue = 1-cdf(t.stat)
  # Return the testing results
  if (alternative == 'greater')
  { 
    cutoff <- CFq(1-alpha, type)
    reject <- (t.stat >= cutoff)
    list(stat = t.stat, cutoff = cutoff, reject = reject)
  }
  else if (alternative=='less')
  {
    cutoff <- CFq(alpha, type)
    reject <-  (t.stat <= cutoff)
    list(stat = t.stat, cutoff = cutoff, reject = reject)
  } 
  else if (alternative=='two.sided')
  {
    reject <- (t.stat >= CFq(1-alpha/2, type) || t.stat <= CFq(alpha/2, type))
    list(stat = t.stat, cutoff = CFq(1-alpha/2, type),
         pvalue = 1-cdf(abs(t.stat), type) + cdf(-abs(t.stat), type), 
         reject = reject)
  }
  else
  {
    stop("Undefined alternative.")
  }
}

################ Transformation test ################
TTi <- function(x1, x2, alternative='greater', alpha=0.05, type = 1)
{
  n1 = length(x1)
  n2 = length(x2)  
  n = n1+n2
  lambda1 = n1/n
  lambda2 = n2/n
  # Estimate the moments and compute the test statistic
  mean1 = drop(matrix(rep(1, n1), nrow=1) %*% matrix(x1)/n1)
  mean2 = drop(matrix(rep(1, n2), nrow=1) %*% matrix(x2)/n2)
  xbar.diff = mean1 - mean2
  var1 = drop(matrix(rep(1, n1), nrow=1) %*% matrix((x1-mean1)^2)/(n1-1))
  var2 = drop(matrix(rep(1, n2), nrow=1) %*% matrix((x2-mean2)^2)/(n2-1))
  s1 = sqrt(var1)
  s2 = sqrt(var2)
  Sp = sqrt( var1/n1+var2/n2 ) 
  y1 = (x1-mean1)/sqrt(var1)
  y2 = (x2-mean2)/sqrt(var2)
  gamma1 = (n1/((n1-1)*(n1-2))) * sum(y1^3)
  gamma2 = (n2/((n2-1)*(n2-2))) * sum(y2^3)
  tau1 = n1* rep(1, n1) %*% ((x1-mean1)^4)/ ((n1-1)*var1)^2 +6/(n1+1)
  tau2 = n2* rep(1, n2) %*% ((x2-mean2)^4)/ ((n2-1)*var2)^2 +6/(n2+1)
  A_hat = ( s1^3*gamma1/(lambda1^2)- s2^3*gamma2/(lambda2^2) ) / ( s1^2/lambda1 + s2^2/lambda2 )^(3/2)
  # Transformations
  if (type == 1)
  {
    TransFun <- function(x) 
    {
      u = x/sqrt(n)
      sqrt(n)*( u + A_hat*u^2/3 + A_hat^2*u^3/27 + A_hat/(6*n) )
    }
  }
  else if (type == 2)
  {
    TransFun <- function(x) 
    {
      u = x/sqrt(n)
      sqrt(n)*( (2/3/sqrt(n)*A_hat )^(-1)*(exp(2/3/sqrt(n)*A_hat*u)-1) + A_hat/6/n )
      
    }
  }
  else if (type == 3)
  {
    TransFun <- function(x) 
    {
      u = x/sqrt(n)
      sqrt(n)*( u + u^2 + u^3/3 +  A_hat/6/n )
    }
  }
  else if (type == 4)
  {
    TransFun <- function(x) 
    {
      u = x/sqrt(n)
      sqrt(n)*( (2/3*A_hat )^(-1)*(exp(2/3*A_hat*u)-1) + A_hat/6/n ) #TT4
    }
  }
  else
  {
    stop("Undefined transformation type")
  }
  t.stat = xbar.diff/Sp
  transformedStat = TransFun(t.stat)
  # Testing result
  if (alternative == 'greater')
  { 
    reject <-  (transformedStat >= qnorm(1-alpha))
    list(reject=reject, stat=transformedStat, pvalue = 1-pnorm(transformedStat))
  }
  else if (alternative == 'less')
  {
    reject <- (transformedStat <= qnorm(alpha))
    list(reject=reject, stat=transformedStat, pvalue = pnorm(transformedStat))
    
  } 
  else if (alternative == 'two.sided')
  {
    reject <- (transformedStat >= qnorm(1-alpha/2) || transformedStat <= qnorm(alpha/2))
    list(reject=reject, stat=transformedStat, 
         pvalue = 2*(1-pnorm(abs(transformedStat))))
  }
  else
  {
    stop("Undefined alternative.")
  }
  
}

# Compute the estimate for A in the sample
computeAB <- function(x1, x2, computeB = F)
{
  n1 = length(x1)
  n2 = length(x2)  
  n = n1+n2
  lambda1 = n1/n
  lambda2 = n2/n
  # Estimate the moments and compute the test statistic
  mean1 = drop(matrix(rep(1, n1), nrow=1) %*% matrix(x1)/n1)
  mean2 = drop( matrix(rep(1, n2), nrow=1) %*% matrix(x2)/n2 )
  var1 = drop(matrix(rep(1, n1), nrow=1) %*% matrix((x1-mean1)^2)/(n1-1))
  var2 = drop(matrix(rep(1, n2), nrow=1) %*% matrix((x2-mean2)^2)/(n2-1))
  s1 = sqrt(var1)
  s2 = sqrt(var2)
  Sp = sqrt( var1/n1+var2/n2 ) 
  y1 = (x1-mean1)/sqrt(var1)
  y2 = (x2-mean2)/sqrt(var2)
  gamma1 = (n1/((n1-1)*(n1-2))) * sum(y1^3)
  gamma2 = (n2/((n2-1)*(n2-2))) * sum(y2^3)
  tau1 = n1* rep(1, n1) %*% ((x1-mean1)^4)/ ((n1-1)*var1)^2 +6/(n1+1)
  tau2 = n2* rep(1, n2) %*% ((x2-mean2)^4)/ ((n2-1)*var2)^2 +6/(n2+1)
  if (computeB){
    ( s1^4*(tau1-3)/(lambda1^3)- s2^4*(tau2-3)/(lambda2^3) ) / ( s1^2/lambda1 + s2^2/lambda2 )^2
    
  } 
  else{
    ( s1^3*gamma1/(lambda1^2)- s2^3*gamma2/(lambda2^2) ) / ( s1^2/lambda1 + s2^2/lambda2 )^(3/2)
    
  }
}


# Define the two-sample t-statistic using vectorization
t_stat <- function(x1,x2){
  n1 = length(x1)
  n2 = length(x2)  
  mean1 = drop(matrix(rep(1, n1), nrow=1) %*% matrix(x1)/n1)
  mean2 = drop(matrix(rep(1, n2), nrow=1) %*% matrix(x2)/n2)
  var1 = drop(matrix(rep(1, n1), nrow=1) %*% matrix((x1-mean1)^2)/(n1-1))
  var2 = drop(matrix(rep(1, n2), nrow=1) %*% matrix((x2-mean2)^2)/(n2-1))
  Sp = sqrt( var1/n1+var2/n2 ) 
  (mean1 - mean2)/Sp
}

# Bootstrap two-sample t-test
# input: data1, data2, number of bootstrap runs, alternative
# output: p-value
bootstrap_test <- function(x1, x2, B = 1000, alternative = 'greater'){
  computed_stat <- t_stat(x1, x2)
  t.vect <- rep(NA, B)
  data_combine <- c(x1, x2)
  n <- length(data_combine)
  n1 <- length(x1)
  n2 <- length(x2)
  for(i in 1:B){
    boot.1 <- sample(1:n, n1, replace=T)
    # boot.2 <- sample(1:n, n2, replace=T)
    boot.2 <- setdiff(1:n, boot.1)
    t.vect[i] <- t_stat(data_combine[boot.1], data_combine[boot.2])
  }
  if(alternative == 'greater'){
    pvalue = mean(t.vect > computed_stat)
  }
  else if (alternative == 'less'){
    pvalue = mean(t.vect < computed_stat)
  }
  else{
    pvalue = mean((t.vect < -abs(computed_stat) | t.vect > abs(computed_stat) ))
  }
  return(pvalue)
}
# x1 <- rnorm(100)
# x2 <- rnorm(40, mean = -0.2, 1)
# bootstrap_test(x1, x2, alternative = 'greater')


# -------------------------------------------------------
# size adjust pauc
pauc <- function(stat_h0, stat_ha, 
                         target_range_lower = 0.01, 
                         target_range_upper = 0.2){
    flag <- stat_h0 > quantile(stat_h0,probs = target_range_lower) &
      stat_h0 < quantile(stat_h0,probs = target_range_upper)
    stat_h0 <- stat_h0[flag]
    rank_sum <- length(stat_h0) * length(stat_ha)
    adj_power <- wilcox.test(stat_ha, stat_h0)$statistic/rank_sum
    return(adj_power)
}