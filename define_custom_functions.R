# Compute the negative log-likelihood of exponentially distributed data
NegLogLikExp <- function(x, rate){
  
  neg_log_lik <- -sum(log(rate) - rate*x) 
  return(neg_log_lik)

  }

# Compute the gradient of the negative log-likelihood function of exponentially distributed data
GradNegLogLikExp <- function(x, rate){
  
  grad_neg_log_lik <- -length(x)/rate + sum(x)
  return(grad_neg_log_lik)

  }

# Compute the negative log-likelihood of data with log-normal distribution
NegLogLikLogNorm <- function(x, par){
  
  meanlog <- par[1]
  sdlog <- par[2]
  
  neg_log_lik <- -sum(dlnorm(meanlog, sdlog, log=TRUE)) 
  return(neg_log_lik)
  
}

# Compute the negative log-likelihood of right-censored, left-truncated data with exponential distribution
NegLogLikCensTruncExp <- function(par, y, c, t_y, t_c){

	# par: parameter of the exponential distribution cf. ANLIM L3 slide 11
	# y: observations with exact values
	# c: observations with censored values
	# t_y,t_c: the observations' corresponding truncation treshold

  # The general form of the log-likelihood function for independent right-censored and left-truncated samples 
  #{(y1,c1,delta1),...,(y1,c1,delta1)} is given on slide 24 of ANLIM L3
  
  rate <- par
  E <- length(y)
  
  neg_log_lik <- -sum(dexp(y, rate, log=TRUE)) + sum(log(1-pexp(t_y, rate))) - sum(log(1-pexp(c, rate))) + sum(log(1-pexp(t_c,rate)))
  
  # neg_log_lik_alt <- -E*log(rate) + rate*(sum(y) - sum(t_y) + sum(c) - sum(t_c))  
  # out <- c(neg_log_lik,neg_log_lik_alt)
  
  return(neg_log_lik)
  
}

# Compute the negative log-likelihood of right-censored, left-truncated data with log-normal distribution
NegLogLikCensTruncLogNorm <- function(par, y, c, t_y, t_c){
  
  # par: parameters of the log-normal distribution cf. ANLIM L3 slide 11
  # y: observations with exact values
  # c: observations with censored values
  # t_y,t_c: the observations' corresponding truncation treshold
  
  # The general form of the log-likelihood function for independent right-censored and left-truncated samples 
  #{(y1,c1,delta1),...,(y1,c1,delta1)} is given on slide 24 of ANLIM L3
  
  meanlog <- par[1]
  sdlog <- par[2]
  
  neg_log_lik <- -sum(dlnorm(y, meanlog, sdlog, log=TRUE)) + sum(log(1-plnorm(t_y, meanlog, sdlog))) - sum(log(1-plnorm(c, meanlog, sdlog))) + sum(log(1-plnorm(t_c, meanlog, sdlog)))
  
  return(neg_log_lik)
  
}

# Compute the negative log-likelihood of right-censored, left-truncated data with inverse Gaussian distribution
NegLogLikCensTruncInvGauss <- function(par, y, c, t_y, t_c){
  
  # par: parameters of the inverse Gaussian distribution cf. ANLIM L3 slide 11
  # y: observations with exact values
  # c: observations with censored values
  # t_y,t_c: the observations' corresponding truncation treshold
  
  # The general form of the log-likelihood function for independent right-censored and left-truncated samples 
  #{(y1,c1,delta1),...,(y1,c1,delta1)} is given on slide 24 of ANLIM L3
  
  mu <- par[1]
  lambda <- par[2]
  
  neg_log_lik <- -sum(dinvgauss(y, mu, lambda, log=TRUE)) + sum(log(1-pinvgauss(t_y, mu, lambda))) - sum(log(1-pinvgauss(c, mu, lambda))) + sum(log(1-pinvgauss(t_c, mu, lambda)))
  
  return(neg_log_lik)
  
}

# Compute the negative log-likelihood of right-censored, left-truncated data with Burr distribution
NegLogLikCensTruncBurr <- function(par, y, c, t_y, t_c){
  
  # par: parameters of the inverse Gaussian distribution cf. ANLIM L3 slide 11
  # y: observations with exact values
  # c: observations with censored values
  # t_y,t_c: the observations' corresponding truncation treshold
  
  # The general form of the log-likelihood function for independent right-censored and left-truncated samples 
  #{(y1,c1,delta1),...,(y1,c1,delta1)} is given on slide 24 of ANLIM L3
  
  shape1 <- par[1]
  shape2 <- par[2]
  rate <- par[3]
  
  neg_log_lik <- -sum(dburr(y, shape1, shape2, rate, log=TRUE)) + sum(log(1-pburr(t_y, shape1, shape2, rate))) - sum(log(1-pburr(c, shape1, shape2, rate))) + sum(log(1-pburr(t_c, shape1, shape2, rate)))
  
  return(neg_log_lik)
  
}

# Compute the negative log-likelihood of right-censored, left-truncated data with Erlang Mixture distribution
NegLogLikCensTruncEM <- function(par, y, c, t_y, t_c){
  
  # par: parameters of the inverse Gaussian distribution cf. ANLIM L3 slide 11
  # y: observations with exact values
  # c: observations with censored values
  # t_y,t_c: the observations' corresponding truncation treshold
  
  # The general form of the log-likelihood function for independent right-censored and left-truncated samples 
  #{(y1,c1,delta1),...,(y1,c1,delta1)} is given on slide 24 of ANLIM L3
  
  theta <- par$theta
  shape <- par$shape
  alpha <- par$alpha
  
  neg_log_lik <- -sum(ME_density(y, theta, shape, alpha, trunclower=0, truncupper=Inf, log=TRUE)) + sum(log(1-ME_cdf(t_y, theta, shape, alpha, trunclower=0, truncupper=Inf, lower.tail=TRUE, log.p=FALSE))) - sum(log(1-ME_cdf(c, theta, shape, alpha, trunclower=0, truncupper=Inf, lower.tail=TRUE, log.p=FALSE))) + sum(log(1-ME_cdf(t_c, theta, shape, alpha, trunclower=0, truncupper=Inf, lower.tail=TRUE, log.p=FALSE)))
  
  return(neg_log_lik)
  
}


# Compute gradient of the negative log-likelihood of right-censored, left-truncated data with exponential distribution
GradNegLogLikCensTruncExp <- function(par, y, c, t_y, t_c){
  
  # par: parameter of the exponential distribution cf. ANLIM L3 slide 11
  # y: observations with exact values
  # c: observations with censored values
  # t_y,t_c: the observations' corresponding truncation treshold
  
  # The general form of the log-likelihood function for independent right-censored and left-truncated samples 
  #{(y1,c1,delta1),...,(y1,c1,delta1)} is given on slide 24 of ANLIM L3
  
  rate <- par
  E <- length(y)
  
  # neg_log_lik <- sum(log(dexp(y,rate))) + sum(log(1-pexp(t_y,rate))) - sum(log(1-pexp(c,rate))) + sum(log(1-pexp(t_c)))
  # the above computation becomes numerically unreliable in case rate->0, resulting in neg_log_lik =  -Inf or NaN
  grad_neg_log_lik <- -E*1/rate + sum(y) - sum(t_y) + sum(c) - sum(t_y)  
  
  return(grad_neg_log_lik)
  
}

# Compute the value of the Akaike information criterion cf. ANLIM L2 slide 79
AICriterion <- function(ll_max, par_dim){
  
  	# ll_max: maximum log-likelihood
  	# par_dim: dimension of the parameter vector
  	
  aic_val <- -2*ll_max + 2*par_dim
  return(aic_val)

}

PdfSplice <- function(par, x){
  
  lambda <- par[1]
  x_shift <- par[2]
  n <- par[3]
  k <- par[4]
  X_treshold <- par[5]
  
  x_exp <- x[x <= X_treshold]
  x_par <- x[x > X_treshold]
  
  splice_exp <- (n-k)/n*(lambda*exp(-lambda*(x_exp-x_shift))/(1-exp(-lambda*(X_treshold-x_shift))))
  splice_par <- (k/n)*(alpha*(x+1)^(-alpha-1))/((X_treshold+1)^-alpha)
  pdf_splice <- c(splice_exp, splice_par)
  
  return(pdf_splice)
  
}

CdfSplice <- function(par, x){
  
  lambda <- par[1]
  alpha <- par[2]
  x_shift <- par[3]
  n <- par[4]
  k <- par[5]
  X_treshold <- par[6]
  
  x_exp <- x[x <= X_treshold]
  x_par <- x[x > X_treshold]
  
  splice_exp <- (n-k)/n*pexp(x_exp-x_shift, lambda)/pexp(X_treshold-x_shift,lambda)
  splice_par <- 1-(k/n)*((x_par+1)/(X_treshold+1))^-alpha
  splice_exp <- (n-k)/n*(1-exp(-lambda*(x_exp-x_shift)))/(1-exp(-lambda*(X_treshold-x_shift)))
  
  x <- c(x_exp, x_par)
  value <- c(splice_exp,splice_par)
  cdf_splice <- data.frame(x, value)
  
  return(cdf_splice)
  
}

# Compute the negative log-li
NegLogLikSplice <- function(par, X, X_shift, k){
  
  lambda <- par[1]
  alpha <- par[2]
  
  
  X <- X[order(X)]
  n <- length(X)
  X_treshold <- X[n-k]
  
  X_exp <- X[X <= X_treshold]
  X_par <- X[X > X_treshold]
  
  N_e <- length(X_exp)
  N_p <- length(X_par)
  
  neg_log_lik_exp <- -sum(dexp(X_exp-X_shift, lambda, log=TRUE)) - N_e*(log((n-k)/n)-pexp(X_treshold-X_shift,lambda, log=TRUE))
  
  neg_log_lik_par <- -N_p*(log(k/n)+log(alpha)+alpha*log(X_treshold+1)) + (alpha+1)*sum(log(X_par+1))
  neg_log_lik <-neg_log_lik_exp +  neg_log_lik_par
  
  return(neg_log_lik)
}

