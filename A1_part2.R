# Set working directory
#######################

setwd("C:/Users/Gebruiker/Documents/ANLIM/Assignments/A1")

# Load user-defined functions
#############################
source("define_custom_functions.R")

# Load data set
###############

data = read.table("securaRe.txt", header = TRUE)
claim_amount <- data.frame(data)
claim_amount <- claim_amount$Loss

k_hat <- 95
X <- claim_amount
X_shift <- 1200000

X <- X[order(X)]
n <- length(X)
X_treshold <- X[n-k_hat]

X_exp <- X[X <= X_treshold]
X_par <- X[X > X_treshold]

lambda_init <- 1/mean(X_exp)
alpha_init <- mean(n/k_hat*(X_par+1))/(mean(n/k_hat*(X_par+1))-(X_treshold+1))
                   
par_splice_init <- c(0.000000000000001,alpha_init)

optim_res <- optim(par_splice_init, NegLogLikSplice, X=claim_amount, X_shift=X_shift, k=k_hat, method="Nelder-Mead", control=list(trace=2))
par_splice_mle <- optim_res$par
ll_splice_mle <- -optim_res$value

# Compare cdf and ecdf
######################

# Obtain empirical CDF values
claim_amount.ecdf = ecdf(claim_amount)
summary(claim_amount.ecdf)

# plot the ecdf
plot(claim_amount.ecdf, xlab = 'Sample Quantiles of Claims', ylab = '', main = 'Compare ecdf and cdf')
mtext(text = expression(hat(F)[n](x)), side = 2, line = 2.5)

# Obtain analytical cdf with MLE parameters
par_cdf <- c(par_splice_mle[1], par_splice_mle[2], X_shift, n, k_hat, X_treshold)
cdf_ana <- CdfSplice(par_cdf, claim_amount)

#  Reshuffle data in order to have monotoneously increasing losses. This wil ensure easy plotting with lines()
cdf_ana_sorted <- cdf_ana[order(cdf_ana$x),]
lines(cdf_ana_sorted$x, cdf_ana_sorted$value, col="red")







