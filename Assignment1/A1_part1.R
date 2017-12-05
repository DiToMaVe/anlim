# Set working directory
#######################

setwd("C:/Users/Gebruiker/Documents/ANLIM/Assignments/A1")

# Load libraries
################
library(statmod)
library(actuar)
library(survival)
library(ggplot2)
library(plotly)
# Sys.setenv("plotly_username"="")
# Sys.setenv("plotly_api_key"="")


# Load data set
###############

severity_censoring = read.table("SeverityCensoring.txt")

# Load user-defined functions
#############################
source("define_custom_functions.R")

# Load EM algorithm for Erlang Mixtures
#######################################
source("2014-12-16_ME.R")

# Cast data to appropriate format, cf. ANLIM L3 slide 24
########################################################

# convert NA to Boolean
is_rc <- is.na(severity_censoring$rc) 
# see slide 21 of ANLIM L3 for the expression for delta
delta <- 0^is_rc; 

yy <- subset(severity_censoring$claimAmount, delta==0)
cc <- subset(severity_censoring$claimAmount, delta==1)

tt_y <- subset(severity_censoring$deductible, delta==0)
tt_c <- subset(severity_censoring$deductible, delta==1)

# ML estimate of exponential distribution
#########################################

# Initial parameter estimate with method of moments
y_exp <- severity_censoring$claimAmount
rate_init = 1/mean(y_exp);

# MLE 
optim_res <- optim(rate_init, NegLogLikCensTruncExp, y=yy, c=cc, t_y=tt_y, t_c=tt_c, gr=GradNegLogLikCensTruncExp, method="BFGS", hessian=FALSE, control=list(trace=2))
rate_cens_trunc_exp_mle <- optim_res$par
ll_cens_trunc_exp_mle <- -optim_res$value

# Akaike Information Criterion
AIC_cens_trunc_exp <- AICriterion(ll_cens_trunc_exp_mle,1)

# ML estimate of lognormal distribution
#######################################

# Initial parameter estimate with method of moments
y_inv_gauss <- severity_censoring$claimAmount
n <- length(y_inv_gauss)
meanlog_init <- 2*log(1/n*sum(y_inv_gauss)) - 1/2*log(1/n*sum(y_inv_gauss^2))
varlog_init <- log(1/n*sum(y_inv_gauss^2)) - 2*log(1/n*sum(y_inv_gauss)) # MoM results in a negative estimate of the variance varlog_init

# MLE
par_lnorm_init <- c(meanlog_init,varlog_init)
optim_res <- optim(par_lnorm_init, NegLogLikCensTruncLogNorm, y=yy, c=cc, t_y=tt_y, t_c=tt_c, method="Nelder-Mead", control=list(trace=2))
par_cens_trunc_lnorm_mle <- optim_res$par
ll_cens_trunc_lnorm_mle <- -optim_res$value

# Akaike Information Criterion
AIC_cens_trunc_lnorm <- AICriterion(ll_cens_trunc_lnorm_mle,2)

# ML estimate of inverse Gaussian distribution
##############################################

# Initial parameter estimate with method of moments
y_inv_gauss <- severity_censoring$claimAmount
mu_init <- mean(y_inv_gauss)
lambda_init <- mu_init^3/var(y_inv_gauss)
par_inv_gauss_init <- c(mu_init, lambda_init)

# MLE
optim_res <- optim(par_inv_gauss_init, NegLogLikCensTruncInvGauss, y=yy, c=cc, t_y=tt_y, t_c=tt_c, method="Nelder-Mead", control=list(trace=2))
par_cens_trunc_inv_gauss_mle <- optim_res$par
ll_cens_trunc_inv_gauss_mle <- -optim_res$value

# Akaike Information Criterion
AIC_cens_trunc_inv_gauss <- AICriterion(ll_cens_trunc_inv_gauss_mle,2)

# ML estimate of Burr type XII distribution
###########################################

# The method of moments and maximum likelihood estimators for the Burr distribution can only be evaluated numerically 
# To remediate the influence of the initialization different starting values were tested
par_burr_init <- c(1,1,1*10^-15)

# MLE
optim_res <- optim(par_burr_init, NegLogLikCensTruncBurr, y=yy, c=cc, t_y=tt_y, t_c=tt_c, method="Nelder-Mead", control=list(trace=2))
par_cens_trunc_burr_mle <- optim_res$par
ll_cens_trunc_burr_mle <- -optim_res$value

# Akaike Information Criterion
AIC_cens_trunc_burr <- AICriterion(ll_cens_trunc_burr_mle,3)

# ML estimate of Erlang Mixture with EM algorithm
#################################################

loss <- severity_censoring$claimAmount
nrc <- severity_censoring$claimAmount*(is.na(severity_censoring$rc)) 
nrc[nrc==0] <- NA
fit.ME <- ME_fit(loss, nrc, trunclower=100, M=5, s=3)

# Kaplan-Meier estimate of survival function
############################################
deds <- severity_censoring$deductible 
loss <- severity_censoring$claimAmount
full <- is.na(severity_censoring$rc)
fit <- survfit(Surv(deds, loss, full) ~ 1)
# Kaplan-Meier plot
plot(fit, mark.time=F, conf.int=F, xlab = 'Claim size', ylab = 'Survival')

#  Reshuffle data in order to have monotoneously increasing losses. This wil ensure easy plotting with lines()
df <- data.frame(loss,deds)
df_sorted <- df[order(df$loss),]

loss <- df_sorted$loss
deds <- df_sorted$deds

# Cumulative density functions
cdf_trunc_exp <- (1-pexp(loss, rate_cens_trunc_exp_mle))/(1-pexp(deds, rate_cens_trunc_exp_mle))
cdf_trunc_lnorm <- (1-plnorm(loss, par_cens_trunc_lnorm_mle[1], par_cens_trunc_lnorm_mle[2]))/(1-plnorm(deds,par_cens_trunc_lnorm_mle[1],par_cens_trunc_lnorm_mle[2]))
cdf_trunc_inv_gauss <- (1-pinvgauss(loss, par_cens_trunc_inv_gauss_mle[1],par_cens_trunc_inv_gauss_mle[2]))/(1-pinvgauss(deds, par_cens_trunc_inv_gauss_mle[1], par_cens_trunc_inv_gauss_mle[2]))
cdf_trunc_burr <- (1-pburr(loss, par_cens_trunc_burr_mle[1], par_cens_trunc_burr_mle[2], par_cens_trunc_burr_mle[3]))/(1-pburr(deds, par_cens_trunc_burr_mle[1], par_cens_trunc_burr_mle[2], par_cens_trunc_burr_mle[3]))
alpha <- fit.ME$alpha
shape <- fit.ME$shape
theta <- fit.ME$theta
cdf_trunc_me <- (1-ME_cdf(loss, theta, shape, alpha, trunclower=0, truncupper=Inf, lower.tail=TRUE))/(1-ME_cdf(deds, theta, shape, alpha, trunclower=0, truncupper=Inf, lower.tail=TRUE)) 

# Add fitted distributions to Kaplan-Meier plot
lines(loss, cdf_trunc_exp, col="blue", lwd=1)
lines(loss, cdf_trunc_lnorm, col="yellow", lwd=1)
lines(loss, cdf_trunc_inv_gauss, col="pink", lwd=1)
lines(loss, cdf_trunc_burr, col="green", lwd=1)
lines(loss, cdf_trunc_me, col="red", lwd=0.1)

# Fancier plot with plotly
##########################

# Store results in data frame
results_km_df <- data.frame(fit$time, fit$surv)
results_fit_df <- data.frame(loss, cdf_trunc_exp, cdf_trunc_lnorm, cdf_trunc_inv_gauss, cdf_trunc_burr, cdf_trunc_me)

# Plot
p <- plot_ly(results_fit_df, x=~loss) 
add_trace(p, y=~cdf_trunc_exp, name = "exponential", mode = "lines") %>%
add_trace(p, y=~cdf_trunc_lnorm, name = "log-normal", mode ="lines") %>%
add_trace(p, y=~cdf_trunc_inv_gauss, name = "inverse Gaussian", mode ="lines") %>%
add_trace(p, y=~cdf_trunc_burr, name = "Burr", mode = "lines") %>%
add_trace(p, y=~cdf_trunc_me, name = "Erlang mixture", mode = "lines") %>%
add_trace(p, data = results_km_df, x=~fit.time, y=~fit.surv, name = "Kaplan-Meier", mode = "lines")

# plotly_POST(p)

# AIC summary 
#############
values <- c(AIC_cens_trunc_exp,AIC_cens_trunc_lnorm,AIC_cens_trunc_inv_gauss,AIC_cens_trunc_burr,fit.ME$AIC)
labels <- c("exp","lnorm","inv_gauss","burr","em")
AIC_summary <- data.frame(values,labels)
AIC_summary <- AIC_summary[order(AIC_summary$values),]





