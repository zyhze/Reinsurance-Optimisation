# install.packages("fitdistrplus")
# install.packages("ADGofTest")
# install.packages("goftest")
# install.packages("actuar")
# install.packages("stats")
library(fitdistrplus)
library(ADGofTest)
library(goftest)
library(actuar)
library(stats)

LOSS <- read.csv("C:\\Users\\Leon\\Downloads\\Loss.csv")
LOSS <- as.numeric(LOSS$X11515.2895389355)

# Question 1
# Fitting Log-Normal Distribution
logNormParams <- fitdist(data = LOSS, distr = "lnorm")

# Fitting Exponential Distribution
expParams <- fitdist(data = LOSS/100, distr = "exp")

# Fitting Gamma Distribution
gammaParams <- fitdist(data = LOSS, distr = "gamma", lower = c(0, 0))

# Fitting Pareto Distribution
# Estimated parameters from MLE
paretoLambda = min(LOSS)
paretoAlpha = length(LOSS) / sum(log(LOSS) - log(paretoLambda))


# Question 2
# QQ Plots
# Log-Normal
par(mfrow = c(2, 2), oma = c(0, 0, 2, 0)) 

quants <- seq(0, 1, length = 81)
LOSSQuants <- quantile(LOSS, quants)

logQuants <- qlnorm(quants, logNormParams$estimate['meanlog'], logNormParams$estimate['sdlog'])

plot(logQuants, LOSSQuants, xlab="theoretical quantiles",ylab="sample quantiles",
     main="Log-Normal",cex=0.9, xlim = c(0, 50000), ylim = c(0, 50000))
abline(0,1,col="red")

# Exponential
expQuants <- qexp(quants, expParams$estimate['rate']/100)
plot(expQuants, LOSSQuants, xlab="theoretical quantiles",ylab="sample quantiles",
     main="Exponential",cex=0.9, xlim = c(0, 30000), ylim = c(0, 30000))
abline(0,1,col="red")

# Gamma
gammaQuants <- qgamma(quants, gammaParams$estimate['shape'], gammaParams$estimate['rate'])
plot(gammaQuants, LOSSQuants, xlab="theoretical quantiles",ylab="sample quantiles",
     main="Gamma",cex=0.9, xlim = c(0, 35000), ylim = c(0, 35000))
abline(0,1,col="red")

# Pareto
qpareto<-function(y,alpha,lambda){ 
  rep(lambda,length(y))/(rep(1,length(y))-y)^(1/alpha)
}

empirical<-ecdf(LOSS)
plot(qpareto(empirical(LOSS),paretoAlpha,paretoLambda),LOSS,
     xlab="theoretical quantiles",ylab="sample quantiles",main="Pareto",
     cex=0.9, xlim = c(0, 10000), ylim = c(0, 10000))
abline(0,1,col="red")

mtext("Figure 2.3: QQ Plots", outer = TRUE, cex = 1.5, font = 2)
par(mfrow = c(1, 1)) 

# Histograms
hist(LOSS, breaks = 100, prob = T, xlab = "Claims", main = "Figure 2.1: Histogram versus Densities", ylim = c(0, 0.0002))
curve(dlnorm(x, meanlog = logNormParams$estimate['meanlog'], sdlog = logNormParams$estimate['sdlog']), 
      col = "forestgreen", lwd = 2.5, add = TRUE)
curve(dexp(x, rate = expParams$estimate['rate']/100), col = "orange", lwd = 2.5, add = TRUE)
curve(dgamma(x, shape = gammaParams$estimate['shape'], rate = gammaParams$estimate['rate']), col = "red2", lwd = 2.5, add = TRUE)
curve(dpareto(x, scale = paretoLambda, shape = paretoAlpha), col = "dodgerblue", lwd = 2.5, add = TRUE) 
legend("topright", legend = c("Log-normal", "Exponential", "Gamma", "Pareto"), 
       col = c("forestgreen", "orange", "red2", "dodgerblue"), lwd = 2)

# ECDF Plots
plot(ecdf(LOSS), main = "Figure 2.2: ECDF", 
     xlab = "Claims", ylab = "Cumulative Probability", col = "black", 
     lwd = 2, xlim = c(0, max(LOSS)), pch = 16, cex = 1.2, 
     verticals = FALSE, do.points = TRUE)
x_vals <- seq(min(LOSS), max(LOSS), length.out = 1000)
lines(x_vals, plnorm(x_vals, meanlog = logNormParams$estimate['meanlog'], sdlog = logNormParams$estimate['sdlog']), 
      col = "forestgreen", lwd = 2)
lines(x_vals, pexp(x_vals, rate = expParams$estimate['rate']/100), 
      col = "orange", lwd = 2)
lines(x_vals, pgamma(x_vals, shape = gammaParams$estimate['shape'], rate = gammaParams$estimate['rate']), 
      col = "red2", lwd = 2)
lines(x_vals, ppareto(x_vals, scale = paretoLambda, shape = paretoAlpha), 
      col = "dodgerblue", lwd = 2)
legend("bottomright", legend = c("Data", "Log-normal", "Exponential", "Gamma", "Pareto"), 
       col = c("black","forestgreen", "orange", "red2", "dodgerblue"), lwd = 2)
par(mfrow = c(1, 1)) 
# Hypothesis Tests
# Log-Normal
logNormKS <- ks.test(LOSS, "plnorm", meanlog = logNormParams$estimate['meanlog'],
                     sdlog = logNormParams$estimate['sdlog'])

logNormAD <- ad.test(LOSS, "plnorm", meanlog = logNormParams$estimate['meanlog'],
                     sdlog = logNormParams$estimate['sdlog'])

# Exponential
expKS <- ks.test(LOSS, "pexp", rate = expParams$estimate["rate"]/100)

expAD <- ad.test(LOSS, "pexp", rate = expParams$estimate["rate"]/100)

# Gamma
gammaKS <- ks.test(LOSS, "pgamma", rate = gammaParams$estimate["rate"],
                 shape = gammaParams$estimate["shape"])
gammaAD <- ad.test(LOSS, "pgamma", rate = gammaParams$estimate["rate"],
                   shape = gammaParams$estimate["shape"])

# Pareto
paretoKS <- ks.test(LOSS, "ppareto", shape = paretoAlpha, scale = paretoLambda)

paretoAD <- ad.test(LOSS, "ppareto", shape = paretoAlpha, scale = paretoLambda)


# AIC, BIC
n <- length(LOSS)

logNormLL <- -sum(log(LOSS)) - n * log(logNormParams$estimate['sdlog']) - 
  (n / 2) * log(2 * pi) - 
  sum((log(LOSS) - logNormParams$estimate['meanlog'])^2 / (2 * logNormParams$estimate['sdlog']^2))

logNormAIC <- -2 * logNormLL + 2 * 2
logNormBIC <- -2 * logNormLL+ log(n) * 2

expLL <- n * log(expParams$estimate["rate"]/100) - expParams$estimate["rate"]/100 * sum(LOSS)
expAIC <- -2 * expLL + 2 * 1
expBIC <- -2 * expLL+ log(n) 

gammaLL <- n * gammaParams$estimate["shape"] * log(gammaParams$estimate["rate"]) - 
  n * log(gamma(gammaParams$estimate["shape"])) + (gammaParams$estimate["shape"] - 1) * sum(log(LOSS)) - 
  gammaParams$estimate["rate"] * sum(LOSS)

gammaAIC <- -2 * gammaLL + 2 * 2
gammaBIC <- -2 * gammaLL + log(n) * 2

paretoLL <- n * log(paretoAlpha) + n * paretoAlpha * log(paretoLambda) - 
              (paretoAlpha + 1) * sum(log(LOSS))
paretoAIC <- -2 * paretoLL + 2 * 2
paretoBIC <- -2 * paretoLL + log(n) * 2

# Censored Data
# Question 1
# Exponential
k <- sum(LOSS < 10000)
sum.k <- sum(LOSS[LOSS < 10000])
policy.limit <- 10000

censoredExpParam <- k / (sum.k + policy.limit * (n - k))

# Pareto
censLOSS <- LOSS
censLOSS[censLOSS > 10000] <- 10000
filterLOSS <- LOSS[LOSS <= 10000]

censParetoLambda = min(censLOSS)
censParetoAlpha = k / (sum(log(filterLOSS)) - k * log(censParetoLambda) - (n - k) * log(censParetoLambda/policy.limit))

# QQ Plots
# Exponential
par(mfrow = c(1, 2), oma = c(0, 0, 2, 0)) 
quants <- seq(0, 1, length = 81)
censLOSSQuants <- quantile(LOSS, quants)

expCensQuants <- qexp(quants, censoredExpParam)
plot(expCensQuants, censLOSSQuants, xlab="theoretical quantiles",ylab="sample quantiles",
     main="Exponential",cex=0.9, xlim = c(0, 30000), ylim = c(0, 30000))
abline(0,1,col="red")

# Pareto
empiricalCens<-ecdf(LOSS)
plot(qpareto(empiricalCens(LOSS),censParetoAlpha,censParetoLambda),LOSS,
     xlab="theoretical quantiles",ylab="sample quantiles",main="Pareto",
     cex=0.9, xlim = c(0, 10000), ylim = c(0, 10000))
abline(0,1,col="red")
mtext("Figure 2.4: QQ Plots for Censored data", outer = TRUE, cex = 1.5, font = 2)



