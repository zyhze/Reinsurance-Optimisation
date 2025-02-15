
# Question 1
# Determining ruin probability
theta <- 0.06
lambda <-5

# Exponential
exp.rate <- 3.0597
exp.premium.rate <- (1 + theta) * lambda * (1 / exp.rate)
c0 <- 18
t <- 60
size <- 1000000

expRuin <- function() {
  N <- rpois(1, lambda * t)
  
  arrival.time <- sort(runif(N, 0, t))
  
  claim.size <- rexp(N, rate = exp.rate)
  
  surplus <- c0 + exp.premium.rate * arrival.time - cumsum(claim.size) 
  
  min(surplus) < 0
}
set.seed(2)
expSim <- replicate(size, expRuin())
expRuin <- mean(expSim)


# Gamma ruin probability
gamma.shape <- 0.67
gamma.rate <- 2.05

gamma.premium.rate <- (1 + theta) * lambda * (gamma.shape / gamma.rate)
gammaRuin <- function() {
  N <- rpois(1, lambda * t)
  
  arrival.time <- sort(runif(N, 0, t))
  
  claim.size <- rgamma(N, shape = gamma.shape, rate = gamma.rate)
  
  surplus <- c0 + gamma.premium.rate * arrival.time - cumsum(claim.size)
  
  min(surplus) < 0
}
set.seed(2)
gammaSim <- replicate(size, gammaRuin())
gammaRuin <- mean(gammaSim)

# Question 2 Part a
# Analysis of reinsurance products
# Reinsurance A
loading <- 0.08
retention <- 0.88

gamma.A.premium.rate <- (1 + loading) * lambda * (1 - retention) * (gamma.shape / gamma.rate)
ruinA <- function() {
  N <- rpois(1, lambda * t)
  
  arrival.time <- sort(runif(N, 0, t))
  
  claim.size <- rgamma(N, shape = gamma.shape, rate = gamma.rate)
  
  surplus <- c0 + (gamma.premium.rate - gamma.A.premium.rate) * arrival.time - retention * cumsum(claim.size)
  
  min(surplus) < 0
}
set.seed(2)
ruinSimA <- replicate(size, ruinA())
gammaRuinA <- mean(ruinSimA)

# Reinsurance B
retentionB <- 0.881391

exp.val <- gamma.shape/gamma.rate * (1 - pgamma(retentionB, 1 + gamma.shape, gamma.rate)) - 
  retentionB * (1 - pgamma(retentionB, gamma.shape, gamma.rate))
gamma.B.premium.rate <- (1 + loading) * lambda * exp.val

ruinB <- function() {
  N <- rpois(1, lambda * t)
  
  arrival.time <- sort(runif(N, 0, t))
  
  claim.size <- rgamma(N, shape = gamma.shape, rate = gamma.rate)
  
  claim.size <- pmin(claim.size, retentionB)
  
  surplus <- c0 + (gamma.premium.rate - gamma.B.premium.rate) * arrival.time - cumsum(claim.size)
  
  min(surplus) < 0
}
set.seed(2)
ruinSimB <- replicate(size, ruinB())
gammaRuinB <- mean(ruinSimB)

# Question 2 Part b
# Reinsurance A
alpha <- (loading - theta) / loading
rangeAlpha <- 1 - alpha

# Reinsurance B
claimIntegral <- function(d) {
  1 - pgamma(d, shape = gamma.shape, rate = gamma.rate)  
}

retentionRange <- function(B) {
  (1 + theta) * gamma.shape/gamma.rate - 
    (1 + loading) * (gamma.shape/gamma.rate * (1 - pgamma(B, 1 + gamma.shape, gamma.rate)) - 
    B * (1 - pgamma(B, gamma.shape, gamma.rate))) - 
    integrate(claimIntegral, lower = 0, upper = B)$value 
}
retRange <- uniroot(retentionRange, interval = c(-10000,10000))
retRange$root

dRange <- retRange$root

# Question 2 Part c
# For Proportional Reinsurance
expY <- gamma.shape/gamma.rate * lambda

alphaMaxFunc <- function(A, r) {
  res <- ((1+theta)*expY - (1+loading)*expY*(1-A))*r - lambda*((1-A*r/gamma.rate)^(-gamma.shape)-1)
  return(res)
}
findRootA <- function(x) {
  uniroot(alphaMaxFunc, lower = 0.0001, upper = 1, A = x)$root
}
optimR.prop <- optimise(findRootA, interval = c(0.01, 0.999), maximum = TRUE)

plot <- c()
for (i in 30:99) {
  plot <- c(plot, findRootA(i/100))
}
plot(30:99/100, plot, xlab = "alpha", ylab = "R", main = "Figure 1.1: Maximising R for alpha")

uppBound.prop <- exp(-optimR.prop$objective * c0)

# For Excess of Loss Reinsurance

integral1 <- function(x, r) {
  exp(r*x) * dgamma(x, shape = gamma.shape, rate = gamma.rate)
}

integral2 <- function(x) {
  dgamma(x, shape = gamma.shape, rate = gamma.rate)
}

dMaxFunc <- function(d, r) {
  res <- r * ((1+theta)*gamma.shape/gamma.rate - (1+loading)*(gamma.shape/gamma.rate * (1 - pgamma(d, 1 + gamma.shape, gamma.rate)) - 
                d * (1 - pgamma(d, gamma.shape, gamma.rate)))) - 
    integrate(integral1, lower = 0, upper = d, r = r)$value - exp(r*d)*integrate(integral2, lower = d, upper = Inf)$value + 1
  return(res)
}

findRootD <- function(x) {
  uniroot(dMaxFunc, lower = 0.001, upper = 100, d = x)$root
}

optimR.EoL <- optimise(findRootD, interval = c(0.001, 10), maximum = TRUE)

plot <- c()
for (i in 15:100) {
  plot <- c(plot, findRootD(i/100))
}
plot(15:100/100, plot, xlab = "d", ylab = "R", main = "Figure 1.2: Maximising R for d")

uppBound.EoL <- exp(-optimR.EoL$objective * c0)

# Question 3.
expChangeA <- 1 - (gamma.premium.rate - gamma.A.premium.rate - lambda * retention * gamma.shape/gamma.rate) / 
  (gamma.premium.rate - lambda * gamma.shape/gamma.rate)
expChangeB <- 1 - (gamma.premium.rate - gamma.B.premium.rate - lambda * integrate(claimIntegral, lower = 0, upper = retentionB)$value) /
  (gamma.premium.rate - lambda * gamma.shape/gamma.rate)

# Question 4.
# The following parameters were changed and rerun multiple times to conduct sensitivity analysis
sens.size <- 10000
sens.c0 <- 18
sens.loading <- 0.08
sens.retention <- 0.88
sens.theta <- 0.06 

sensRuinA <- function(c0, loading, retention, theta) {
  sens.premium.rate <- (1 + theta) * lambda * (gamma.shape / gamma.rate)
  sens.A.premium.rate <- (1 + loading) * lambda * (1 - retention) * (gamma.shape / gamma.rate)
  
  N <- rpois(1, lambda * t)
  
  arrival.time <- sort(runif(N, 0, t))
  
  claim.size <- rgamma(N, shape = gamma.shape, rate = gamma.rate)
  
  surplus <- c0 + (sens.premium.rate - sens.A.premium.rate) * arrival.time - retention * cumsum(claim.size)
  
  min(surplus) < 0
}

c0.values <- seq(0, 30, by = 3)
loading.values <- seq(0.05, 0.2, by = 0.05)
retention.values <- seq(0.7, 1, by = 0.05)
theta.values <- seq(0.02, 0.1, by = 0.02)

ruin.c0 <- sapply(c0.values, function(c0) mean(replicate(sens.size, sensRuinA(c0, sens.loading, sens.retention, sens.theta))))
ruin.loading <- sapply(loading.values, function(loading) mean(replicate(sens.size, sensRuinA(sens.c0, loading, sens.retention, sens.theta))))
ruin.retention <- sapply(retention.values, function(retention) mean(replicate(sens.size, sensRuinA(sens.c0, sens.loading, retention, sens.theta))))
ruin.theta <- sapply(theta.values, function(theta) mean(replicate(sens.size, sensRuinA(sens.c0, sens.loading, sens.retention, theta))))

par(mfrow = c(2, 2), oma = c(0, 0, 2, 0), mar = c(4, 4, 2, 1)) 

plot(c0.values, ruin.c0, type = "o", col = "dodgerblue", xlab = "c0", ylab = "Ruin Probability",
     main = "Initial Surplus")
abline(h = 0.01, col = "red", lty = 2)

plot(loading.values, ruin.loading, type = "o", col = "dodgerblue", xlab = "Loading", ylab = "Ruin Probability",
     main = "Reinsurer Premium")
abline(h = 0.01, col = "red", lty = 2)

plot(retention.values, ruin.retention, type = "o", col = "dodgerblue", xlab = "Retention", ylab = "Ruin Probability",
     main = "Retention Rate")
abline(h = 0.01, col = "red", lty = 2)

plot(theta.values, ruin.theta, type = "o", col = "dodgerblue", xlab = "Theta", ylab = "Ruin Probability",
     main = "Premium Loading Factor")
abline(h = 0.01, col = "red", lty = 2)
mtext("Figure 1.3: Ruin Probability Sensitivity Analysis", outer = TRUE, cex = 1.5, font = 2)
par(mfrow = c(1, 1))




