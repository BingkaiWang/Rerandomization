####################
# Simulation design 2
# Outcome type: binary
# Estimand: Relative risk of ATE
# Randomization design: simple, rerandomization, or stratified rerandomization with pi = 0.5
# Outcome missing: no
# sample size: 400 in total
# Estimators (no missing outcomes): Unadjusted, g-comp2, ML
# Measures: bias, ese (empirical se), ase* (averge of se estimators assuming simple randomization is used), 
#           CP_N (coverage prob based on normal distribution), CP_true (coverage prob based on true asymptotic distribution)
###################
rm(list = ls())
set.seed(123)
library(tidyverse)
library(SuperLearner)
library(foreach)
library(doSNOW)
cl <- makeCluster(4)
registerDoSNOW(cl)
package_list <- c("tidyverse", "SuperLearner")
SLmethods <- c("SL.glm", "SL.rpart", "SL.nnet")
expit <- function(x){1/(1+exp(-x))}
logit <- function(x){log(x/(1-x))}

n <- 400
pi <- 0.5
n_sim <- 1000
tt <- 1

tictoc::tic()
results <- foreach(iter = 1:n_sim, .combine = cbind, .packages = package_list) %dopar% {
  
  X1 <- rnorm(n, 1, 1)
  S <- rbinom(n, size = 1, prob = 0.4 + 0.2 * (abs(X1) < 1))
  X2 <- rnorm(n, 0, 1)
  Y1 <- rbinom(n, size = 1, prob = expit(-4 + exp(X1) + 4 * X2^2 * (S+1) - abs(X2))) 
  Y0 <- rbinom(n, size = 1, prob = expit(-4 + exp(X1) - abs(X2))) 
  Xr <- cbind(X1, X2) # rerandomization variables
  
  tryCatch({
  # simple randomization----------------
  A_s <- rbinom(n, 1, pi)
  d_s <- data.frame(Y=Y1 * A_s + Y0 * (1-A_s), A=A_s, X1, X2, S)
  
  # unadjusted
  mu_1hat <- mean(d_s$Y[d_s$A==1])
  mu_0hat <- mean(d_s$Y[d_s$A==0])
  est_unadj_s <- mu_1hat/mu_0hat
  IF_unadj_s <- 1/mu_0hat * d_s$A/pi * (d_s$Y - mu_1hat) - mu_1hat/mu_0hat^2 * (1-d_s$A)/(1-pi) * (d_s$Y - mu_0hat)
  V_unadj_s <- mean(IF_unadj_s^2)
  
  # glm2
  glm2.fit <- glm(Y~A*(X1+X2+S), data = d_s, family = "binomial")
  glm2.pred1 <- predict(glm2.fit, newdata = mutate(d_s, A=1), type = "response")
  glm2.pred0 <- predict(glm2.fit, newdata = mutate(d_s, A=0), type = "response")
  mu_1hat <- mean(glm2.pred1)
  mu_0hat <- mean(glm2.pred0)
  est_glm2_s <- mu_1hat/mu_0hat
  IF_glm2_s <- 1/mu_0hat * ( d_s$A/pi * (d_s$Y - glm2.pred1) + glm2.pred1 - mu_1hat) -
    mu_1hat/mu_0hat^2 * ((1-d_s$A)/(1-pi) * (d_s$Y - glm2.pred0) + glm2.pred0 - mu_0hat)
  V_glm2_s <- mean(IF_glm2_s^2)
  
  # ML
  K <- 5
  sample_split_index <- sample(rep(1:K, each = n/K), size = n)
  sample_split <- map(1:K, ~which(sample_split_index==.))
  temp <- map(1:K, function(k){
    validation_index <- sample_split[[k]]
    training_index1 <- intersect(setdiff(1:n, validation_index), which(d_s$A==1))
    training_index0 <- intersect(setdiff(1:n, validation_index), which(d_s$A==0))
    eta_fit1 <- SuperLearner(d_s$Y[training_index1], X = d_s[training_index1,c("X1","X2","S")], family = "binomial", SL.library = SLmethods)
    eta_fit0 <- SuperLearner(d_s$Y[training_index0], X = d_s[training_index0,c("X1","X2","S")], family = "binomial", SL.library = SLmethods)
    eta.pred1 <- predict(eta_fit1, newdata = d_s[validation_index,c("X1","X2","S")])$pred
    eta.pred0 <- predict(eta_fit0, newdata = d_s[validation_index,c("X1","X2","S")])$pred
    f1 <- d_s$A[validation_index]/pi * (d_s$Y[validation_index] - eta.pred1) + eta.pred1
    f0 <- (1-d_s$A[validation_index])/(1-pi) * (d_s$Y[validation_index] - eta.pred0) + eta.pred0
    data.frame(pred1=f1, pred0=f0, fold=k)
  }) %>% Reduce(rbind,.)
  mu_1hat <- mean(temp$pred1)
  mu_0hat <- mean(temp$pred0)
  est_ml_s <- mu_1hat/mu_0hat
  V_ml_s <- map_dbl(1:K, function(k){
    var(1/mu_0hat * temp$pred1[temp$fold==k]- mu_1hat/mu_0hat^2 * temp$pred0[temp$fold==k])
  }) %>% mean
  
  # rerandomization----------------
  Q <- 100
  while (Q >= tt){
    A_star <- rbinom(n, 1, pi)
    I <- colMeans(Xr[A_star == 1,]) - colMeans(Xr[A_star == 0,])
    VI <- n/sum(A_star)/sum(1-A_star) * var(Xr)
    Q <- as.numeric(t(I) %*% solve(VI) %*% I)
  }
  A_re <- A_star
  d_re <- data.frame(Y=Y1 * A_re + Y0 * (1-A_re), A=A_re, X1, X2, S)
  
  # unadjusted
  mu_1hat <- mean(d_re$Y[d_re$A==1])
  mu_0hat <- mean(d_re$Y[d_re$A==0])
  est_unadj_re <- mu_1hat/mu_0hat
  IF_unadj_re <- 1/mu_0hat * d_re$A/pi * (d_re$Y - mu_1hat) - mu_1hat/mu_0hat^2 * (1-d_re$A)/(1-pi) * (d_re$Y - mu_0hat)
  V_unadj_re <- mean(IF_unadj_re^2)
  VI <- n/sum(A_re)/sum(1-A_re) * var(Xr)
  cIIF <- colMeans(diag((d_re$A-pi)/pi/(1-pi) * IF_unadj_re) %*% (Xr - rep(1, n) %*% t(colMeans(Xr))))
  R2_unadj_re <- t(cIIF) %*% solve(n * VI) %*% cIIF /V_unadj_re
  
  # glm2
  glm2.fit <- glm(Y~A*(X1+X2+S), data = d_re, family = "binomial")
  glm2.pred1 <- predict(glm2.fit, newdata = mutate(d_re, A=1), type = "response")
  glm2.pred0 <- predict(glm2.fit, newdata = mutate(d_re, A=0), type = "response")
  mu_1hat <- mean(glm2.pred1)
  mu_0hat <- mean(glm2.pred0)
  est_glm2_re <- mu_1hat/mu_0hat
  IF_glm2_re <- 1/mu_0hat * ( d_re$A/pi * (d_re$Y - glm2.pred1) + glm2.pred1 - mu_1hat) -
    mu_1hat/mu_0hat^2 * ((1-d_re$A)/(1-pi) * (d_re$Y - glm2.pred0) + glm2.pred0 - mu_0hat)
  V_glm2_re <- mean(IF_glm2_re^2)
  VI <- n/sum(A_re)/sum(1-A_re) * var(Xr)
  cIIF <- colMeans(diag((d_re$A-pi)/pi/(1-pi) * IF_glm2_re) %*% (Xr - rep(1, n) %*% t(colMeans(Xr))))
  R2_glm2_re <- t(cIIF) %*% solve(n * VI) %*% cIIF /V_glm2_re
  
  # ML
  K <- 5
  sample_split_index <- sample(rep(1:K, each = n/K), size = n)
  sample_split <- map(1:K, ~which(sample_split_index==.))
  temp <- map(1:K, function(k){
    validation_index <- sample_split[[k]]
    training_index1 <- intersect(setdiff(1:n, validation_index), which(d_s$A==1))
    training_index0 <- intersect(setdiff(1:n, validation_index), which(d_s$A==0))
    eta_fit1 <- SuperLearner(d_re$Y[training_index1], X = d_re[training_index1,c("X1","X2","S")], family = "binomial", SL.library = SLmethods)
    eta_fit0 <- SuperLearner(d_re$Y[training_index0], X = d_re[training_index0,c("X1","X2","S")], family = "binomial", SL.library = SLmethods)
    eta.pred1 <- predict(eta_fit1, newdata = d_re[validation_index,c("X1","X2","S")])$pred
    eta.pred0 <- predict(eta_fit0, newdata = d_re[validation_index,c("X1","X2","S")])$pred
    f1 <- d_re$A[validation_index]/pi * (d_re$Y[validation_index] - eta.pred1) + eta.pred1
    f0 <- (1-d_re$A[validation_index])/(1-pi) * (d_re$Y[validation_index] - eta.pred0) + eta.pred0
    mutate(d_re[validation_index,], pred1=f1, pred0=f0, fold=k)
  }) %>% Reduce(rbind,.)
  mu_1hat <- mean(temp$pred1)
  mu_0hat <- mean(temp$pred0)
  est_ml_re <- mu_1hat/mu_0hat
  V_ml_re <- map_dbl(1:K, function(k){
    var(1/mu_0hat * temp$pred1[temp$fold==k]- mu_1hat/mu_0hat^2 * temp$pred0[temp$fold==k])
  }) %>% mean
  IF <- as.vector(1/mu_0hat *temp$pred1 - mu_1hat/mu_0hat^2 * temp$pred0-est_ml_re)
  cIIF <- colMeans(diag((temp$A-pi)/pi/(1-pi) * IF) %*% as.matrix(temp[,c("X1","X2")] - rep(1, n) %*% t(colMeans(temp[,c("X1","X2")]))))
  VI <- n/sum(A_re)/sum(1-A_re) * var(Xr)
  R2_ml_re <- t(cIIF) %*% solve(n * VI) %*% cIIF /V_ml_re
  
  # stratified rerandomization----------------
  Q <- 100
  while(Q >= tt){
    A_tilde <- rep(0, n)
    for(s in unique(S)){
      A_tilde[sample(which(S==s), size = pi * sum(S==s))] <- 1
    }
    I <- colMeans(Xr[A_tilde == 1,]) - colMeans(Xr[A_tilde == 0,])
    VI <- n/sum(A_tilde)/sum(1-A_tilde) * map(unique(S), function(s){var(Xr[S==s,]) * mean(S==s)}) %>% Reduce(`+`,.)
    Q <- as.numeric(t(I) %*% solve(VI) %*% I)
  }
  A_stre <- A_tilde
  d_stre <- data.frame(Y=Y1 * A_stre + Y0 * (1-A_stre), A=A_stre, X1, X2, S)
  
  # unadjusted
  mu_1hat <- mean(d_stre$Y[d_stre$A==1])
  mu_0hat <- mean(d_stre$Y[d_stre$A==0])
  est_unadj_stre <- mu_1hat/mu_0hat
  IF_unadj_stre <- 1/mu_0hat * d_stre$A/pi * (d_stre$Y - mu_1hat) - mu_1hat/mu_0hat^2 * (1-d_stre$A)/(1-pi) * (d_stre$Y - mu_0hat)
  V_unadj_stre <- mean(IF_unadj_stre^2)
  VI <- n/sum(A_stre)/sum(1-A_stre) * var(Xr)
  cIIF <- colMeans(diag((d_stre$A-pi)/pi/(1-pi) * IF_unadj_stre) %*% (Xr - rep(1, n) %*% t(colMeans(Xr))))
  R2_unadj_stre <- t(cIIF) %*% solve(n * VI) %*% cIIF /V_unadj_stre
  
  # glm2
  glm2.fit <- glm(Y~A*(X1+X2+S), data = d_stre, family = "binomial")
  glm2.pred1 <- predict(glm2.fit, newdata = mutate(d_stre, A=1), type = "response")
  glm2.pred0 <- predict(glm2.fit, newdata = mutate(d_stre, A=0), type = "response")
  mu_1hat <- mean(glm2.pred1)
  mu_0hat <- mean(glm2.pred0)
  est_glm2_stre <- mu_1hat/mu_0hat
  IF_glm2_stre <- 1/mu_0hat * ( d_stre$A/pi * (d_stre$Y - glm2.pred1) + glm2.pred1 - mu_1hat) -
    mu_1hat/mu_0hat^2 * ((1-d_stre$A)/(1-pi) * (d_stre$Y - glm2.pred0) + glm2.pred0 - mu_0hat)
  V_glm2_stre <- mean(IF_glm2_stre^2)
  VI <- n/sum(A_stre)/sum(1-A_stre) * var(Xr)
  cIIF <- colMeans(diag((d_stre$A-pi)/pi/(1-pi) * IF_glm2_stre) %*% (Xr - rep(1, n) %*% t(colMeans(Xr))))
  R2_glm2_stre <- t(cIIF) %*% solve(n * VI) %*% cIIF /V_glm2_stre
  
  # ML
  K <- 5
  temp <- map(unique(d_stre$S), function(s){
    d_stre_s <- d_stre[d_stre$S==s,]
    sample_split_index <- rep(0, nrow(d_stre_s))
    for(a in 0:1){
      nn <- sum(d_stre_s$A == a)
      sample_split_index[d_stre_s$A == a] <- sample(rep(1:K, each = ceiling(nn/K)), size = nn)
    }
    sample_split <- map(1:K, ~which(sample_split_index==.))
    temp1 <- map(1:K, function(k){
      validation_index <- sample_split[[k]]
      training_index1 <- intersect(setdiff(1:nrow(d_stre_s), validation_index), which(d_stre_s$A==1))
      training_index0 <- intersect(setdiff(1:nrow(d_stre_s), validation_index), which(d_stre_s$A==0))
      eta_fit1 <- SuperLearner(d_stre_s$Y[training_index1], X = d_stre_s[training_index1,c("X1","X2")], family = "binomial", SL.library = SLmethods)
      eta_fit0 <- SuperLearner(d_stre_s$Y[training_index0], X = d_stre_s[training_index0,c("X1","X2")], family = "binomial", SL.library = SLmethods)
      eta.pred1 <- predict(eta_fit1, newdata = d_stre_s[validation_index,c("X1","X2")])$pred
      eta.pred0 <- predict(eta_fit0, newdata = d_stre_s[validation_index,c("X1","X2")])$pred
      f1 <- d_stre_s$A[validation_index]/pi * (d_stre_s$Y[validation_index] - eta.pred1) + eta.pred1
      f0 <- (1-d_stre_s$A[validation_index])/(1-pi) * (d_stre_s$Y[validation_index] - eta.pred0) + eta.pred0
      mutate(d_stre_s[validation_index,], pred1=f1, pred0=f0, fold=k)
    }) %>% Reduce(rbind,.)
  }) %>% Reduce(rbind,.)
  mu_1hat <- mean(temp$pred1)
  mu_0hat <- mean(temp$pred0)
  est_ml_stre <- mu_1hat/mu_0hat
  V_ml_stre <- map_dbl(1:K, function(k){
    map_dbl(unique(d_stre$S), function(s){
      var(1/mu_0hat * temp$pred1[temp$fold==k & temp$S==s]- mu_1hat/mu_0hat^2 * temp$pred0[temp$fold==k & temp$S==s]) * mean(d_stre$S==s)
    }) %>% sum
  }) %>% mean
  IF <- as.vector(1/mu_0hat *temp$pred1 - mu_1hat/mu_0hat^2 * temp$pred0-est_ml_re)
  cIIF <- colMeans(diag((temp$A-pi)/pi/(1-pi) * IF) %*% as.matrix(temp[,c("X1","X2")] - rep(1, n) %*% t(colMeans(temp[,c("X1","X2")]))))
  VI <- n/sum(A_re)/sum(1-A_re) * var(Xr)
  R2_ml_stre <- t(cIIF) %*% solve(n * VI) %*% cIIF /V_ml_stre
  
  # result
  c(est_unadj_s, V_unadj_s, 0, est_glm2_s, V_glm2_s, 0, est_ml_s, V_ml_s, 0,
    est_unadj_re, V_unadj_re, R2_unadj_re, est_glm2_re, V_glm2_re, R2_glm2_re, est_ml_re, V_ml_re, R2_ml_re,
    est_unadj_stre, V_unadj_stre, R2_unadj_stre, est_glm2_stre, V_glm2_stre, R2_glm2_stre, est_ml_stre, V_ml_stre, R2_ml_stre  )
  }, error = function(e) { return(rep(NA, 27))})
}
  
tictoc::toc()
stopCluster(cl)

saveRDS(results, file = "simulation/sim2-1.rds")

# summarizing results
rm(list = ls())
results <- readRDS(file = "simulation/sim2-1.rds")
n_sim <- ncol(results)
tt <- 1
n <- 400
expit <- function(x){1/(1+exp(-x))}
nn <- 1000000
pi <- 0.5
X1 <- rnorm(nn, 1, 1)
S <- rbinom(nn, size = 1, prob = 0.4 + 0.2 * (abs(X1) < 1))
X2 <- rnorm(nn, 0, 1)
Y1 <- rbinom(nn, size = 1, prob = expit(-4 + exp(X1) + 4 * X2^2 * (S+1) - abs(X2))) 
Y0 <- rbinom(nn, size = 1, prob = expit(-4 + exp(X1) - abs(X2))) 
delta <- mean(Y1)/mean(Y0)

est_index <- seq(1, 27, by = 3)
ase_index <- est_index + 1
r2_index <- est_index + 2
z <- rnorm(10000)
d1 <- rnorm(10000)
d2 <- rnorm(10000)
index <- which(d1^2 + d2^2 < tt)
z <- z[index]
d1 <- d1[index]


cp_true <- map_dbl(1:6, function(j){
  map_dbl(1: n_sim, function(i){
    emp <- (sqrt(1-results[r2_index[j], i]) * z + sqrt(results[r2_index[j], i]) * d1) * sqrt(results[ase_index[j],i]/n) + delta
    (results[est_index[j], i] >= quantile(emp, 0.025, na.rm=T)) & (results[est_index[j], i] <= quantile(emp, 0.975, na.rm=T))
  }) %>% mean(na.rm = T)
})



summary_results <- cbind(bias = apply(results[est_index,], 1, mean, na.rm = T) - delta,
                         ese = apply(results[est_index,], 1, sd, na.rm = T), 
                         ase = apply(results[ase_index,], 1, function(x){mean(sqrt(x/n),na.rm=T)}),
                         cp_normal = apply(abs(results[est_index,] - delta)/sqrt(results[ase_index,]/n) <= qnorm(0.975), 1, mean, na.rm = T),
                         cp_true = cp_true
)


rownames(summary_results) <- rep(c("Unadj", "GLM2", "ML"), 3)

round(summary_results, 2)
xtable::xtable(summary_results)

