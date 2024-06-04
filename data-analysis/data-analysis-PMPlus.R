rm(list = ls())
set.seed(123)
library(tidyverse)
library(lme4)
library(SuperLearner)
SLmethods <- c("SL.glm", "SL.rpart", "SL.nnet")
d <- read.csv("data-analysis/PMPlus_Data.csv")

# preprocessing
d <- d %>% filter(visit == 4) %>%
  select(cluster, A=INT, Y=GHQ_total, gender_cluster, mh_access, urban, disaster_risk, BL_GHQ_total)
d$cluster <- as.character(d$cluster)
d$gender_cluster <- d$gender_cluster == 1
d$mh_access <- d$mh_access == "High access"
d$urban <- d$urban == "Urban"
d$disaster_risk <- d$disaster_risk == "High"
cluster_d <- d %>% group_by(cluster) %>%summarise_all(mean) %>% 
  mutate(n=as.data.frame(group_by(d, cluster) %>% summarise(n = n()))[,2])
pi <- mean(cluster_d$A)
n <- nrow(cluster_d)

# unadjusted
mu_1hat <- mean(cluster_d$Y[cluster_d$A==1])
mu_0hat <- mean(cluster_d$Y[cluster_d$A==0])
est_unadj <- mu_1hat - mu_0hat
IF_unadj <- cluster_d$A/pi * (cluster_d$Y - mu_1hat) - (1-cluster_d$A)/(1-pi) * (cluster_d$Y - mu_0hat)
V_unadj <- mean(IF_unadj^2)/n

# cluster-level ANCOVA
ancova.fit <- lm(Y~.-cluster, data = cluster_d)
ancova.pred1 <- predict(ancova.fit, newdata = mutate(cluster_d, A=1), type = "response")
ancova.pred0 <- predict(ancova.fit, newdata = mutate(cluster_d, A=0), type = "response")
mu_1hat <- mean(ancova.pred1)
mu_0hat <- mean(ancova.pred0)
est_ancova <- mu_1hat - mu_0hat
IF_ancova <- ( cluster_d$A/pi * (cluster_d$Y - ancova.pred1) + ancova.pred1 - mu_1hat) -
  ((1-cluster_d$A)/(1-pi) * (cluster_d$Y - ancova.pred0) + ancova.pred0 - mu_0hat)
V_ancova <- mean(IF_ancova^2)/n

# mixed-ANCOVA
mixed.fit <- lmer(Y~.-cluster + (1|cluster), data=d)
est_mixed <- summary(mixed.fit)$coefficients['A', 'Estimate']
V_mixed <- summary(mixed.fit)$coefficients['A', 'Std. Error']^2

# ML
K <- 3
covariates_names <- c("gender_cluster", "mh_access", "urban", "disaster_risk", "BL_GHQ_total")
sample_split_index <- sample(rep(1:K, each = n/K), size = n)
sample_split <- map(1:K, ~which(sample_split_index==.))
temp <- map(1:K, function(k){
  validation_index <- sample_split[[k]]
  training_index1 <- intersect(setdiff(1:n, validation_index), which(cluster_d$A==1))
  training_index0 <- intersect(setdiff(1:n, validation_index), which(cluster_d$A==0))
  eta_fit1 <- SuperLearner(cluster_d$Y[training_index1], X = cluster_d[training_index1,covariates_names], family = "gaussian", SL.library = SLmethods)
  eta_fit0 <- SuperLearner(cluster_d$Y[training_index0], X = cluster_d[training_index0,covariates_names], family = "gaussian", SL.library = SLmethods)
  eta.pred1 <- predict(eta_fit1, newdata = cluster_d[validation_index,covariates_names])$pred
  eta.pred0 <- predict(eta_fit0, newdata = cluster_d[validation_index,covariates_names])$pred
  f1 <- cluster_d$A[validation_index]/pi * (cluster_d$Y[validation_index] - eta.pred1) + eta.pred1
  f0 <- (1-cluster_d$A[validation_index])/(1-pi) * (cluster_d$Y[validation_index] - eta.pred0) + eta.pred0
  data.frame(pred1=f1, pred0=f0, fold=k)
}) %>% Reduce(rbind,.)
est_ml <- mean(temp$pred1 - temp$pred0)
V_ml <- map_dbl(1:K, function(k){
  var(temp$pred1[temp$fold==k]-temp$pred0[temp$fold==k])/n
}) %>% mean

summary_results <- data.frame(est = c(est_unadj, est_ancova, est_mixed, est_ml),
                              se = sqrt(c(V_unadj, V_ancova, V_mixed, V_ml))) %>%
  mutate(cp.lower = est - 1.96 * se, cp.upper = est + 1.96 * se)
rownames(summary_results) <- c("unadj", "ancova", "mixed", "ML")
xtable::xtable(summary_results)
 