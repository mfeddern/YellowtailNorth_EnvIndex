mod<-gam
Vb <- vcov(mod)
pred2 <- predict(mod, se.fit = TRUE)
N <- 10000
BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)),V = Vb)
Cg <- predict(mod, type = "lpmatrix")
simDev <- Cg %*% t(BUdiff)
absDev <- abs(sweep(simDev, 1, pred2$se.fit, FUN = "/"))
masd <- apply(absDev, 2L, max)
crit2 <- quantile(masd, prob = 0.95, type = 8)

predData2 <- transform(cbind(data.frame(pred2), train_data),
                       uprP = fit + (crit2 * se.fit),
                       lwrP = fit - (crit2 * se.fit),
                       uprCI = fit + (2 * se.fit),
                       lwrCI = fit - (2 * se.fit))


ggplot(predData2) +
  geom_ribbon(aes(x = year, ymin = lwrP, ymax = uprP),
              alpha = 0.2, fill = "red") +
  geom_ribbon(aes(x = year, ymin = lwrCI, ymax = uprCI),
              alpha = 0.2, fill = "red") +
  geom_line(aes(x = year, y = fit)) 



### calculating prediction intervals ###
beta <- coef(gam)
V <- vcov(gam)
num_beta_vecs <- 10000
Cv <- chol(V)
set.seed(1)
nus <- rnorm(num_beta_vecs * length(beta))
beta_sims <- beta + t(Cv) %*% matrix(nus, nrow = length(beta), ncol = num_beta_vecs)
dim(beta_sims)
d_beta <- cbind(summary(gam)$se, apply(beta_sims, 1, sd))
head(d_beta)
plot(d_beta[, 1], d_beta[, 2], 
     xlab = "Calculated SE", 
     ylab = "Simulated SE")
abline(0, 1)
n_obs <- 100
sim_idx <- sample.int(nrow(train_data), size = n_obs, replace = TRUE)
sim_dat <-  train_data[sim_idx, c("CutiSTI", "DDegg", "LSTpjuv", "ONIpjuv")]
dim(sim_dat)
covar_sim <- predict(gam, newdata = sim_dat, type = "lpmatrix")
linpred_sim <- covar_sim %*% beta_sims
invlink <- function(x) x
exp_val_sim <- invlink(linpred_sim)
y_sim <- matrix(rnorm(n = prod(dim(exp_val_sim)), 
                      mean = exp_val_sim, 
                      sd = sqrt(summary(gam)$scale)), 
                nrow = nrow(exp_val_sim), 
                ncol = ncol(exp_val_sim))
dim(y_sim)
pred_int_sim <- apply(y_sim, 1, quantile, prob = c(.025, 0.975))
dim(pred_int_sim)
plot(Y_rec ~ CutiSTI, data = train_data,ylim=c(-1.5,1.5))

sim_dat_x0ord <- order(sim_dat$CutiSTI)
lines(sim_dat$CutiSTI[sim_dat_x0ord], pred_int_sim[1L, sim_dat_x0ord], col = 2, lwd = 2)
lines(sim_dat$CutiSTI[sim_dat_x0ord], pred_int_sim[2L, sim_dat_x0ord], col = 2, lwd = 2)

x0<-CutiSTI
plot(Y_rec ~ CutiSTI, data = train_data,ylim=c(-1.5,1.5))

sim_dat_x0ord <- order(sim_dat$CutiSTI)
lines(sim_dat$CutiSTI[sim_dat_x0ord], pred_int_sim[1L, sim_dat_x0ord], col = 2, lwd = 2)
lines(sim_dat$CutiSTI[sim_dat_x0ord], pred_int_sim[2L, sim_dat_x0ord], col = 2, lwd = 2)

sim_dat_pred <- predict(gam, newdata = sim_dat, se = TRUE)
lines(sim_dat$CutiSTI[sim_dat_x0ord], sim_dat_pred$fit[sim_dat_x0ord], col = 1, lwd = 2)

upr_ci <- invlink(sim_dat_pred$fit + 2*sim_dat_pred$se.fit)
lwr_ci <- invlink(sim_dat_pred$fit - 2*sim_dat_pred$se.fit)
lines(sim_dat$CutiSTI[sim_dat_x0ord], upr_ci[sim_dat_x0ord], col = 3, lwd = 2)
lines(sim_dat$CutiSTI[sim_dat_x0ord], lwr_ci[sim_dat_x0ord], col = 3, lwd = 2)
