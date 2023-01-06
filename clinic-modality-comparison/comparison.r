library(TMB)
library(RColorBrewer)
## Sorting out data.
dwm.df <- as.matrix(read.csv("data/mo-dwm.csv", header = FALSE))[, -1]
ef.df <- as.matrix(read.csv("data/mo-ef.csv", header = FALSE))[, -1]
egs.df <- as.matrix(read.csv("data/mo-egs.csv", header = FALSE))[, -1]
emb.df <- as.matrix(read.csv("data/mo-emb.csv", header = FALSE))[, -1]
gs.df <- as.matrix(read.csv("data/mo-gs.csv", header = FALSE))[, -1]
times <- c(0, 3, 4, 7, 9, 10, 34)
all.times <- seq(0, 34, by = 0.5)
n.all.times <- length(all.times)
time.indices <- which(all.times %in% times) - 1
end.treatment <- c(20, 10)
n.birds <- nrow(dwm.df)
n.times <- length(times)
n.measurements <- c(rep(7, 6), rep(6, 2), rep(7, 8))
dwm.df[is.na(dwm.df)] <- 9999
ef.df[is.na(ef.df)] <- 9999
egs.df[is.na(egs.df)] <- 9999
emb.df[is.na(emb.df)] <- 9999
gs.df[is.na(gs.df)] <- 9999

## Plotting the raw data.
dwm.na.df <- dwm.df
dwm.na.df[dwm.na.df == 9999] <- NA
ef.na.df <- ef.df
ef.na.df[ef.na.df == 9999] <- NA
egs.na.df <- egs.df
egs.na.df[egs.na.df == 9999] <- NA
emb.na.df <- emb.df
emb.na.df[emb.na.df == 9999] <- NA
gs.na.df <- gs.df
gs.na.df[gs.na.df == 9999] <- NA

pdf(width = 8, height = 10, file = "data-plot.pdf")
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1), las = 1)
dfs <- c("dwm.na.df", "gs.na.df", "egs.na.df", "emb.na.df", "ef.na.df")
method.names <- c("DWM", "GS", "MSGS", "MSMB", "MST")
for (i in 1:5){
    use.df <- get(dfs[i])
    plot.new()
    plot.window(xlim = range(times), ylim = range(use.df, na.rm = TRUE))
    box()
    axis(1)
    axis(2)
    for (j in 1:n.birds){
        lines(times, use.df[j, ])
    }
    title(main = method.names[i])
    if (i %in% c(4, 5)){
        title(xlab = "Days since first treatment")
    }
    if (i %in% c(1, 3, 5)){
        title(ylab = "MO organisms detected")
    }
}
dev.off()


## Figuring out which animals were still infected at the end of treatment.
check <- 6
all.df <- cbind(dwm.df[, check], ef.df[, check], egs.df[, check], emb.df[, check], gs.df[, check])
logical.infected.birds.et <- apply(all.df, 1, function(x) any(x > 0 & x < 9999))
infected.birds.et <- which(logical.infected.birds.et)
n.infected.birds.et <- length(infected.birds.et)

## Figuring out which animals were still infected at follow-up.
check <- 7
all.df <- cbind(dwm.df[, check], ef.df[, check], egs.df[, check], emb.df[, check], gs.df[, check])
infected.birds.fu <- which(apply(all.df, 1, function(x) any(x > 0 & x < 9999)))
n.infected.birds.fu <- length(infected.birds.fu)

## Loading estimated sickness levels from infected birds.
load("s-infected.RData")
test.cases <- s.infected.vec
n.test.cases <- length(test.cases)

compile("comparison.cpp")
dyn.load(dynlib("comparison"))

mo.data <- list(n_birds = n.birds,
                times = times,
                all_times = all.times,
                time_indices = time.indices,
                end_treatment = end.treatment,
                cov_function = 1, ## 1 for exponential, 2 for squared exponential.
                time_relationship = 3, ## 1 for linear, 2 for quadratic, 3 for piecewise linear.
                n_times = n.times,
                n_measurements = n.measurements,
                max_measurements = max(n.measurements),
                dwm = dwm.df,
                ef = ef.df,
                egs = egs.df,
                emb = emb.df,
                gs = gs.df,
                infected_birds = as.numeric(logical.infected.birds.et),
                test_cases = test.cases)
mo.parameters <-  list(beta_0_base = log(50),
                       beta_0_diff = rep(0, 4),
                       beta_1_base = 1,
                       beta_1_diff = rep(0, 4),
                       alpha_1 = 0,
                       alpha_2 = 0,
                       log_theta_base = log(5),
                       log_theta_diff = rep(0, 4),
                       log_sigma_s = 0.1,
                       log_rho_s = log(1),
                       log_sigma_u1 = log(1),
                       log_sigma_u2 = log(1),
                       s = matrix(1, nrow = n.birds, ncol = n.all.times),
                       u_1 = rep(1, n.birds), u_2 = rep(1, n.birds))

## Object for full model.
mo.obj <- MakeADFun(data = mo.data,
                    parameters = mo.parameters,
                    map = list(##beta_0_diff = factor(rep(NA, 4)),
                               beta_1_base = factor(NA),
                               beta_1_diff = factor(rep(NA, 4))
                               ##alpha_1 = factor(NA),
                               ##alpha_2 = factor(NA)
                               ##log_theta_diff = factor(rep(NA, 4))
                               ##log_rho_s = factor(NA)
                               ),
                    random = c("s", "u_1", "u_2"), DLL = "comparison")

## Object for model with alpha_1 = 0.
mo.obj.a1.null <- MakeADFun(data = mo.data,
                            parameters = mo.parameters,
                            map = list(##beta_0_diff = factor(rep(NA, 4)),
                                       beta_1_base = factor(NA),
                                       beta_1_diff = factor(rep(NA, 4)),
                                       alpha_1 = factor(NA)
                                       ##alpha_2 = factor(NA)
                                       ##log_theta_diff = factor(rep(NA, 4))
                                       ##log_rho_s = factor(NA)
                                       ),
                            random = c("s", "u_1", "u_2"), DLL = "comparison")

## Object for model with alpha_2 = 0.
mo.obj.a2.null <- MakeADFun(data = mo.data,
                            parameters = mo.parameters,
                            map = list(##beta_0_diff = factor(rep(NA, 4)),
                                       beta_1_base = factor(NA),
                                       beta_1_diff = factor(rep(NA, 4)),
                                       alpha_2 = factor(NA)
                                       ##alpha_2 = factor(NA)
                                       ##log_theta_diff = factor(rep(NA, 4))
                                       ##log_rho_s = factor(NA)
                                       ),
                            random = c("s", "u_1", "u_2"), DLL = "comparison")

## Object for model with alpha_1 = alpha_2.
mo.data.a1.a2 <- mo.data
mo.data.a1.a2$time_relationship <- 1
mo.obj.a1.a2.null <- MakeADFun(data = mo.data.a1.a2,
                               parameters = mo.parameters,
                               map = list(##beta_0_diff = factor(rep(NA, 4)),
                                   beta_1_base = factor(NA),
                                   beta_1_diff = factor(rep(NA, 4)),
                                   alpha_2 = factor(NA)
                                   ##log_theta_diff = factor(rep(NA, 4))
                                   ##log_rho_s = factor(NA)
                               ),
                               random = c("s", "u_1", "u_2"), DLL = "comparison")

## Object for model with the same beta0 for all treatments.
mo.obj.b0.null <- MakeADFun(data = mo.data,
                            parameters = mo.parameters,
                            map = list(beta_0_diff = factor(rep(NA, 4)),
                                       beta_1_base = factor(NA),
                                       beta_1_diff = factor(rep(NA, 4))
                                       ##alpha_1 = factor(NA),
                                       ##alpha_2 = factor(NA)
                                       ##log_theta_diff = factor(rep(NA, 4))
                                       ##log_rho_s = factor(NA)
                                       ),
                            random = c("s", "u_1", "u_2"), DLL = "comparison")

mo.fit <- nlminb(mo.obj$par, mo.obj$fn, mo.obj$gr)
mo.fit.a1.null <- nlminb(mo.obj.a1.null$par, mo.obj.a1.null$fn, mo.obj.a1.null$gr)
mo.fit.a2.null <- nlminb(mo.obj.a2.null$par, mo.obj.a2.null$fn, mo.obj.a2.null$gr)
mo.fit.a1.a2.null <- nlminb(mo.obj.a1.a2.null$par, mo.obj.a1.a2.null$fn, mo.obj.a1.a2.null$gr)
mo.fit.b0.null <- nlminb(mo.obj.b0.null$par, mo.obj.b0.null$fn, mo.obj.b0.null$gr)

## AIC (remember objective is already the NLL). 1740.224 is current best model.
2*mo.fit$objective + 2*length(mo.fit$par)

## Choosing which model to inspect.
mo.obj <- mo.obj
## Creating report.
mo.rep <- sdreport(mo.obj)
## Looking at fixed parameters.
fixed.summary <- summary(mo.rep, "fixed")
fixed.summary
## Some p-values.
cbind(summary(mo.rep, "fixed"),
      1 - pnorm(abs(summary(mo.rep, "fixed")[, 1]/summary(mo.rep, "fixed")[, 2])))
## Looking at reported parameters.
summary(mo.rep, "report")
## Looking at random effects.
summary(mo.rep, "random")
## Collecting random effects.
s <- matrix(summary(mo.rep, "random")[rownames(summary(mo.rep, "random")) == "s", ][, 1], nrow = 16)
## Collecting mean of random effects.
mu.s <- matrix(summary(mo.rep, "report")[rownames(summary(mo.rep, "report")) == "mu_s", 1], nrow = n.birds)
## Collecting lambda.y for dwm.
lambda.y.rep <- matrix(summary(mo.rep, "report")[rownames(summary(mo.rep, "report")) == "lambda_y_mat", 1],
                       nrow = n.birds, ncol = max(n.measurements))
## Collecting beta0 comparisons.
beta.0.comparisons.ests <- summary(mo.rep, "report")[rownames(summary(mo.rep, "report")) ==
                                                "beta_0_comparison", 1]
beta.0.comparisons.se <- summary(mo.rep, "report")[rownames(summary(mo.rep, "report")) ==
                                                   "beta_0_comparison", 2]
beta.0.comparisons.z <- beta.0.comparisons.ests/beta.0.comparisons.se
beta.0.comparisons.p <- 2*pnorm(-abs(beta.0.comparisons.z))
beta.0.comparisons.ests.mat <-  matrix(beta.0.comparisons.ests,
                                       nrow = 5, ncol = 5)
beta.0.comparisons.z.mat <-  matrix(beta.0.comparisons.z,
                                    nrow = 5, ncol = 5)
beta.0.comparisons.p.mat <-  matrix(beta.0.comparisons.p,
                                    nrow = 5, ncol = 5)
## Collecting p_zero comparisons.
p.zero.comparisons.ests <- summary(mo.rep, "report")[rownames(summary(mo.rep, "report")) ==
                                                     "p_zero_comparison", 1]
p.zero.comparisons.se <- summary(mo.rep, "report")[rownames(summary(mo.rep, "report")) ==
                                                   "p_zero_comparison", 2]
p.zero.comparisons.z <- p.zero.comparisons.ests/p.zero.comparisons.se
p.zero.comparisons.p <- 2*pnorm(-abs(p.zero.comparisons.z))
p.zero.comparisons.ests.mat <- vector(mode = "list", length = n.test.cases)
p.zero.comparisons.se.mat <- vector(mode = "list", length = n.test.cases)
p.zero.comparisons.z.mat <- vector(mode = "list", length = n.test.cases)
p.zero.comparisons.p.mat <- vector(mode = "list", length = n.test.cases)
for (i in 1:n.test.cases){
    indices <- ((i - 1)*25 + 1):((i - 1)*25 + 25)
    p.zero.comparisons.ests.mat[[i]] <- matrix(p.zero.comparisons.ests[indices], nrow = 5, ncol = 5)
    p.zero.comparisons.se.mat[[i]] <- matrix(p.zero.comparisons.se[indices], nrow = 5, ncol = 5)
    p.zero.comparisons.z.mat[[i]] <- matrix(p.zero.comparisons.z[indices], nrow = 5, ncol = 5)
    p.zero.comparisons.p.mat[[i]] <- matrix(p.zero.comparisons.p[indices], nrow = 5, ncol = 5)
}
## P-values seem suspiciously small to me, in light of the CIs plotted
## below. I think it might be the skewness of the logits' sampling
## distributions messing the delta method up.
sig.ps <- lapply(p.zero.comparisons.p.mat, function(x) which(x < 0.05, arr.ind = TRUE))
lapply(sig.ps, function(x) sum(x[, 2] == 4))

## Plotting mu.s for each bird.
plot(all.times, mu.s[1, ], type = "l", ylim = range(mu.s))
for (i in 1:n.birds){
    lines(all.times, mu.s[i, ])
}


## Plotting covariance function.
rho.s <- summary(mo.rep, "report")[rownames(summary(mo.rep, "report")) == "rho_s", 1]
sigma.s <- summary(mo.rep, "report")[rownames(summary(mo.rep, "report")) == "sigma_s", 1]
if (mo.data$cov_function == 1){
    cov <- sigma.s^2*exp(-seq(0, 35, by = 0.1)/sigma.s)
} else if (mo.data$cov_function == 2){
    cov <- sigma.s^2*exp(-(seq(0, 35, by = 0.1)^2)/sigma.s^2) 
}
plot(seq(0, 35, by = 0.1), cov, ylim = c(0, max(cov)), type = "l")

## Transforming bird sickness latent variables to percentage change.
s.perc <- 100*(exp(s) - 1)

## Plot of bird sicknesses.
pdf(height = 4, width = 6, file = "bird-sickness.pdf")
par(mar = c(4, 4, 0, 0), oma = rep(0.1, 4), xaxs = "i", las = 1)
cols <- brewer.pal(3, "Set1")[c(2, 3)]
plot(all.times, s.perc[1, ],  ylim = c(-100, max(s.perc)), type = "n",
     ylab = "Percentage difference in organisms shed",
     xlab = "Days since first treatment", axes = FALSE)
axis(1)
axis(2, at = seq(-100, 800, by = 100))
abline(h = -100, col = "grey")
box()
for (i in 1:n.birds){
    lines(all.times, s.perc[i, ], col = cols[logical.infected.birds.et[i] + 1])
}
abline(v = c(10, 20), lty = "dotted")
legend("topright", c("Ten-day treatment", "Twenty-day treatment"), lty = rep(1, 2),
       col = cols, bty = "n", cex = 0.8)
dev.off()

## Estimating probabilities of nonzero counts for infected birds
## post-treatment and at follow-up.

## Collecting estimated sickness values for infected birds for the
## last two tests.
s.infected <- s[infected.birds.fu, time.indices[6:7] + 1]
s.infected.vec <- c(s.infected)
save(s.infected.vec, file = "s-infected.RData")

## Collecting p.zero for end of treatment.
logit.p.zero <- matrix(summary(mo.rep, "report")[rownames(summary(mo.rep, "report")) == "logit_p_zero", 1],
                       nrow = length(mo.data$test_cases))
logit.p.zero.se <- matrix(summary(mo.rep, "report")[rownames(summary(mo.rep, "report")) == "logit_p_zero", 2],
                             nrow = length(mo.data$test_cases))
logit.p.zero.lower <- logit.p.zero - qnorm(0.975)*logit.p.zero.se
logit.p.zero.upper <- logit.p.zero + qnorm(0.975)*logit.p.zero.se

## Old and new method names.
methods <- c("DWM", "EF", "EGS", "EMB", "GS")
methods.new <- c("DWM", "MST", "MSGS", "MSMB", "GS")
methods.order <- rev(order(methods.new))
logit.p.zero <- logit.p.zero[, methods.order]
logit.p.zero.se <- logit.p.zero.se[, methods.order]
logit.p.zero.lower <- logit.p.zero.lower[, methods.order]
logit.p.zero.upper <- logit.p.zero.upper[, methods.order]
methods.new <- methods.new[methods.order]

## Transforming to probability scale.
p.zero <- plogis(logit.p.zero)
p.zero.lower <- plogis(logit.p.zero.lower)
p.zero.upper <- plogis(logit.p.zero.upper)

n.methods <- length(methods.new)
colnames(p.zero) <- colnames(p.zero.lower) <- colnames(p.zero.upper) <-
  methods.new
rownames(p.zero) <- rownames(p.zero.lower) <- rownames(p.zero.upper) <-
    rep(infected.birds.fu, 2)
## Probability of detecting infection.
p.nonzero <- 1 - p.zero
## Lower- and upper-limits for confidence intervals.
p.nonzero.lower <- 1 - p.zero.upper
p.nonzero.upper <- 1 - p.zero.lower

## Plotting the confidence intervals.
pdf(width = 8, height = 12, file = "ci-plot.pdf")
par(mar = c(4, 4, 0, 0), oma = rep(0.1, 4), yaxs = "i")
plot.new()
bird.spacing <- 20
time.spacing <- 10
method.spacing <- 1
bird.y.pos <- bird.spacing*(1:n.infected.birds.fu) - bird.spacing/2
plot.window(xlim = c(0, 1), ylim = c(0, max(bird.y.pos) + bird.spacing/2))
cols <- brewer.pal(10, "Paired")
cols.et <- cols[c(TRUE, FALSE)]
cols.fu <- cols[c(FALSE, TRUE)]
abline(v = 1, lty = "dotted")
lines(rep(0, 2), y = c(0, bird.y.pos[n.infected.birds.fu]), lty = "dotted")
for (i in 1:n.infected.birds.fu){
    for (j in 1:2){
        for (k in 1:n.methods){
            y.val <- bird.y.pos[i] + c(-1, 1)[j]*time.spacing/2 + (k - 3)*method.spacing 
            p <- p.nonzero[rownames(p.nonzero) == infected.birds.fu[i], k][j]
            p.lower <- p.nonzero.lower[rownames(p.nonzero) == infected.birds.fu[i], k][j]
            p.upper <- p.nonzero.upper[rownames(p.nonzero) == infected.birds.fu[i], k][j]
            col <- c(cols.et[k], cols.fu[k])[j]
            points(p, y.val, pch = 16, col = col)
            lines(c(p.lower, p.upper), rep(y.val, 2), col = col)
        }
    }
}
abline(h = bird.y.pos[-n.infected.birds.fu] + bird.spacing/2)
box()
axis(1)
axis(2, bird.y.pos, paste("Bird", infected.birds.fu, "\n", "\n"), tick = FALSE)
axis(2, rep(bird.y.pos, each = 2) + rep(c(-1, 1)*time.spacing/2, n.infected.birds.fu),
     rep(c("ET", "FU"), n.infected.birds.fu), tick = FALSE)
title(xlab = "Probability of detecting an MO organism")
legend("topleft", rev(methods.new), lty = rep(1, n.methods), pch = 16, col = rev(cols.fu), bty = "n",
       bg = "white", cex = 0.8)
dev.off()

## Likelihood-ratio test for treatment effect.
lrts.a1 <- 2*(mo.fit.a1.null$objective - mo.fit$objective)
lrts.a1.df <- length(mo.fit$par) - length(mo.fit.a1.null$par)
1 - pchisq(lrts.a1, lrts.a1.df)

## Likelihood-ratio test for after-treatment effect.
lrts.a2 <- 2*(mo.fit.a2.null$objective - mo.fit$objective)
lrts.a2.df <- length(mo.fit$par) - length(mo.fit.a2.null$par)
1 - pchisq(lrts.a2, lrts.a2.df)

## Likelihood-ratio test for after-treatment effect being the same as
## treatment effect.
lrts.a1.a2 <- 2*(mo.fit.a1.a2.null$objective - mo.fit$objective)
lrts.a1.a2.df <- length(mo.fit$par) - length(mo.fit.a1.a2.null$par)
1 - pchisq(lrts.a1.a2, lrts.a1.a2.df)

## Likelihood-ratio test for differences between treatments.
lrts.b0 <- 2*(mo.fit.b0.null$objective - mo.fit$objective)
lrts.b0.df <- length(mo.fit$par) - length(mo.fit.b0.null$par)
1 - pchisq(lrts.b0, lrts.b0.df)

save.image("out.RData")
