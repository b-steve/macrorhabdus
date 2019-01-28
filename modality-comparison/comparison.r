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
n.chickens <- nrow(dwm.df)
n.times <- length(times)
n.measurements <- c(rep(7, 6), rep(6, 2), rep(7, 8))
dwm.df[is.na(dwm.df)] <- 9999
ef.df[is.na(ef.df)] <- 9999
egs.df[is.na(egs.df)] <- 9999
emb.df[is.na(emb.df)] <- 9999
gs.df[is.na(gs.df)] <- 9999

## Figuring out which animals were still infected at the end of treatment.
check <- 6
all.df <- cbind(dwm.df[, check], ef.df[, check], egs.df[, check], emb.df[, check], gs.df[, check])
logical.infected.chickens.et <- apply(all.df, 1, function(x) any(x > 0 & x < 9999))
infected.chickens.et <- which(logical.infected.chickens.et)
n.infected.chickens.et <- length(infected.chickens.et)

## Figuring out which animals were still infected at follow-up.
check <- 7
all.df <- cbind(dwm.df[, check], ef.df[, check], egs.df[, check], emb.df[, check], gs.df[, check])
infected.chickens.fu <- which(apply(all.df, 1, function(x) any(x > 0 & x < 9999)))
n.infected.chickens.fu <- length(infected.chickens.fu)

## Loading estimated sickness levels from infected birds.
load("s-infected.RData")
test.cases <- s.infected.vec

compile("comparison.cpp")
dyn.load(dynlib("comparison"))

mo.data <- list(n_chickens = n.chickens,
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
                infected_chickens = as.numeric(logical.infected.chickens.et),
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
                       s = matrix(1, nrow = n.chickens, ncol = n.all.times),
                       u_1 = rep(1, n.chickens), u_2 = rep(1, n.chickens))

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

mo.fit <- nlminb(mo.obj$par, mo.obj$fn, mo.obj$gr)
## AIC (remember objective is already the NLL). 1740.224 is current best model.
2*mo.fit$objective + 2*length(mo.fit$par)

## Choosing which model to inspect.
mo.obj <- mo.obj
## Creating report.
mo.rep <- sdreport(mo.obj)
## Looking at fixed parameters.
summary(mo.rep, "fixed")
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
mu.s <- matrix(summary(mo.rep, "report")[rownames(summary(mo.rep, "report")) == "mu_s", 1], nrow = n.chickens)
## Collecting lambda.y for dwm.
lambda.y.rep <- matrix(summary(mo.rep, "report")[rownames(summary(mo.rep, "report")) == "lambda_y_mat", 1],
                       nrow = n.chickens, ncol = max(n.measurements))

## Plotting mu.s for each chicken.
plot(all.times, mu.s[1, ], type = "l", ylim = range(mu.s))
for (i in 1:n.chickens){
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

## Plot of s for ith chicken.
#pdf(height = 4, width = 6, file = "chicken-sickness.pdf")
par(mar = c(4, 4, 0, 0), oma = rep(0.1, 4))
plot(all.times, s[1, ],  ylim = range(s), type = "n", ylab = "Chicken sickness",
     xlab = "Days since first treatment")
cols <- c("blue", "red")
for (i in 1:n.chickens){
    lines(all.times, s[i, ], col = cols[logical.infected.chickens.et[i] + 1])
}
#dev.off()

## Estimating probabilities of nonzero counts for infected chickens
## post-treatment and at follow-up.

## Collecting estimated sickness values for infected chickens for the
## last two tests.
s.infected <- s[infected.chickens.fu, time.indices[6:7] + 1]
s.infected.vec <- c(s.infected)
save(s.infected.vec, file = "s-infected.RData")

## Collecting p.zero for end of treatment.
logit.p.zero <- matrix(summary(mo.rep, "report")[rownames(summary(mo.rep, "report")) == "logit_p_zero", 1],
                       nrow = length(mo.data$test_cases))
logit.p.zero.se <- matrix(summary(mo.rep, "report")[rownames(summary(mo.rep, "report")) == "logit_p_zero", 2],
                             nrow = length(mo.data$test_cases))
logit.p.zero.lower <- logit.p.zero - qnorm(0.975)*logit.p.zero.se
logit.p.zero.upper <- logit.p.zero + qnorm(0.975)*logit.p.zero.se

p.zero <- plogis(logit.p.zero)
p.zero.lower <- plogis(logit.p.zero.lower)
p.zero.upper <- plogis(logit.p.zero.upper)
methods <- c("DWM", "EF", "EGS", "EMB", "GS")
n.methods <- length(methods)
colnames(p.zero) <- colnames(p.zero.lower) <- colnames(p.zero.upper) <-
  methods
rownames(p.zero) <- rownames(p.zero.lower) <- rownames(p.zero.upper) <-
    rep(infected.chickens.fu, 2)
## Probability of detecting infection.
p.nonzero <- 1 - p.zero
## Lower- and upper-limits for confidence intervals.
p.nonzero.lower <- 1 - p.zero.upper
p.nonzero.upper <- 1 - p.zero.lower

## Plotting the confidence intervals.
#pdf(width = 8, height = 12, file = "ci-plot.pdf")
par(mar = c(4, 4, 0, 0), oma = rep(0.1, 4), yaxs = "i")
plot.new()
bird.spacing <- 20
time.spacing <- 10
method.spacing <- 1
bird.y.pos <- bird.spacing*(1:n.infected.chickens.fu) - bird.spacing/2
plot.window(xlim = c(0, 1), ylim = c(0, max(bird.y.pos) + bird.spacing/2))
cols <- brewer.pal(10, "Paired")
cols.et <- cols[c(TRUE, FALSE)]
cols.fu <- cols[c(FALSE, TRUE)]
abline(v = 1, lty = "dotted")
lines(rep(0, 2), y = c(0, bird.y.pos[n.infected.chickens.fu]), lty = "dotted")
for (i in 1:n.infected.chickens.fu){
    for (j in 1:2){
        for (k in 1:n.methods){
            y.val <- bird.y.pos[i] + c(-1, 1)[j]*time.spacing/2 + (k - 3)*method.spacing 
            p <- p.nonzero[rownames(p.nonzero) == infected.chickens.fu[i], k][j]
            p.lower <- p.nonzero.lower[rownames(p.nonzero) == infected.chickens.fu[i], k][j]
            p.upper <- p.nonzero.upper[rownames(p.nonzero) == infected.chickens.fu[i], k][j]
            col <- c(cols.et[k], cols.fu[k])[j]
            points(p, y.val, pch = 16, col = col)
            lines(c(p.lower, p.upper), rep(y.val, 2), col = col)
        }
    }
}
abline(h = bird.y.pos[-n.infected.chickens.fu] + bird.spacing/2)
box()
axis(1)
axis(2, bird.y.pos, paste("Bird", infected.chickens.fu, "\n", "\n"), tick = FALSE)
axis(2, rep(bird.y.pos, each = 2) + rep(c(-1, 1)*time.spacing/2, n.infected.chickens.fu),
     rep(c("ET", "FU"), n.infected.chickens.fu), tick = FALSE)
title(xlab = "Probability of detecting an MO organism")
legend("topleft", rev(methods), lty = rep(1, n.methods), pch = 16, col = rev(cols.fu), bty = "n",
       bg = "white", cex = 0.8)
#dev.off()
