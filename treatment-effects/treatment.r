library(TMB)
## Sorting out data.
dwm.df <- as.matrix(read.csv("../modality-comparison/data/mo-dwm.csv", header = FALSE))[, -1]
ef.df <- as.matrix(read.csv("../modality-comparison/data/mo-ef.csv", header = FALSE))[, -1]
egs.df <- as.matrix(read.csv("../modality-comparison/data/mo-egs.csv", header = FALSE))[, -1]
emb.df <- as.matrix(read.csv("../modality-comparison/data/mo-emb.csv", header = FALSE))[, -1]
gs.df <- as.matrix(read.csv("../modality-comparison/data/mo-gs.csv", header = FALSE))[, -1]
times <- c(0, 3, 4, 7, 9, 10, 34)
all.times <- seq(0, 34, by = 0.5)
n.all.times <- length(all.times)
time.indices <- which(all.times %in% times) - 1
n.birds <- nrow(dwm.df)
n.times <- length(times)
n.measurements <- c(rep(7, 6), rep(6, 2), rep(7, 8))
## Combining counts from all modalities.
combined.df <- dwm.df + ef.df + egs.df + emb.df + gs.df
combined.na.df <- combined.df
combined.df[is.na(combined.df)] <- 9999
## Figuring out which birds are still infected at the end of treatment.
infected.birds.et <- combined.df[, 6] > 0

## Loading in weight data.
weight.df <- as.matrix(read.csv("../modality-comparison/data/weight.csv", header = FALSE))
weight.times <- c(0, 7, 10, 34)
weight.start <- weight.df[, 1]
weight.end <- ifelse(is.na(weight.df[, 4]), weight.df[, 3], weight.df[, 4])

## Getting a vector of treatment endings.
end.treatment <- ifelse(infected.birds.et, 20, 10)

## Data to pass to TMB.
treatment.data <- list(n_birds = n.birds,
                       times = times,
                       all_times = all.times,
                       time_indices = time.indices,
                       end_treatment = end.treatment,
                       mo_counts = combined.df,
                       n_measurements = n.measurements,
                       weight_start = weight.start,
                       weight_end = weight.end,
                       cov_function = 1 ## 1 for exponential, 2 for squared exponential.
                       )

## Compiling TMB template.
compile("treatment.cpp")
dyn.load(dynlib("treatment"))

## NULL MODEL.
null.pars <- list(beta_0 = log(mean(combined.df[, 1])), # Mean at time = 0.
                  beta_1 = 0, # Effect of treatment.
                  beta_2 = 0, # Effect of coming off treatment for early finishers.
                  beta_3 = 0, # Effect of coming off treatment for late finishers.
                  beta_w = 0, # Fixed effect of weight.
                  beta_wt = 0, # Fixed effect of weight-treatment interaction.
                  beta_wc = 0, # Effect of weight change on final expected count.
                  alpha_1 = 0,
                  log_theta = 0,
                  log_sigma_s = 0,
                  log_rho_s = 0,
                  s = matrix(0, nrow = n.birds, ncol = n.all.times))
## Object for null model.
null.obj <- MakeADFun(data = treatment.data,
                      parameters = null.pars,
                      map = list(beta_1 = factor(NA),
                                 beta_2 = factor(NA),
                                 beta_3 = factor(NA),
                                 beta_w = factor(NA),
                                 beta_wt = factor(NA),
                                 beta_wc = factor(NA),
                                 alpha_1 = factor(NA),
                                 log_sigma_s = factor(NA),
                                 log_rho_s = factor(NA),
                                 s = factor(matrix(NA, nrow = n.birds, ncol = n.all.times))
                                 ),
                      DLL = "treatment", silent = TRUE)
null.obj$fn.use <- function(x = last.par[-random], ...){
    out <- null.obj$fn(x, ...)
    cat(x, out, "\n")
    out
}
## Fitting null model.
null.fit <- nlminb(null.obj$par, null.obj$fn.use, null.obj$gr)
null.sdrep <- sdreport(null.obj)
summary(null.sdrep, "fixed")
null.ests <- summary(null.sdrep, "fixed")[, 1]

## INCLUDING RANDOM EFFECTS.
re.pars <- list(beta_0 = null.ests["beta_0"], # Mean at time = 0.
                beta_1 = 0, # Effect of treatment.
                beta_2 = 0, # Effect of coming off treatment for early finishers.
                beta_3 = 0, # Effect of coming off treatment for late finishers.
                beta_w = 0, # Fixed effect of weight.
                beta_wt = 0, # Fixed effect of weight-treatment interaction.
                beta_wc = 0, # Effect of weight change on final expected count.
                alpha_1 = 0,
                log_theta = null.ests["log_theta"],
                log_sigma_s = 0,
                log_rho_s = 0,
                s = matrix(0, nrow = n.birds, ncol = n.all.times))
## Object for the random-effects model.
re.obj <- MakeADFun(data = treatment.data,
                    parameters = re.pars,
                    map = list(beta_1 = factor(NA),
                               beta_2 = factor(NA),
                               beta_3 = factor(NA),
                               beta_w = factor(NA),
                               beta_wt = factor(NA),
                               beta_wc = factor(NA),
                               alpha_1 = factor(NA)
                               ),
                    random = "s",
                    DLL = "treatment", silent = TRUE)
re.obj$fn.use <- function(x = last.par[-random], ...){
    out <- re.obj$fn(x, ...)
    cat(x, out, "\n")
    out
}
## Fitting random effects model.
re.fit <- nlminb(re.obj$par, re.obj$fn.use, re.obj$gr)
re.sdrep <- sdreport(re.obj)
summary(re.sdrep, "fixed")
re.ests <- summary(re.sdrep, "fixed")[, 1]

## INCLUDING FIXED TREATMENT EFFECTS.
tr.pars <-  list(beta_0 = re.ests["beta_0"], # Mean at time = 0.
                 beta_1 = 0, # Effect of treatment.
                 beta_2 = 0, # Effect of coming off treatment for early finishers.
                 beta_3 = 0, # Effect of coming off treatment for late finishers.
                 beta_w = 0, # Fixed effect of weight.
                 beta_wt = 0, # Fixed effect of weight-treatment interaction.
                 beta_wc = 0, # Effect of weight change on final expected count.
                 alpha_1 = 0,
                 log_theta = 0,
                 log_sigma_s = re.ests["log_sigma_s"],
                 log_rho_s = re.ests["log_rho_s"],
                 s = matrix(0, nrow = n.birds, ncol = n.all.times))
## Object for the treatment-effects model.
tr.obj <- MakeADFun(data = treatment.data,
                    parameters = tr.pars,
                    map = list(beta_w = factor(NA),
                               beta_wt = factor(NA),
                               beta_wc = factor(NA),
                               alpha_1 = factor(NA)
                               ),
                    random = "s",
                    DLL = "treatment", silent = TRUE)
tr.obj$fn.use <- function(x = last.par[-random], ...){
    out <- tr.obj$fn(x, ...)
    cat(x, out, "\n")
    out
}
## Fitting treatment-effects model.
tr.fit <- nlminb(tr.obj$par, tr.obj$fn.use, tr.obj$gr)
tr.sdrep <- sdreport(tr.obj)
summary(tr.sdrep, "fixed")
tr.ests <- summary(tr.sdrep, "fixed")[, 1]

## INCLUDING INTERCEPT INTERACTION EFFECT.
ii.pars <- list(beta_0 = tr.ests["beta_0"], # Mean at time = 0.
                beta_1 = tr.ests["beta_1"], # Effect of treatment.
                beta_2 = tr.ests["beta_2"], # Effect of coming off treatment for early finishers.
                beta_3 = tr.ests["beta_3"], # Effect of coming off treatment for late finishers.
                alpha_1 = 0,
                beta_w = 0, # Fixed effect of weight.
                beta_wt = 0, # Fixed effect of weight-treatment interaction.
                beta_wc = 0, # Effect of weight change on final expected count.
                log_theta = tr.ests["log_theta"],
                log_sigma_s = tr.ests["log_sigma_s"],
                log_rho_s = tr.ests["log_rho_s"],
                s = matrix(0, nrow = n.birds, ncol = n.all.times))
## Object for the intercept-interaction-effect model.
ii.obj <- MakeADFun(data = treatment.data,
                    parameters = ii.pars,
                    map = list(beta_w = factor(NA),
                               beta_wt = factor(NA),
                               beta_wc = factor(NA)
                               ),
                    random = "s",
                    DLL = "treatment", silent = TRUE)
ii.obj$fn.use <- function(x = last.par[-random], ...){
    out <- ii.obj$fn(x, ...)
    cat(x, out, "\n")
    out
}
## Fitting intercept-interaction model.
ii.fit <- nlminb(ii.obj$par, ii.obj$fn.use, ii.obj$gr)
ii.sdrep <- sdreport(ii.obj)
summary(ii.sdrep, "fixed")
ii.ests <- summary(ii.sdrep, "fixed")[, 1]

## INCLUDING WEIGHT-DIFFERENCE EFFECT FOR FINAL EXPECTATION.
wc.pars <- list(beta_0 = ii.ests["beta_0"], # Mean at time = 0.
                beta_1 = ii.ests["beta_1"], # Effect of treatment.
                beta_2 = ii.ests["beta_2"], # Effect of coming off treatment for early finishers.
                beta_3 = ii.ests["beta_3"], # Effect of coming off treatment for late finishers.
                alpha_1 = ii.ests["alpha_1"],
                beta_w = 0, # Fixed effect of weight.
                beta_wt = 0, # Fixed effect of weight-treatment interaction.
                beta_wc = 0, # Effect of weight change on final expected count.
                log_theta = ii.ests["log_theta"],
                log_sigma_s = ii.ests["log_sigma_s"],
                log_rho_s = ii.ests["log_rho_s"],
                s = matrix(0, nrow = n.birds, ncol = n.all.times))
## Object for the weight-difference effect model.
wc.obj <- MakeADFun(data = treatment.data,
                    parameters = wc.pars,
                    map = list(beta_w = factor(NA),
                               beta_wt = factor(NA)
                               ),
                    random = "s",
                    DLL = "treatment", silent = TRUE)
wc.obj$fn.use <- function(x = last.par[-random], ...){
    out <- wc.obj$fn(x, ...)
    cat(x, out, "\n")
    out
}
## Fitting weight-difference model.
wc.fit <- nlminb(wc.obj$par, wc.obj$fn.use, wc.obj$gr)
wc.sdrep <- sdreport(wc.obj)
summary(wc.sdrep, "fixed")
wc.ests <- summary(wc.sdrep, "fixed")[, 1]

## Holding beta_1 at 0 for likelihood-ratio test.
wc.hold.b1.pars <- wc.pars
wc.hold.b1.pars$beta_1 <- 0
wc.hold.b1.obj <- wc.obj <- MakeADFun(data = treatment.data,
                                      parameters = wc.hold.b1.pars,
                                      map = list(beta_1 = factor(NA),
                                                 beta_w = factor(NA),
                                                 beta_wt = factor(NA)
                                                 ),
                                      random = "s",
                                      DLL = "treatment", silent = TRUE)
wc.hold.b1.fit <- nlminb(wc.hold.b1.obj$par, wc.hold.b1.obj$fn, wc.hold.b1.obj$gr)

## Holding alpha_1 at 0 for likelihood-ratio test.
wc.hold.a1.pars <- wc.pars
wc.hold.a1.pars$alpha_1 <- 0
wc.hold.a1.obj <- wc.obj <- MakeADFun(data = treatment.data,
                                      parameters = wc.hold.a1.pars,
                                      map = list(beta_w = factor(NA),
                                                 beta_wt = factor(NA),
                                                 alpha_1 = factor(NA)
                                                 ),
                                      random = "s",
                                      DLL = "treatment", silent = TRUE)
wc.hold.a1.fit <- nlminb(wc.hold.a1.obj$par, wc.hold.a1.obj$fn, wc.hold.a1.obj$gr)

## INCLUDING WEIGHT EFFECT.
w.pars <- list(beta_0 = wc.ests["beta_0"], # Mean at time = 0.
               beta_1 = wc.ests["beta_1"], # Effect of treatment.
               beta_2 = wc.ests["beta_2"], # Effect of coming off treatment for early finishers.
               beta_3 = wc.ests["beta_3"], # Effect of coming off treatment for late finishers.
               alpha_1 = wc.ests["alpha_1"],
               beta_w = 0, # Fixed effect of weight.
               beta_wt = 0, # Fixed effect of weight-treatment interaction.
               beta_wc = wc.ests["beta_wc"], # Effect of weight change on final expected count.
               log_theta = wc.ests["log_theta"],
               log_sigma_s = wc.ests["log_sigma_s"],
               log_rho_s = wc.ests["log_rho_s"],
               s = matrix(0, nrow = n.birds, ncol = n.all.times))
## Object for the weight-effect model.
w.obj <- MakeADFun(data = treatment.data,
                   parameters = w.pars,
                   map = list(beta_wt = factor(NA)
                              ),
                   random = "s",
                   DLL = "treatment", silent = TRUE)
w.obj$fn.use <- function(x = last.par[-random], ...){
    out <- w.obj$fn(x, ...)
    cat(x, out, "\n")
    out
}
## Fitting weight-effects model.
w.fit <- nlminb(w.obj$par, w.obj$fn.use, w.obj$gr)
w.sdrep <- sdreport(w.obj)
summary(w.sdrep, "fixed")
w.ests <- summary(w.sdrep, "fixed")[, 1]

## INCLUDING WEIGHT INTERACTION EFFECT.
wi.pars <- list(beta_0 = w.ests["beta_0"], # Mean at time = 0.
                beta_1 = w.ests["beta_1"], # Effect of treatment.
                beta_2 = w.ests["beta_2"], # Effect of coming off treatment for early finishers.
                beta_3 = w.ests["beta_3"], # Effect of coming off treatment for late finishers.
                alpha_1 = w.ests["alpha_1"],
                beta_w = w.ests["beta_w"], # Fixed effect of weight.
                beta_wt = 0, # Fixed effect of weight-treatment interaction.
                beta_wc = w.ests["beta_wc"], # Effect of weight change on final expected count.
                log_theta = w.ests["log_theta"],
                log_sigma_s = w.ests["log_sigma_s"],
                log_rho_s = w.ests["log_rho_s"],
                s = matrix(0, nrow = n.birds, ncol = n.all.times))
## Object for the weight-interaction-effect model.
wi.obj <- MakeADFun(data = treatment.data,
                    parameters = wi.pars,
                    random = "s",
                    DLL = "treatment", silent = TRUE)
wi.obj$fn.use <- function(x = last.par[-random], ...){
    out <- wi.obj$fn(x, ...)
    cat(x, out, "\n")
    out
}
## Fitting weight-interaction model.
wi.fit <- nlminb(wi.obj$par, wi.obj$fn.use, wi.obj$gr)
wi.sdrep <- sdreport(wi.obj)
summary(wi.sdrep, "fixed")
wi.ests <- summary(wi.sdrep, "fixed")[, 1]

## Calculating log-likelihoods.
null.ll <- -null.fit$objective
re.ll <- -re.fit$objective
tr.ll <- -tr.fit$objective
ii.ll <- -ii.fit$objective
wc.ll <- -wc.fit$objective
wc.hold.b1.ll <- -wc.hold.b1.fit$objective
wc.hold.a1.ll <- -wc.hold.a1.fit$objective
w.ll <- -w.fit$objective
wi.ll <- -wi.fit$objective
## Calculating number of parameters.
null.npar <- length(null.fit$par)
re.npar <- length(re.fit$par)
tr.npar <- length(tr.fit$par)
ii.npar <- length(ii.fit$par)
wc.npar <- length(wc.fit$par)
wc.hold.b1.npar <- length(wc.hold.b1.fit$par)
wc.hold.a1.npar <- length(wc.hold.a1.fit$par)
w.npar <- length(w.fit$par)
wi.npar <- length(wi.fit$par)

## Calculating AICs.
null.AIC <- -2*null.ll + 2*null.npar
re.AIC <- -2*re.ll + 2*re.npar
tr.AIC <- -2*tr.ll + 2*tr.npar
ii.AIC <- -2*ii.ll + 2*ii.npar
wc.AIC <- -2*wc.ll + 2*wc.npar
wc.hold.b1.AIC <- -2*wc.hold.b1.ll + 2*wc.hold.b1.npar
wc.hold.a1.AIC <- -2*wc.hold.a1.ll + 2*wc.hold.a1.npar
w.AIC <- -2*w.ll + 2*w.npar
wi.AIC <- -2*wi.ll + 2*wi.npar

## Tabulating AICs.
AICs <- c(null.AIC, re.AIC, tr.AIC, ii.AIC, wc.AIC, w.AIC, wi.AIC)
npars <- c(null.npar, re.npar, tr.npar, ii.npar, wc.npar, w.npar, wi.npar)
mod.names <- c("null", "re", "tr", "ii", "wc", "w", "wi")
## Best model by AIC is wc.
data.frame(mod.names, AICs, npars)

## Hypothesis tests via likelihood ratio.

## Testing for presence of an original weight effect on expected
## number of MO shedded. Comopares models W and WC.
weight.lrts <- 2*(w.ll - wc.ll)
weight.p <- 1 - pchisq(weight.lrts, w.npar - wc.npar) # 0.585.

## Testing for presence of an original weight effect at all (both its
## direct effect on shedding, and its interaction with
## treatment). Compares models WI and WC.
joint.weight.lrts <- 2*(wi.ll - wc.ll)
joint.weight.p <- 1 - pchisq(joint.weight.lrts, wi.npar - wc.npar) # 0.165.

## Testing for a presence of initial shedding on the efficacy of the
## treatment. Compares models WC and WC.HOLD.A1.
shed.lrts <- 2*(wc.ll - wc.hold.a1.ll)
shed.p <- 1 - pchisq(shed.lrts, wc.npar - wc.hold.a1.npar) # 0.014

## Testing for a presence of a treatment effect on expected number of
## MO shedded. Compares models WC and WC.HOLD.B1.
treat.lrts <- 2*(wc.ll - wc.hold.b1.ll)
treat.p <- 1 - pchisq(treat.lrts, wc.npar - wc.hold.b1.npar) # <0.0001.

## Testing for presence of a weight-difference effect on the final
## reading at follow-up. Compares models WC and II.
weightchange.lrts <- 2*(wc.ll - ii.ll)
weightchange.p <- 1 - pchisq(weightchange.lrts, wc.npar - ii.npar) # 0.003.


## Collecting expectations.
model.fit <- wc.fit # Choose which model to use here. (Set at best model by AIC, wc.fit)
model.obj <- wc.obj # And here.
model.sdrep <- wc.sdrep # And here.
model.rep <- summary(model.sdrep, "report")
log.mu <- matrix(model.rep[rownames(model.rep) == "log_mu", 1], nrow = 16)

## Calculating the treatment effects for the birds.
t.effects <-  model.rep[rownames(model.rep) == "t_effect", 1]
t.effects.se <- model.rep[rownames(model.rep) == "t_effect", 2]
t.effects.ci.lower <- t.effects - qnorm(0.975)*t.effects.se
t.effects.ci.upper <- t.effects + qnorm(0.975)*t.effects.se
t.effects.ci <- cbind(t.effects.ci.lower, t.effects.ci.upper)
t.effects.perc <- 100*(exp(t.effects) - 1)
t.effects.perc.ci <- 100*(exp(t.effects.ci) - 1)

## Interpreting the effect of weight change on final reading.
beta.wc <- summary(model.sdrep, "fixed")["beta_wc", 1]
beta.wc.se <- summary(model.sdrep, "fixed")["beta_wc", 2]
beta.wc.ci <- beta.wc + c(-1, 1)*qnorm(0.975)*beta.wc.se
wc.effect <- 100*(exp(beta.wc) - 1)
wc.effect.ci <- 100*(exp(beta.wc.ci) - 1)


## Plotting fixed treatment effects.
t.effect <- model.rep[rownames(model.rep) == "t_effect", 1]
plot.new()
plot.window(xlim = range(all.times), ylim = (range(log.mu)))
box()
axis(1)
axis(2)
for (i in 1:n.birds){
    lines(all.times, (log.mu[i, ]))
}


## Plotting the data with expected values.
plot.new()
plot.window(xlim = range(all.times), ylim = c(0, max(combined.na.df, na.rm = TRUE)))
box()
axis(1)
axis(2)

for (i in 1:n.birds){
    lines(times, combined.na.df[i, ])
    lines(all.times, exp(log.mu[i, ]), col = "red")
}

## Choose a bird to just plot one.
i <- 1
plot(times, combined.na.df[i, ], xlim = range(all.times),
         ylim = c(0, max(c(exp(log.mu[i, ]), combined.na.df[i, ]), na.rm = TRUE)),
     col = "blue")
lines(all.times, exp(log.mu[i, ]), col = "red")

## Plotting the estimated expectations for each bird.
plot.new()
plot.window(xlim = range(all.times), ylim = exp(range(log.mu)))
box()
axis(1)
axis(2)
for (i in 1:n.birds){
    lines(all.times, exp(log.mu[i, ]))
}
save.image(file = "treatment.RData")
