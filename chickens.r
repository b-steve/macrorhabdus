## Importing the data.
chickens.df <- read.csv("chickens-all.csv")
colnames(chickens.df) <- c("group", "cage", "bird", "weight", "slide", "mo")

## Sorting out data.
n <- length(unique(chickens.df$bird))
mo <- matrix(0, nrow = n, ncol = 10)
group <- numeric(n)
weight <- numeric(n)
for (i in 1:n){
    for (j in 1:10){
        mo[i, j] <- (chickens.df$mo[chickens.df$bird == unique(chickens.df$bird)[i]])[j]
    }
    if ((chickens.df$group[chickens.df$bird == unique(chickens.df$bird)[i]])[1] == "25mg/kg"){
        group[i] <- 1
    } else if ((chickens.df$group[chickens.df$bird == unique(chickens.df$bird)[i]])[1] == "100mg/kg"){
        group[i] <- 2
    }
    weight[i] <- (chickens.df$weight[chickens.df$bird == unique(chickens.df$bird)[i]])[1]/1000
}
## Average number of organisms per slide for each bird.
mean.mo <- apply(mo, 1, mean)
## Points at which the fitted lines are calculated for the plot.
weight.points <- seq(0.3, 0.6, length.out = 1000)

library(TMB)
## Compiling and loading TMB DLLs.
compile("tmb/chickens.cpp")
dyn.load(dynlib("chickens"))
compile("tmb/chickens-hold.cpp")
dyn.load(dynlib("chickens-hold"))

## Data to pass to TMB.

## Full model.
f.full <- MakeADFun(data = list(n = n,
                                mo_count = mo,
                                weight = weight,
                                group = group,
                                weight_points = weight.points),
                    parameters = list(beta_0 = 1,
                                      beta_weight = 0,
                                      beta_weight_2 = 0,
                                      beta_25 = 0,
                                      beta_100 = 0,
                                      log_sigma_u = log(1),
                                      log_sigma_v = log(1),
                                      u = numeric(n),
                                      v = matrix(0, nrow = n, ncol = 10)),
                    random = c("u", "v"), DLL = "chickens")

## Model with no quadratic weight term.
f.no_weight2 <- MakeADFun(data = list(n = n,
                                     mo_count = mo,
                                     weight = weight,
                                     group = group,
                                     weight_points = weight.points),
                         parameters = list(beta_0 = 1,
                                           beta_weight = 0,
                                           beta_weight_2 = 0,
                                           beta_25 = 0,
                                           beta_100 = 0,
                                           log_sigma_u = log(1),
                                           log_sigma_v = log(1),
                                           u = numeric(n),
                                           v = matrix(0, nrow = n, ncol = 10)),
                         map = list(beta_weight_2 = factor(NA)),
                         random = c("u", "v"), DLL = "chickens")

## Model with no weight terms.
f.no_weight <- MakeADFun(data = list(n = n,
                                     mo_count = mo,
                                     weight = weight,
                                     group = group,
                                     weight_points = weight.points),
                         parameters = list(beta_0 = 1,
                                           beta_weight = 0,
                                           beta_weight_2 = 0,
                                           beta_25 = 0,
                                           beta_100 = 0,
                                           log_sigma_u = log(1),
                                           log_sigma_v = log(1),
                                           u = numeric(n),
                                           v = matrix(0, nrow = n, ncol = 10)),
                         map = list(beta_weight_2 = factor(NA), beta_weight = factor(NA)),
                         random = c("u", "v"), DLL = "chickens")

## Model with no treatment group terms.
f.no_group <- MakeADFun(data = list(n = n,
                                    mo_count = mo,
                                    weight = weight,
                                    group = group,
                                    weight_points = weight.points),
                        parameters = list(beta_0 = 1,
                                          beta_weight = 0,
                                          beta_weight_2 = 0,
                                          beta_25 = 0,
                                          beta_100 = 0,
                                          log_sigma_u = log(1),
                                          log_sigma_v = log(1),
                                          u = numeric(n),
                                          v = matrix(0, nrow = n, ncol = 10)),
                        map = list(beta_25 = factor(NA), beta_100 = factor(NA)),
                        random = c("u", "v"), DLL = "chickens")

## Model restricting the effects of NTC and 25mg/kg groups to be the same.
f.hold_ntc.25 <- MakeADFun(data = list(n = n,
                                       mo_count = mo,
                                       weight = weight,
                                       group = group,
                                       weight_points = weight.points),
                           parameters = list(beta_0 = 1,
                                             beta_weight = 0,
                                             beta_weight_2 = 0,
                                             beta_25 = 0,
                                             beta_100 = 0,
                                             log_sigma_u = log(1),
                                             log_sigma_v = log(1),
                                             u = numeric(n),
                                             v = matrix(0, nrow = n, ncol = 10)),
                           map = list(beta_25 = factor(NA)),
                           random = c("u", "v"), DLL = "chickens")

## Model restricting the effects of NTC and 100mg/kg groups to be the same.
f.hold_ntc.100 <- MakeADFun(data = list(n = n,
                                        mo_count = mo,
                                        weight = weight,
                                        group = group,
                                        weight_points = weight.points),
                            parameters = list(beta_0 = 1,
                                              beta_weight = 0,
                                              beta_weight_2 = 0,
                                              beta_25 = 0,
                                              beta_100 = 0,
                                              log_sigma_u = log(1),
                                              log_sigma_v = log(1),
                                              u = numeric(n),
                                              v = matrix(0, nrow = n, ncol = 10)),
                            map = list(beta_100 = factor(NA)),
                            random = c("u", "v"), DLL = "chickens")

## Model restricting the effects of the 25mg/kg and 100mg/kg groups to be the same.
f.hold_25.100 <- MakeADFun(data = list(n = n,
                                        mo_count = mo,
                                        weight = weight,
                                        group = group,
                                        weight_points = weight.points),
                            parameters = list(beta_0 = 1,
                                              beta_weight = 0,
                                              beta_weight_2 = 0,
                                              beta_g = 0,
                                              log_sigma_u = log(1),
                                              log_sigma_v = log(1),
                                              u = numeric(n),
                                              v = matrix(0, nrow = n, ncol = 10)),
                           random = c("u", "v"), DLL = "chickens_hold")

## Fitting full and restricted models.
fit.full <- nlminb(f.full$par, f.full$fn, f.full$gr)
fit.no_weight2 <- nlminb(f.no_weight2$par, f.no_weight2$fn, f.no_weight2$gr)
fit.no_weight <- nlminb(f.no_weight$par, f.no_weight$fn, f.no_weight$gr)
fit.no_group <-  nlminb(f.no_group$par, f.no_group$fn, f.no_group$gr)
fit.hold_ntc.25 <- nlminb(f.hold_ntc.25$par, f.hold_ntc.25$fn, f.hold_ntc.25$gr)
fit.hold_ntc.100 <- nlminb(f.hold_ntc.100$par, f.hold_ntc.100$fn, f.hold_ntc.100$gr)
fit.hold_25.100 <- nlminb(f.hold_25.100$par, f.hold_25.100$fn, f.hold_25.100$gr)

## AICs
2*length(fit.full$par) + 2*fit.full$objective
2*length(fit.no_weight2$par) + 2*fit.no_weight2$objective
2*length(fit.no_weight$par) + 2*fit.no_weight$objective
2*length(fit.no_group$par) + 2*fit.no_group$objective

## Saving maximised log-likelihoods.
loglik.full <- -fit.full$objective
loglik.no_weight2 <- -fit.no_weight2$objective
loglik.no_weight <- -fit.no_weight$objective
loglik.no_group <- -fit.no_group$objective
loglik.hold_ntc.25 <- -fit.hold_ntc.25$objective
loglik.hold_ntc.100 <- -fit.hold_ntc.100$objective
loglik.hold_25.100 <- -fit.hold_25.100$objective

## Likelihood-ratio tests:
## H0: No weight effect at all.
lrts <- 2*(loglik.full - loglik.no_weight)
1 - pchisq(lrts, 2)
## H0: No quadratic (i.e., linear) weight effect.
lrts <- 2*(loglik.full - loglik.no_weight2)
1 - pchisq(lrts, 1)
## H0: No group effect.
lrts <- 2*(loglik.full - loglik.no_group)
1 - pchisq(lrts, 2)
## H0: NTC = 25mg/kg.
lrts <- 2*(loglik.full - loglik.hold_ntc.25)
1 - pchisq(lrts, 1)
## H0: NTC = 100mg/kg.
lrts <- 2*(loglik.full - loglik.hold_ntc.100)
1 - pchisq(lrts, 1)
## H0: 25mg/kg = 100mg/kg.
lrts <- 2*(loglik.full - loglik.hold_25.100)
1 - pchisq(lrts, 1)

## Getting fitted lines for the plot.
summary(rep.full, c("fixed"))
report.summary <- summary(rep.full, "report")
ntc.summary <- report.summary[rownames(report.summary) == "ntc_preds", ]
d25.summary <- report.summary[rownames(report.summary) == "d25_preds", ]
d100.summary <- report.summary[rownames(report.summary) == "d100_preds", ]

## Making the plot. The commented lines provide 95% CIs, but makes the plot messy.
par(xaxs = "i")
plot.new()
plot.window(xlim = c(0.35, 0.60), ylim = c(0, 11.1))
box()
axis(1)
axis(2)
abline(h = 0, col = "grey")
lines(weight.points, exp(ntc.summary[, 1]))
#lines(weight.points, exp(ntc.summary[, 1] + 1.96*ntc.summary[, 2]), lty = "dotted")
#lines(weight.points, exp(ntc.summary[, 1] - 1.96*ntc.summary[, 2]), lty = "dotted")
lines(weight.points, exp(d25.summary[, 1]), lty = "dashed")
#lines(weight.points, exp(d25.summary[, 1] + 1.96*d25.summary[, 2]), lty = "dotted", col = "blue")
#lines(weight.points, exp(d25.summary[, 1] - 1.96*d25.summary[, 2]), lty = "dotted", col = "blue")
lines(weight.points, exp(d100.summary[, 1]), lty = "dotted")
#lines(weight.points, exp(d100.summary[, 1] + 1.96*d100.summary[, 2]), lty = "dotted", col = "red")
#lines(weight.points, exp(d100.summary[, 1] - 1.96*d100.summary[, 2]), lty = "dotted", col = "red")
points(weight, mean.mo, pch = c(1, 2, 5)[group + 1])
title(xlab = "Bird weight (kg)", ylab = "Expected MO count per slide for the average bird")
legend(x = "topright", legend = c("NTC",  "25mg/kg", "100mg/kg"), col = "black", lty = c("solid", "dashed", "dotted"), pch = c(1, 2, 5))

## Getting some profile-likelihood CIs.
prof.beta_0 <- tmbprofile(f.full, "beta_0")
prof.beta_weight <- tmbprofile(f.full, "beta_weight")
prof.beta_weight_2 <- tmbprofile(f.full, "beta_weight_2")
prof.beta_NTC.25 <- tmbprofile(f.full, "beta_25")
prof.beta_NTC.100 <- tmbprofile(f.full, "beta_100")
prof.beta_25.100 <- tmbprofile(f.full, name = "beta_compare", lincomb = c(0, 0, 0, -1, 1, 0, 0))
prof.log_sigma_u <- tmbprofile(f.full, "log_sigma_u")
prof.log_sigma_v <- tmbprofile(f.full, "log_sigma_v")

## Displaying the profile-likelihood CIs.
rbind(confint(prof.beta_0),
      confint(prof.beta_weight),
      confint(prof.beta_weight_2),
      confint(prof.beta_NTC.25),
      confint(prof.beta_NTC.100),
      confint(prof.beta_25.100),
      exp(confint(prof.log_sigma_u)),
      exp(confint(prof.log_sigma_v)))

## Note that it should be possible to fit the above models using
## glmer(), but the full model suffers from nonconvergence for some
## reason. TMB appears to be a more stable optimiser.
library(lme4)
mo <- chickens.df$mo
group <- factor(chickens.df$group)
group <- relevel(group, ref = "NTC")
bird <- as.character(chickens.df$bird)
slide <- paste(chickens.df$bird, ".", chickens.df$slide, sep = "")
weight <- chickens.df$weight/1000
## This is equivalent to fit.no_weight2.
fit.glmer.linear <- glmer(mo ~ group + weight + (1 | bird) + (1 | slide), family = "poisson")
## This is equivalent to fit.full---doesn't converge.
fit.glmer.quad <- glmer(mo ~ group + weight + I(weight^2) + (1 | bird) + (1 | slide), family = "poisson")

## Comparing weights of infected and uninfected chickens.
chickens_ind.df <- within(chickens_ind.df, {infected <- ifelse(group == "PC", "no", "yes")})
weight <- chickens_ind.df$weight
infected <- chickens_ind.df$infected
## Plotting data.
boxplot(weight ~ infected)
## Equality of variance is a little questionable.
tapply(weight, infected, sd)
## Fit a model assuming equal variance.
fit <- lm(weight ~ infected)
summary(fit)
confint(fit)

## Do a randomisation test---this is robust to unequal variances.
## Getting the observed difference.
obs <- diff(tapply(weight, infected, mean))
## Getting 10000 differences after randomisation.
diffs <- numeric(10000)
for (i in 1:10000){
    diffs[i] <- diff(tapply(weight, sample(infected), mean))
}
## The below p-value very similar to the linear model; we'll just
## report that.
mean(abs(diffs) >= abs(obs))
