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
                       cov_function = 1  ## 1 for exponential, 2 for squared exponential.
                       )
## Parameter starting values.
treatment.parameters <- list(beta_0 = log(mean(combined.df[, 1])), # Mean at time = 0.
                             beta_1 = 0.1, # Effect of treatment.
                             beta_2 = 0.1, # Effect of coming off treatment for early finishers.
                             beta_3 = 0.1, # Effect of coming off treatment for late finishers.
                             beta_w = 0.1, # Fixed effect of weight.
                             beta_wt = 0.1, # Fixed effect of weight-treatment interaction.
                             log_theta = log(5),
                             log_sigma_s = 0.1,
                             log_rho_s = log(1),
                             s = matrix(1, nrow = n.birds, ncol = n.all.times))

## Compiling TMB template.
compile("treatment.cpp")
dyn.load(dynlib("treatment"))

## Object for full model.
treatment.obj <- MakeADFun(data = treatment.data,
                           parameters = treatment.parameters,
                           #map = list(beta_w = factor(NA)),
                           random = "s",
                           DLL = "treatment")

## Fitting model.
treatment.fit <- nlminb(treatment.obj$par, treatment.obj$fn, treatment.obj$gr)
treatment.rep <- sdreport(treatment.obj)
summary(treatment.rep, "fixed")

## Collecting expectations.
summary.rep <- summary(treatment.rep, "report")
log.mu <- matrix(summary.rep[rownames(summary.rep) == "log_mu", 1], nrow = 16)

## Plotting the data.
plot.new()
plot.window(xlim = range(all.times), ylim = c(0, max(combined.na.df, na.rm = TRUE)))
box()
axis(1)
axis(2)

for (i in 1:n.birds){
    lines(times, combined.na.df[i, ])
    lines(all.times, exp(log.mu[i, ]), col = "red")
}

i <- 10
plot(times, combined.na.df[i, ], xlim = range(all.times),
         ylim = c(0, max(c(exp(log.mu[i, ]), combined.na.df[i, ]), na.rm = TRUE)),
     col = "blue")
lines(all.times, exp(log.mu[i, ]), col = "red")




plot.new()
plot.window(xlim = range(all.times), ylim = range(log.mu))
box()
axis(1)
axis(2)
for (i in 1:n.birds){
    lines(all.times, log.mu[i, ])
}
