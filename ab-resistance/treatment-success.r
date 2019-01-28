## Importing the data.
treatment.df <- read.csv("data/treatment-success.csv")
colnames(treatment.df)[4] <- "success"
treatment.df <- treatment.df[treatment.df$Dose %in% c(25, 100), ]
success <- treatment.df$success
dose <- factor(treatment.df$Dose)
s.time <- treatment.df$Survival.Time
f.time <- treatment.df$Followup
time <- apply(cbind(s.time, f.time), 1, function(x) x[!is.na(x)])
event <- ifelse(is.na(s.time), 0, 1)

## Logistic regression.
fit.logistic <- glm(success ~ dose, family = "binomial")
summary(fit.logistic)
anova(fit.logistic, test = "Chisq")

## Survival analysis.
library(survival)
## Creating survival object.
bird.surv <- Surv(time = time, event = event)
dose <- factor(dose, levels = c(25, 100))
## Best-fitting parametric survival regression model.
fit <- survreg(bird.surv ~ dose, dist = "loglogistic")
summary(fit)
