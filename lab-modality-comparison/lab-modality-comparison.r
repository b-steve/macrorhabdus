library(epitools)
library(stringr)
## Loading in data frame.
df <- read.csv("lab-modality-comparison.csv")
n.birds <- nrow(df)

## Effect of sex on macro-positive.

## A loop that analyses effects of the following variables on
## macro-positive:
## - sex (cock or hen)
## - age (adult or juvenile)
## - gross findings (abnormal or normal feathers)
## - circovirus liver status (positive or negaive)
vars <- c("sex", "age", "gf", "cl")
varnames <- c("Sex", "Age", "Gross findings", "Circovirus liver status")
## Lists to store output.
tabs <- props <- chisqs <- fits <- ors <- vector("list", length(vars))
names(tabs) <- names(props) <- names(chisqs) <- names(fits) <- names(ors) <- vars
for (i in vars){
    var <- df[[i]]
    ## Creating contingency table.
    tabs[[i]] <- table(var, df$mp)
    ## Calculating sample proportions.
    props[[i]] <- tapply(df$mp == "positive", var, mean)
    ## Carrying out chi-squared test.
    chisqs[[i]] <- chisq.test(df$mp, var, simulate.p.value = TRUE, B = 10000)
    ## Fitting logistic regression model.
    fits[[i]] <- glm(df$mp == "positive" ~ var, family = "binomial")
    ## Estimating odds ratio.
    ors[[i]] <- suppressWarnings(oddsratio.small(table(var, df$mp)))
}

## A plot comparing the sample proportions.
pdf(file = "prop-mp.pdf", width = 9)
n.vars <- length(vars)
par(las = 1, yaxs = "i")
plot.new()
plot.window(xlim = c(1, 2*n.vars), ylim = c(0, 1))
box()
axis(2)
title(ylab = "Proportion macro positive")
for (i in 1:n.vars){
    ## Plotting sample proportions.
    points((i*2 - 1):(i*2), props[[i]], pch = 16)
    ## Separating varialbes with a vertical line.
    abline(v = i*2 + 0.5, lty = "dotted")
    ## Adding variable level names to x-axis.
    axis(1, at = c(i*2 - 1, i*2), labels = str_to_sentence(names(props[[i]])))
    ## Adding the variable name below the level names.
    mtext(varnames[i], 1, 3, at = i*2 - 0.5)
    ## Calculating a confidence interval for each level.
    newdata <- data.frame(var = names(props[[i]]))
    preds <- predict(fits[[i]], newdata = data.frame(var = names(props[[i]])), se.fit = TRUE)
    cis <- matrix(0, nrow = length(preds$fit), ncol = 2)
    for (j in 1:length(preds$fit)){
        cis[j, ] <- preds$fit[j] + c(-1, 1)*qnorm(0.975)*preds$se.fit[j]
    }
    cis <- plogis(cis)
    ## Plotting confidence intervals using lines.
    for (j in 1:nrow(cis)){
        lines(rep(i*2 - 2 + j, 2), cis[j, ])
    }
}
dev.off()

## Extracting p-values from the logistic regression models for each of the effects.
sapply(fits, function(x) summary(x)$coefficients[2, 4])
## We can also get p-values from odds-ratio tests, which are similar.
sapply(ors, function(x) x$p.value[2, 1])

## Comparisons between modalities.

## Extracting just the columns for the modalities being compared.
mod.df <- df[, c("mft", "pvf", "pcri", "pcrf", "hist")]
mod.names <- c("Macro float technique", "Proventricular float", "PCR (isthmus)", "PCR (faeces)", "Histopathology")
n.mods <- ncol(mod.df)
## Calculating sample positive proportions.
samp.props <- apply(mod.df, 2, function(x) mean(x == "positive"))
## Extracting just the columns for the modalities being compared.
mod.df <- df[, c("mft", "pvf", "pcri", "pcrf", "hist")]
mod.names <- c("Macro float technique", "Proventricular float", "PCR (isthmus)", "PCR (faeces)", "Histopathology")
n.mods <- ncol(mod.df)
## Calculating sample positive proportions.
samp.props <- apply(mod.df, 2, function(x) mean(x == "positive"))
## Calculating confidence intervals for the proportions.
prop.cis <- matrix(0, nrow = n.mods, ncol = 2)
for (i in 1:n.mods){
    prop.cis[i, ] <- plogis(confint(glm(mod.df[, i] == "positive" ~ 1, family = "binomial")))
}

## A plot for sample positive proportions.
pdf("mod-comparison.pdf", width = 9)
par(las = 1, yaxs = "i")
plot.new()
plot.window(xlim = c(1, n.mods), ylim = c(0, 1))
box()
axis(2)
axis(1, at = 1:n.mods, labels = mod.names)
points(1:n.mods, samp.props, pch = 16)
segments(1:5, prop.cis[, 1], 1:5, prop.cis[, 2])
dev.off()

## Functions for exact- and mid-p values for McNemar's test.
mcnemar.exact.p <- function(x){
    if (x[1, 2] == x[2, 1]){
        out <- 1 - 0.5*dbinom(x[1, 2], x[1, 2] + x[2, 1], 0.5)
    } else {
        out <- 2*pbinom(min(c(x[1, 2], x[2, 1])), x[1, 2] + x[2, 1], 0.5)
    }
    out
}

mcnemar.mid.p <- function(x){
    exact.p <- mcnemar.exact.p(x)
    exact.p - dbinom(min(c(x[1, 2], x[2, 1])), x[1, 2] + x[2, 1], 0.5)
}

## McNemar tests for each pair of modalities.
standard.p <- corrected.p <- exact.p <- mid.p <- boot.p <- matrix(NA, n.mods, n.mods)
rownames(standard.p) <- rownames(corrected.p) <- rownames(exact.p) <- rownames(mid.p) <- rownames(boot.p) <- names(mod.df)
colnames(standard.p) <- colnames(corrected.p) <- colnames(exact.p) <- colnames(mid.p) <- colnames(boot.p) <- names(mod.df)
for (i in 1:(n.mods - 1)){
    for (j in (i + 1):n.mods){
        ct <- table(mod.df[, i], mod.df[, j])
        standard.p[i, j] <- mcnemar.test(ct, correct = FALSE)$p.value
        corrected.p[i, j] <- mcnemar.test(ct, correct = TRUE)$p.value
        exact.p[i, j] <- mcnemar.exact.p(ct)
        mid.p[i, j] <- mcnemar.mid.p(ct)
    }
}

standard.p
corrected.p
exact.p
mid.p
