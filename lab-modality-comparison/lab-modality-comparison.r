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

## Extracting p-values from the logistic regression models for each of the effects.
sapply(fits, function(x) summary(x)$coefficients[2, 4])
## We can also get p-values from odds-ratio tests, which are similar.
sapply(ors, function(x) x$p.value[2, 1])

## Comparisons between modalities.

## Extracting just the columns for the modalities being compared.
mod.df <- df[, c("mft", "pvf", "pcri", "pcrf", "hist")]
n.mods <- ncol(mod.df)
## Calculating sample positive proportions.
samp.props <- apply(mod.df, 2, function(x) mean(x == "positive"))
## Creating bootstrap data sets.
n.boots <- 10000
boots.df <- vector("list", n.boots)
for (i in 1:n.boots){
    ## Sampling rows with replacement.
    boots.df[[i]] <- mod.df[sample(n.birds, replace = TRUE), ]
}

## Positive proportions for each modality for each bootstrap data set.
boot.props <- sapply(boots.df, function(x) apply(x, 2, function(x) mean(x == "positive")))
## Bootstrap percentile CIs for each modality's positive proportion.
prop.cis <- t(apply(boot.props, 1, quantile, probs = c(0.025, 0.975)))

## A plot for sample positive proportions.
par(las = 1, yaxs = "i")
plot.new()
plot.window(xlim = c(1, n.mods), ylim = c(0, 1))
box()
axis(2)
axis(1, at = 1:n.mods, labels = names(mod.df))
points(1:n.mods, samp.props, pch = 16)
segments(1:5, prop.cis[, 1], 1:5, prop.cis[, 2])

## A function for pairwise differences.
pairwise.diffs <- function(x){
    out <- outer(x, x, `-`)
    out <- out[lower.tri(out)]
    if (!is.null(names(x))){
        out.names <- outer(names(x), names(x), function(x1, x2) paste0(x1, " - ", x2))
        names(out) <- out.names[lower.tri(out.names)]
    }
    out
}
## Calculating differences between sample proportions.
samp.diffs <- pairwise.diffs(samp.props)
## Differences between modalities for each boostrap data set.
boot.diffs <- apply(boot.props, 2, pairwise.diffs)
## Bootstrap percentile CIs for each difference.
diff.cis <- t(apply(boot.diffs, 1, quantile, probs = c(0.025, 0.975)))
diff.cis
## Proportion of bootstrapped estimates that are negative.
boot.diff.neg <- apply(boot.diffs, 1, function(x) mean(x < 0))
## Proportion of bootstrapped estimates that are equal to zero.
boot.diff.zero <- apply(boot.diffs, 1, function(x) mean(x == 0))
## Proportion of bootstrapped estimates that are positive.
boot.diff.pos <- apply(boot.diffs, 1, function(x) mean(x > 0))
## The p-value is twice the tail probability equal to or beyond zero,
## then doubled for a two-tailed test.
boot.diff.ps <- 2*(1 - apply(cbind(boot.diff.neg, boot.diff.pos), 1, max))
## Calculation gets it wrong for pcri vs pvf because there is no
## "tail", given that both vectors are exactly the same.
boot.diff.ps[boot.diff.ps > 1] <- 1
boot.diff.ps
