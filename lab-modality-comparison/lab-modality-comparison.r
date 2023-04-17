## Loading in data frame.
df <- read.csv("lab-modality-comparison.csv")
n.birds <- nrow(df)

## Comparisons between modalities.

## Extracting just the columns for the modalities being compared.
mod.df <- df[, c("mft", "pvf", "pcri", "pcrf", "hist")]
mod.names <- c("Macro float technique", "Proventricular float", "PCR (isthmus)", "PCR (faeces)", "Histopathology")
n.mods <- ncol(mod.df)

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

## Estimates of sensitivity relative to gold standard.
pos.df <- mod.df[apply(mod.df, 1, function(x) any(x == "positive")), ]
pos.df <- data.frame(ifelse(pos.df == "positive", 1, 0))
library(DescTools)
BinomCI(apply(pos.df, 2, sum), nrow(pos.df), method = "wilson")
