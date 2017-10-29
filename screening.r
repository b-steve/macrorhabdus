## Importing the data.
screening.df <- read.csv("data/screening.csv")

library(lubridate)
## Extracting variables.
date <- dmy(screening.df$date)
lovebird <- screening.df$lovebird
budgie <- screening.df$budgie

## Making the plot.
par(mfrow = c(2, 1), mar = c(0, 4, 0.5, 0), oma = c(4, 0, 0, 0) + 0.1, las = 1)
## Making the lovebird half.
plot.new()
plot.window(xlim = range(date), ylim = c(0, max(lovebird, na.rm = TRUE)))
## Limits of grey box.
treat.start <- date[18] - 1
treat.end <- date[27] + 1
## Making the grey box.
rect(treat.start, -20, treat.end, 100, col = "grey85", border = NA)
box()
## Horizontal line at 0.
abline(h = 0, lty = "dashed")
## Plotting lovebird data.
lines(date[!is.na(lovebird)], lovebird[!is.na(lovebird)])
points(date[!is.na(lovebird)], lovebird[!is.na(lovebird)], pch = 20)
## Adding lovebird label.
text(par("usr")[1] + grconvertX(0.1, from = "inches", to = "user") - grconvertX(0, from = "inches", to = "user"),
     par("usr")[4] - (grconvertY(0.1, from = "inches", to = "user") - grconvertY(0, from = "inches", to = "user")),
     labels = "Lovebird", adj = c(0, 1))
axis(2)
## Adding y-axis label.
par(las = 0)
mtext("MO/HPF", side = 2, line = 2.5)
par(las = 1)
## Making the budgie half.
plot.new()
plot.window(xlim = range(date), ylim = c(0, max(budgie, na.rm = TRUE)))
## Making the grey box.
rect(treat.start, -20, treat.end, 100, col = "grey85", border = NA)
box()
## Horizontal line at 0.
abline(h = 0, lty = "dashed")
## Plotting budgie data.
lines(date, budgie)
points(date, budgie, pch = 20)
axis(1, at = pretty(date)[2:6], labels = c("Feb", "Mar", "Apr", "May", "Jun"))
axis(2)
## Adding x-axis label.
mtext("Date", side = 1, line = 2.5)
## Adding budgie label.
text(par("usr")[1] + grconvertX(0.1, from = "inches", to = "user") - grconvertX(0, from = "inches", to = "user"),
     par("usr")[4] - (grconvertY(0.1, from = "inches", to = "user") - grconvertY(0, from = "inches", to = "user")),
     labels = "Budgerigar", adj = c(0, 1))
## Adding y-axis label.
par(las = 0)
mtext("MO/HPF", side = 2, line = 2.5)

## Permutation test

## Turning dates into day progression.
day <- as.vector(date - date[1]) + 1
## Lovebird variables.
lovebird <- screening.df$lovebird
day.l <- day[!is.na(lovebird)]
lovebird <- lovebird[!is.na(lovebird)]
lovebird.std <- lovebird/sd(lovebird)
## Budgie variables.
budgie <- screening.df$budgie
budgie.std <- budgie/sd(budgie)
day.b <- day[!is.na(budgie)]
## Treatment day limits.
treat.lims <- day[c(18, 27)]
## All possible start- and end-date pairings.
start.perm <- 1:(max(day.l) - 31)
end.perm <- start.perm + 31

## Observed means for lovebird.
obs.l <- mean(lovebird[day.l >= treat.lims[1] & day.l <= treat.lims[2]])
obs.l.std <- mean(lovebird.std[day.l >= treat.lims[1] & day.l <= treat.lims[2]])
## Doing permutations for lovebird.
n.l <- length(start.perm)
means.l <- numeric(n.l)
means.l.std <- numeric(n.l)
for (i in 1:n.l){
    means.l[i] <- mean(lovebird[day.l >= start.perm[i] & day.l <= end.perm[i]])
    means.l.std[i] <- mean(lovebird.std[day.l >= start.perm[i] & day.l <= end.perm[i]])
}
## Getting rid of periods with no observations.
means.l <- means.l[is.finite(means.l)]
means.l.std <- means.l.std[is.finite(means.l.std)]
## Individual rank and p-value for lovebird.
sum(means.l <= obs.l)
mean(means.l <= obs.l)
## Standardised results should be the same.
sum(means.l.std <= obs.l.std)
mean(means.l.std <= obs.l.std)

## Observed means for budgie.
obs.b <- mean(budgie[day.b >= treat.lims[1] & day.b <= treat.lims[2]])
obs.b.std <- mean(budgie.std[day.b >= treat.lims[1] & day.b <= treat.lims[2]])
## Doing permutations for budgie.
n.b <- length(start.perm)
means.b <- numeric(n.b)
means.b.std <- numeric(n.b)
for (i in 1:n.b){
    means.b[i] <- mean(budgie[day.b >= start.perm[i] & day.b <= end.perm[i]])
    means.b.std[i] <- mean(budgie.std[day.b >= start.perm[i] & day.b <= end.perm[i]])
}
## Getting rid of periods with no observations.
means.b <- means.b[is.finite(means.b)]
means.b.std <- means.b.std[is.finite(means.b.std)]
## Individual rank and p-value for budgie.
sum(means.b <= obs.b) 
mean(means.b <= obs.b)
## Standardised results should be the same.
sum(means.b.std <= obs.b.std) 
mean(means.b.std <= obs.b.std)

## Combining via standardised sums of means. Using standardised values
## is better---otherwise the lovebird data dominates the analysis, as
## these have higher variance.

## Outer product of possible sums across all 31-day periods.
m <- outer(means.l.std, means.b.std, `+`)
## Exact p-value.
mean(m <= obs.l.std + obs.b.std)
