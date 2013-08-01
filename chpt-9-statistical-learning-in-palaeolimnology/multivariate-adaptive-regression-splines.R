## Chapter 9 DPER Numerical Book
##
## MARS example

## load the 'earth' package
require("earth")

## load the Ozone data
data(ozone1)
## look at the first few data points
head(ozone1)

## fit a MARS model allowing one-way interactions
mod <- earth(O3 ~ . - doy, data = ozone1, degree = 2,
             glm = list(family = Gamma))
mod

## model summary, including the coefficients
summary(mod)

## Model diagnostics plot
plot(mod, legend.pos = "bottomright", caption = "")

## plot each fitted functions holding other terms at their medians
plotmo(mod, nrug = -1)

## plot variable importance
plot(evimp(mod), type.nsubsets = "b", x.legend = "topright",
     type.gcv = "b", type.rss = "b", col.gcv = "lightblue")

## Figure for paper: plot model terms
plotmo(mod, nrug = -1, caption = "",
       main = c("Temperature","DPG","Visibility","Wind speed:IBH",
       "Humidity:Temperature","Humidity:IBH","Temperature:DPG",
       "IBT:Visibility"), trace = FALSE)

## Figure for paper: Variable importance
varimp <- evimp(mod)
class(varimp) <- c(class(varimp), "matrix")
xvals <- seq_len(nrow(varimp))
varimp <- data.frame(varimp)
max.subsets <- varimp[1, "nsubsets"]
op <- par(mar = c(7.5,4,2,4) + 0.1)
plot(nsubsets ~ xvals, data = varimp, type = "n",
     ylim = c(0, max.subsets), axes = FALSE, ann = FALSE)
grid()
lines(nsubsets ~ xvals, data = varimp, type = "b", pch = 21,
      col = "black", bg = "white", lty = "dashed")
lines(max.subsets * rss/100 ~ xvals, data = varimp, type = "b", pch = 21,
      col = "black", bg = "grey", lty = "solid")
lines(max.subsets * gcv/100 ~ xvals, data = varimp, type = "b", pch = 21,
      col = "black", bg = "black", lty = "solid")
## axes
axis(side = 2, at = 0:9, las = 2)
labs <- c("Temperature", "Visibility", "Humidity", "IBH", "IBT", "DPG",
          "Wind speed")
labs <- paste(labs, varimp$col, sep = " - ")
axis(side = 1, at = xvals, labels = labs, las = 3)
axis(side = 4, at = seq(0, 1, by = 0.2) * max.subsets,
     labels = seq(0, 100, by = 20), las = 2)
box()
## titles
title(ylab = "Number of subsets")
mtext(side = 4, srt = 0, text = "Relative GCV or RSS", line = 2.8,
      cex = par("cex.axis"))
## legend
legend("topright", legend = c("N subsets", "RSS", "GCV"),
       pch = 21, pt.bg = c("white","grey","black"),
       lty = c("dashed", rep("solid", 2)),
       bty = "n")
par(op)

## Figure showing reflected pairs
pos.pair <- function(x, t) ifelse(x >= t, x - t, 0)
neg.pair <- function(x, t) ifelse(x <= t, t - x, 0)
x <- seq(0, 1, by = 0.1)
op <- par(mar = c(5,4,2,2) + 0.1)
plot(pos.pair(x, t = 0.5) ~ x, type = "l", ylim = c(0,0.6),
     ylab = expression(y), xlab = expression(x), axes = FALSE)
lines(neg.pair(x, t = 0.5) ~ x, type = "l", lty = "dashed")
axis(side = 1, at = seq(0, 1, by = 0.25))
axis(side = 2, at = seq(0, 0.7, by = 0.2), las = 2)
rug(x = 0.5, lwd = 1, ticksize = 0.02)
box()
legend("topright",
       legend = c(expression(group("(", x - t, ")")["+"]),
       expression(group("(", t - x, ")")["+"])),
       lty = c("solid", "dashed"), bty = "n")
par(op)
