## Shrinkage methods script

## load Chapter 8 functions
source("dper-chapter9-functions.R")

## Selection bias illustration
## data generating function
f <- function(n, beta.est = 0.8, mean = 0, sd = 1) {
    1 + (beta.est * seq(-1, 1, length = n)) +
        rnorm(n, mean = mean, sd = sd)
}

## generate 1000 data sets of size 10 accoring to f
set.seed(11)
dat <- lapply(seq_len(5000), function(x, ...) f(...), n = 10)

## fit a model to each data set
mods <- lapply(dat, function(y, x) lm(y ~ x),
               x = seq(-1, 1, length = 10))

## Get the coefficient for each model
coefs <- sapply(mods, function(x) coef(x)[2])

## get p values of slopes
p.vals <- sapply(mods, function(m) coef(summary(m))[2,4])

## select \hat{beta} if p-value <= 0.05, else 0
beta.hat <- ifelse(p.vals <= 0.05, coefs, 0)

## plot the coefficients
layout(1:2)
op <- par(mar = c(3,4,2,2) + 0.1, las = 1)
hist(coefs, freq = FALSE, main = "", xlab = "",
     breaks = seq(floor(min(coefs)), ceiling(max(coefs)), by = 0.1) + 0.05)
lines(density(coefs), lwd = 2)
abline(v = 0.8, lty = "dashed", col = "black")
text(x = 0.8, y = 0.9, labels = expression(beta == 0.8), xpd = NA)
## and the significant betas
hist(beta.hat, freq = FALSE, main = "", xlab = "",
     breaks = seq(floor(min(coefs)), ceiling(max(coefs)), by = 0.1) + 0.05)
lines(density(beta.hat), lwd = 2)
abline(v = 0.8, lty = "dashed", col = "black")
text(x = 0.8, y = 7.75, labels = expression(beta == 0.8), xpd = NA)
par(op)
layout(1)

## Some examples of shrinkage using the Ozone data from the MARS example
## load the 'earth' package
require("earth")

## load the Ozone data
data(ozone1)
## look at the first few data points
head(ozone1)

## load the glmnet package
require("glmnet")

## ridge regression - fudge this with a tiny alpha
mod.r <- glmnet(data.matrix(ozone1[,-1]), log(ozone1[,1]), alpha = 0.00001)
cv.r <- cv.glmnet(data.matrix(ozone1[,-1]), log(ozone1[,1]), alpha = 0.00001)
##mod.r

## lasso
mod.l <- glmnet(data.matrix(ozone1[,-1]), log(ozone1[,1]), alpha = 1)
cv.l <- cv.glmnet(data.matrix(ozone1[,-1]), log(ozone1[,1]), alpha = 1)
##mod.l

## elastic net
mod.el <- glmnet(data.matrix(ozone1[,-1]), log(ozone1[,1]), alpha = 0.5)
cv.el <- cv.glmnet(data.matrix(ozone1[,-1]), log(ozone1[,1]), alpha = 0.5)
##mod.el

## figures 9.27
layout(matrix(1:6, ncol = 2, byrow = TRUE))
mar <- c(5,4,4,6.7) + 0.1
## ridge
op <- par(mar = mar)
coefPlot(mod.r, label = TRUE,
         varNames = c("Pressure Height","Wind Speed","Humidity",
         "Temperature","Inversion BH", "Pressure Gradient",
         "Inversion Temp.", "Visibility", "Day of Year"),
         air = 2, col = "black", lty = 1, allDF = FALSE)
figCaption("a", "topleft", cex = 1.2, inset = 0.05, font = 2)
par(op)
cvPlot(cv.r)
## lasso
op <- par(mar = mar)
coefPlot(mod.l, label = TRUE,
         varNames = c("Pressure Height","Wind Speed","Humidity",
         "Temperature","Inversion BH", "Pressure Gradient",
         "Inversion Temp.", "Visibility", "Day of Year"),
         air = 2, col = "black", lty = 1, allDF = TRUE)
figCaption("b", "topleft", cex = 1.2, inset = 0.05, font = 2)
par(op)
cvPlot(cv.l)
## elastic net
op <- par(mar = mar)
coefPlot(mod.el, label = TRUE,
         varNames = c("Pressure Height","Wind Speed","Humidity",
         "Temperature","Inversion BH", "Pressure Gradient",
         "Inversion Temp.", "Visibility", "Day of Year"),
         air = 2, col = "black", lty = 1, allDF = TRUE)
figCaption("c", "topleft", cex = 1.2, inset = 0.05, font = 2)
par(op)
cvPlot(cv.el)
layout(1)
