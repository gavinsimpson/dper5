## Script to reproduce the principal curves analysis of the Abernethy
## Forest data set in Chapter 9

## Abernethy data
require(analogue)
data(abernethy)

## Stratigraphic plot of the Abernethy data, main taxa only
(BAR <- Stratiplot(Age ~ . - Depth, data =
                   chooseTaxa(abernethy, max.abun = 15, n.occ = 10),
                   type = c("h","g","l"), sort = "wa"))

## Drop the dating-related information
abernethy2 <- abernethy[, -(37:38)]

## Fit PCA...
aber.pca <- rda(abernethy2)

## ...and CA, and ...
aber.ca <- cca(abernethy2)

## ...a PrCurve to the data. First curve has same DF for each species,
aber.pc <- prcurve(abernethy2, method = "ca", trace = TRUE, plotit = TRUE,
                   vary = FALSE, penalty = 1.4)
## ...whilst this curve allows for different complexities per species, and
## this is the curve that we'll use for the rest of the script
aber.pc2 <- prcurve(abernethy2, method = "ca", trace = TRUE, plotit = TRUE,
                    vary = TRUE, penalty = 1.4)

## Plot the fitted curve
plot(aber.pc2)

## Plot the fitted curve - BW as per chapter
plot(aber.pc2, col.seg = "darkgrey", col = "black")

## left p[anel of Figure 9.24
layout(matrix(1:2, ncol = 2))
plot(aber.pc2)
#plot(gradientDist(aber.pc2))
plot(x = seq_len(nrow(abernethy2) - 1),
     y = unclass(abs(diff(gradientDist(aber.pc2)))), type = "h")
layout(1)

dis <- gradientDist(aber.pc)

## Reproduces Figure 9.24 right panel - correctly!
## note we plot 1 - gradientDist() to get all curves pointing in same
## direction, though the vagaries of eigenvector signs may mean
## this works one my machine and not yours!
Depth <- abernethy$Depth
Age <- abernethy$Age
plot(1 - gradientDist(aber.pc2), orderBy = Age,
     xlim = rev(range(Age)),
     type = "o", flipAxes = TRUE, xlab = "Age (Radiocarbon years BP)",
     main = "Abernethy Forest",
     cex = 0.8, pch = 21, col = "black", bg = "black")
lines(gradientDist(aber.pca), orderBy = Age,
      lty = "dashed", flipAxes = TRUE, type = "o",
      cex = 0.8, pch = 22, col = "red", bg = "red")
lines(1 - gradientDist(aber.ca), orderBy = Age,
      lty = "dotdash", flipAxes = TRUE, type = "o",
      cex = 0.8, pch = 23, col = "forestgreen", bg = "forestgreen")
legend("bottomleft", bty = "n", pch = 21:23,
       legend = c("PCurve","PCA","CA"),
       lty = c("solid","dashed","dotdash"),
       col = c("black","red","forestgreen"),
       pt.bg = c("black","red","forestgreen"), inset = 0.01)

## Reproduces Figure 9.24 - correctly!
Depth <- abernethy$Depth
Age <- abernethy$Age
layout(matrix(c(1,2,2), ncol = 3))
nseg <- nrow(abernethy2) - 1
RoC <- unclass(abs(diff(gradientDist(aber.pc2)))) / diff(Age)
plot(y = Age[-length(Age)],
     x = RoC * 1000, type = "n",
     ylim = rev(range(Age)),
     ylab = "Age (Radiocarbon years BP)",
     xlab  = expression("Rate of Change" ~ (kyr^{-1})),
     main = "a", cex.main = 1.5)
segments(x0 = rep(0, nseg), y0 = Age[-length(Age)],
         x1 = RoC * 1000, y1 = Age[-length(Age)])
plot(1 - gradientDist(aber.pc2), orderBy = Age,
     xlim = rev(range(Age)),
     type = "o", flipAxes = TRUE, xlab = "Age (Radiocarbon years BP)",
     main = "b",
     cex = 0.8, pch = 21, col = "black", bg = "black", cex.main = 1.5)
lines(gradientDist(aber.pca), orderBy = Age,
      lty = "dashed", flipAxes = TRUE, type = "o",
      cex = 0.8, pch = 22, col = "black", bg = "black")
lines(1 - gradientDist(aber.ca), orderBy = Age,
      lty = "dotdash", flipAxes = TRUE, type = "o",
      cex = 0.8, pch = 23, col = "black", bg = "black")
legend("topright", bty = "n", pch = 21:23,
       legend = c(expression(PCurve), expression(PCA[1]), expression(CA[1])),
       lty = c("solid","dashed","dotdash"),
       col = c("black","black","black"),
       pt.bg = c("black","black","black"), inset = 0.01, cex = 1.4,
       seg.len = 4)
layout(1)

## reproduces Figure 9.25
taxaWant <- tran(chooseTaxa(abernethy2, max.abun = 25, n.occ = 10),
                 method = "pcent2prop")
ylim <- c(0, max(apply(taxaWant, 2, max)))
p <- with(aber.pc2, data.frame(gradient = seq(min(lambda), max(lambda),
                               length = 1000)))
interpCurve <- approxfun(with(aber.pc2, lambda[tag]), Age)
p <- transform(p, AgeGradient = interpCurve(gradient))
layout(matrix(1:NCOL(taxaWant), ncol = 3))
op <- par(mar = c(2,1,4.1,2.1), oma = c(3.1,3.1,0,0), las = 1)
for(i in seq_len(NCOL(taxaWant))) {
    d <- with(aber.pc2, data.frame(gradient = lambda[tag],
                                   abundance = taxaWant[,i],
                                   AgeGradient = interpCurve(lambda[tag])))
    mod <- with(d,
                smooth.spline(AgeGradient,##gradient,
                              abundance,
                              control.spar = list(low = 0.2),
                              penalty = 1.4,
                              df = aber.pc2$complexity[names(taxaWant)[i]]))
    with(aber.pc2, plot(interpCurve(lambda[tag]), #lambda[tag],
                        taxaWant[,i],
                        main = colnames(taxaWant)[i],
                        ylim = ylim, ylab = NA, xlab = NA,
                        type = "n"))
    fit <- pred <- predict(mod, p[,"AgeGradient"])$y##[,1]
    with(aber.pc2, points(interpCurve(lambda[tag]), #lambda[tag],
                                      taxaWant[,i]))
    with(p, lines(AgeGradient, fit, lwd = 2, col = "black"))
}
title(xlab = "Distance along the Principle Curve", outer = TRUE,
      cex.lab = 1.3, line = 1)
title(ylab = "Proportional Abundance", outer = TRUE,
      cex.lab = 1.3, line = 2)
layout(1)
par(op)

## Or you can plot the fitted response curves in terms of the principal
## curve
want2 <- chooseTaxa(abernethy2, max.abun = 25, n.occ = 10, value = FALSE)
aber.curves <- sppResponse(aber.pc2)

layout(matrix(1:9, ncol = 3))
plot(aber.curves, which = want2)
layout(1)

## Alternative to Figure 9.24 - right panel.
## Here after placing each extracted gradient on a [0,1] scale,
## I scale each methods curve by the total variance in the data that
## that method's first axis explains.
## Hence, whilst the shapes of the curves for CA and PrCurves are similar
## the PrCurve is clearly doing a better job.
Depth <- abernethy$Depth
Age <- abernethy$Age
layout(matrix(c(1,2,2), ncol = 3))
nseg <- nrow(abernethy2) - 1
RoC <- unclass(abs(diff(gradientDist(aber.pc2)))) / diff(Age)
plot(y = Age[-length(Age)],
     x = RoC * 1000, type = "n",
     ylim = rev(range(Age)),
     ylab = "Age (calendar years BP)",
     xlab  = expression("Rate of Change" ~ (kyr^{-1})),
     main = "a", cex.main = 2,
     cex.lab = 1.5, cex.axis = 1.3)
segments(x0 = rep(0, nseg), y0 = Age[-length(Age)],
         x1 = RoC * 1000, y1 = Age[-length(Age)], lwd = 1)
plot((1 - gradientDist(aber.pc2)) * with(aber.pc, (totalDist - dist) / totalDist),
     orderBy = Age,
     xlim = rev(range(Age)),
     type = "o", flipAxes = TRUE, xlab = "Age (calendar years BP)",
     main = "b",
     cex = 1.5, pch = 21, col = "black", bg = "black", cex.main = 2, lwd = 1,
     cex.lab = 1.5, cex.axis = 1.3)
lines(gradientDist(aber.pca) * eigenvals(aber.pca)[1] / aber.pca$tot.chi, orderBy = Age,
      lty = "solid", flipAxes = TRUE, type = "o",
      cex = 1.5, pch = 22, col = "red", bg = "red", lwd = 1)
lines((1 - gradientDist(aber.ca)) * eigenvals(aber.pca)[1] / aber.pca$tot.chi,
      orderBy = Age,
      lty = "solid", flipAxes = TRUE, type = "o",
      cex = 1.5, pch = 23, col = "forestgreen", bg = "forestgreen", lwd = 1)
legend("topright", bty = "n", pch = 21:23,
       legend = c(expression(PCurve), expression(PCA[1]), expression(CA[1])),
       lty = c("solid","solid","solid"), lwd = 1, seg.len = 4,
       col = c("black","red","forestgreen"),
       pt.bg = c("black","red","forestgreen"), inset = 0.01, pt.cex = 2, cex = 1.5)
layout(1)
