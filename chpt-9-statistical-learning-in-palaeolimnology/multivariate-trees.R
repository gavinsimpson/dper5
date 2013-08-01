## Example of multivariate trees as constrained clustering
## uses abernethy forest data from analogue package

## load analogue and mvpart
require("mvpart")
require("analogue")

## source chpt 9 Functions
source("dper-chapter9-functions.R")

## load the Abernethy data
data(abernethy)

## Pull out the just the species data
aber.spp <- data.matrix(abernethy[, -(37:38)])

## Fit the mvpart
set.seed(5)
aber.mrt <- mvpart(aber.spp ~ Age, data = abernethy, xv = "none",
                   plot.add = FALSE)

## plot the cp table
cpPlot(aber.mrt, legend.loc = "bottomleft")

## prune the MRT to 1se - 6 nodes == 5 splits
aber.pmrt <- prune(aber.mrt, cp = 0.032)

## Plot the fitted tree
plot(aber.pmrt)
mytext.rpart(aber.pmrt, cex = 0.8, legend = TRUE, digits = 5)
## reset plotting par not being turned off
par(xpd = FALSE)

## MRT biplot - do my own version as this is not very good!
## rpart.pca(aber.pmrt, wgt.ave = TRUE, cex = 0.8)
myrpart.pca(aber.pmrt, wgt.ave = TRUE, cex = 1)
