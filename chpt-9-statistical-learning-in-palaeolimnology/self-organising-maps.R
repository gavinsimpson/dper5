#######################################################################
##                                                                   ##
## Script to analyse the SWAP diatom-pH training set via SOMs        ##
##                                                                   ##
## Gavin L. Simpson                                                  ##
##                                                                   ##
## Version 0.9-0                                                     ##
##                                                                   ##
## History:                                                          ##
##                                                                   ##
## 28-Apr-2011 -  0.9-0 - Script to accompany draft chapter          ##
##                                                                   ##
#######################################################################

## load packages
require("kohonen")  ## package to fits SOMs
require("vegan")    ## reading cep files
require("analogue") ## selecting taxa

## load chpater 8 function
source("dper-chapter9-functions.R")

## load the example data
##  + supplsp.cep contains the SWAP 138 and RLGH core diatom data
##  + swap138.log contains the SWAP 138 Log10 transformed chemistry
diatom <- read.cep("../data/supplsp.cep", force = TRUE, maxdata = 20000)
chem <- read.cep("../data/swap138.log", force = TRUE)
## better variable names for chem
names(chem) <- c("pH", "Conductivity", "TOC", "Ca", "Mg", "K",
                 "SO4", "Cl", "Alkalinity", "Total Aluminium")

## subset the diatom data into SWAP and RLGH sets
## rows 1-138 == SWAP
## rows 138-  == RLGH
swap <- diatom[1:138, ]
swap <- chooseTaxa(swap, n.occ = 5, max.abun = 2)
rlgh <- diatom[139:239, ]

## SOM of the chemistry data ##########################################
chem.sc <- scale(chem) ## scale to 0 mean and unit variance
set.seed(42)           ## uses random starts
## fit SOM without Conductivity
chem.som <- som(chem.sc[, -2], grid = somgrid(6,5, "hexagonal"))

## Plot the training progress to check for convergence
plot(chem.som, type = "changes", main = "SOM: SWAP Lake-water Chemistry")

## Plot the codebook vectors for the chemistry SOM
plot(chem.som, main = "SWAP Lake-water Chemistry",
     palette.name = ggHueColours)
Gray <- colorRampPalette(c("white", "grey20"))
plot(chem.som, main = "SWAP Lake-water Chemistry",
     palette.name = Gray)

## Plot the SOM grid showing the mapping quality
plot(chem.som, type = "quality", main = "Mean distance to unit",
     palette.name = Gray)

## Plot the number of observations in each unit
plot(chem.som, type = "counts", palette.name = Gray)
scatter <- matrix(rnorm(nrow(chem) * 2, sd = 0.1), ncol = 2)
points(chem.som$grid$pts[chem.som$unit.classif, ] + scatter,
       pch = 21, bg = "white", col = "black")
## End example ########################################################

## XYF SOM relating diatoms and chemistry #############################
set.seed(11)
swap.xyf <- xyf(data = chem.sc, Y = sqrt(data.matrix(swap)),
                xweight = 0.5,
                grid = somgrid(6, 5, "hexagonal"), rlen = 100)

## Plot the training progress to check for convergence
plot(swap.xyf, type = "changes", main = "XYF: SWAP Diatoms & Chemistry",
     col = "black")

## Plot the SOM grid showing the mapping quality
plot(swap.xyf, type = "quality", main = "Mean distance to unit",
     palette.name = Gray)

## Plot the number of observations in each unit
plot(swap.xyf, type = "counts", palette.name = Gray)
points(swap.xyf$grid$pts[swap.xyf$unit.classif, ] + scatter,
       pch = 21, bg = "white", col = "black")

## Plot the codebook vectors for the chemistry map X
plot(swap.xyf, main = "SWAP Diatoms & Chemistry", whatmap = 1,
     palette.name = ggHueColours)
plot(swap.xyf, main = "SWAP Diatoms & Chemistry", whatmap = 1,
     palette.name = Gray)

## predict for some key diatom taxa
diatom.pred <- predict(swap.xyf)

layout(matrix(1:4, ncol = 2, byrow = TRUE))
plot(swap.xyf, type = "property", palette.name = Gray,
     property = diatom.pred$unit.predictions[, "AC013A"]^2,
     main = "Achnanthese minutissima")
plot(swap.xyf, type = "property", palette.name = Gray,
     property = diatom.pred$unit.predictions[, "BR006A"]^2,
     main = "Brachysira brebissonii")
plot(swap.xyf, type = "property", palette.name = Gray,
     property = diatom.pred$unit.predictions[, "TA003A"]^2,
     main = "Tabellaria binalis")
plot(swap.xyf, type = "property", palette.name = Gray,
     property = diatom.pred$unit.predictions[, "EU047A"]^2,
     main = "Eunotia incisa")
layout(1)
## End example ########################################################

## Calibration example ################################################
set.seed(11)
pH.xyf <- xyf(data = sqrt(data.matrix(swap)), Y = chem$pH, xweight = 0.5,
              grid = somgrid(6, 5, "hexagonal"), rlen = 100)

## Plot the training progress to check for convergence
plot(pH.xyf, type = "changes", main = "XYF: SWAP pH", col = "black")

## Plot the codebook vectors for the pH map Y
pH.pred <- predict(pH.xyf)
plot(pH.xyf, type = "property", palette.name = Gray,
     property = pH.pred$unit.predictions, main = "Predicted pH")

## RMSE
(rmsep <- sqrt(mean((chem$pH - pH.pred$prediction)^2)))

## Predict for RLGH
## select species we want
rlgh <- rlgh[, colnames(swap)]
## read in the ages
age <- read.table("rlgh_age.txt", row.names = 1)[, 1]
rlgh.pred <- predict(pH.xyf, newdata = sqrt(data.matrix(rlgh)))

## plot
layout(matrix(c(1,2), ncol = 1))
op <- par(mar = c(5.5,4,4,2) + 0.1)
plot(age, rlgh.pred$prediction, xlim = rev(range(age)), type = "b",
     main = "a", ylab = "Predicted pH", xlab = "Calendar Years BP",
     cex = 0.7, las = 1)
par(op)
rlgh.map <- map(pH.xyf, newdata = sqrt(data.matrix(rlgh)), whatmap = 2)
plot(pH.xyf, type = "property", property = pH.pred$unit.predictions,
     keepMargins = FALSE, palette.name = Gray, main = "b")
set.seed(13)
scatter <- matrix(rnorm(nrow(rlgh) * 2, sd = 0.11), ncol = 2)
points(pH.xyf$grid$pts[rlgh.map$unit.classif, ] + scatter,
       pch = 21, bg = "white", col = "black", cex = 0.8, xpd = NA)
layout(1)

## End example ########################################################

## Figures for chapter ################################################

## BW palette
Gray <- colorRampPalette(c("white", "grey20"))

## Figure 8.X
layout(matrix(1:4, ncol = 2, byrow = TRUE))
## Plot the training progress to check for convergence
plot(chem.som, type = "changes", main = "a")
## Plot the codebook vectors for the chemistry SOM
plot(chem.som, main = "b", palette.name = Gray)
## Plot the SOM grid showing the
Gray <- colorRampPalette(c("white", "black"))
plot(chem.som, type = "quality", main = "c", palette.name = Gray)
## Plot the number of observations in each unit
plot(chem.som, type = "counts", main = "d", palette.name = Gray)
set.seed(42)
scatter <- matrix(rnorm(nrow(chem) * 2, sd = 0.13), ncol = 2)
points(chem.som$grid$pts[chem.som$unit.classif, ] + scatter,
       pch = 21, bg = "white", col = "black", cex = 0.8)
layout(1)

## Figure 8.X
layout(matrix(1:2, ncol = 2, byrow = TRUE))
op <- par(las = 1, cex.axis = 0.8)
## Plot the training progress to check for convergence
plot(swap.xyf, type = "changes", main = "a", col = "black")
plot(swap.xyf, main = "b", whatmap = 1,
     palette.name = Gray)
par(op)
layout(1)

## Figure 8.X
layout(matrix(1:4, ncol = 2, byrow = TRUE))
plot(swap.xyf, type = "property", palette.name = Gray,
     property = diatom.pred$unit.predictions[, "AC013A"]^2,
     main = "Achnanthes minutissima")
plot(swap.xyf, type = "property", palette.name = Gray,
     property = diatom.pred$unit.predictions[, "BR006A"]^2,
     main = "Brachysira brebissonii")
plot(swap.xyf, type = "property", palette.name = Gray,
     property = diatom.pred$unit.predictions[, "TA003A"]^2,
     main = "Tabellaria binalis")
plot(swap.xyf, type = "property", palette.name = Gray,
     property = diatom.pred$unit.predictions[, "EU047A"]^2,
     main = "Eunotia incisa")
layout(1)

## Figure 8.X
layout(matrix(c(1,2), ncol = 1))
op <- par(mar = c(5.5,4,4,2) + 0.1)
plot(age, rlgh.pred$prediction, xlim = rev(range(age)), type = "b",
     main = "a", ylab = "Predicted pH", xlab = "Radiocarbon Years BP",
     cex = 0.7, las = 1)
par(op)
rlgh.map <- map(pH.xyf, newdata = sqrt(data.matrix(rlgh)), whatmap = 2)
plot(pH.xyf, type = "property", property = pH.pred$unit.predictions,
     keepMargins = FALSE, palette.name = Gray, main = "b")
set.seed(13)
scatter <- matrix(rnorm(nrow(rlgh) * 2, sd = 0.13), ncol = 2)
points(pH.xyf$grid$pts[rlgh.pred$unit.classif, ] + scatter,
       pch = 21, bg = "white", col = "black", cex = 0.8, xpd = NA)
layout(1)


## Colour Figures for chapter ################################################

## colour palette
## ggHueColours() from above

## Figure 8.X
## Plot the training progress to check for convergence
plot(chem.som, type = "changes", main = "a")
## Plot the codebook vectors for the chemistry SOM
plot(chem.som, main = "b", palette.name = ggHueColours)
## Plot the SOM grid showing the
plot(chem.som, type = "quality", main = "c", palette.name = Gray)
## Plot the number of observations in each unit
plot(chem.som, type = "counts", main = "d", palette.name = Gray)
set.seed(42)
scatter <- matrix(rnorm(nrow(chem) * 2, sd = 0.13), ncol = 2)
points(chem.som$grid$pts[chem.som$unit.classif, ] + scatter,
       pch = 21, bg = "white", col = "black", cex = 0.8)
layout(1)

## Figure 8.X
layout(matrix(1:2, ncol = 2, byrow = TRUE))
op <- par(las = 1, cex.axis = 0.8)
## Plot the training progress to check for convergence
plot(swap.xyf, type = "changes", main = "a")
plot(swap.xyf, main = "b", whatmap = 1,
     palette.name = ggHueColours)
par(op)
layout(1)

## Figure 8.X
layout(matrix(1:4, ncol = 2, byrow = TRUE))
plot(swap.xyf, type = "property", palette.name = Gray,
     property = diatom.pred$unit.predictions[, "AC013A"]^2,
     main = "Achnanthes minutissima")
plot(swap.xyf, type = "property", palette.name = Gray,
     property = diatom.pred$unit.predictions[, "BR006A"]^2,
     main = "Brachysira brebissonii")
plot(swap.xyf, type = "property", palette.name = Gray,
     property = diatom.pred$unit.predictions[, "TA003A"]^2,
     main = "Tabellaria binalis")
plot(swap.xyf, type = "property", palette.name = Gray,
     property = diatom.pred$unit.predictions[, "EU047A"]^2,
     main = "Eunotia incisa")
layout(1)

## Figure 8.X
layout(matrix(c(1,2), ncol = 1))
op <- par(mar = c(5.5,4,4,2) + 0.1)
plot(age, rlgh.pred$prediction, xlim = rev(range(age)), type = "b",
     main = "a", ylab = "Predicted pH", xlab = "Radiocarbon Years BP",
     cex = 0.7, las = 1)
par(op)
rlgh.map <- map(pH.xyf, newdata = sqrt(data.matrix(rlgh)), whatmap = 2)
plot(pH.xyf, type = "property", property = pH.pred$unit.predictions,
     keepMargins = FALSE, palette.name = Gray, main = "b")
set.seed(13)
scatter <- matrix(rnorm(nrow(rlgh) * 2, sd = 0.13), ncol = 2)
points(pH.xyf$grid$pts[rlgh.pred$unit.classif, ] + scatter,
       pch = 21, bg = "white", col = "black", cex = 0.8, xpd = NA)
layout(1)

