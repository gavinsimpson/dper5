## Data Analysis Script for DPER Vol. 5
## Chapter 9
##
## Gavin L. Simpson
##
## Version 0.1 10 Dec 2010

## Univariate Trees
require(rpart)

## source Chpt9 Functions
source("dper-chapter9_functions.R")

## Load the SCP data ~/work/data/scp/SCP5.csv
scp <- read.csv("~/work/data/scp/SCP5.csv")
## process the names:
names(scp) <- c("Site","Particle","Fuel","FuelID","Fuel2","Fuel2ID","Index",
                "Na","Mg","Al","Si","P","S","Cl","Cd","K","Ca","Ti","V",
                "Ba","Cr","Mn","Fe","Co","Ni","Cu","Zn","Pb")
## fix up some levels on some variables
scp <- within(scp, levels(Fuel) <- c("Black Coal","Brown Coal","Coal","Oil",
                                     "Peat","Oil Shale"))
scp <- within(scp, levels(Fuel2) <- c("Coal","Oil","Oil Shale"))
## set rownames
rownames(scp) <- with(scp, as.character(Particle))
## drop some variables we won't use here
scp <- scp[, -c(1,2,4,6,7)]

## fit a classification tree:
set.seed(1)
scp2.rt <- rpart(Fuel2 ~ ., data = scp[, -1],
                control = rpart.control(cp = 0.0001))

scp2.prt <- prune(scp2.rt, cp = 0.0022)

op <- par(mar = rep(0,4) + 0.1)
plot(scp2.prt, branch = 0.8, uniform = FALSE)
text(scp2.prt, cex = 0.8)
par(op)

cpPlot(scp2.rt)

## Confusion matrix - apparent ability!
confus.scp.ct <- table(predict(scp2.rt, type = "class"), with(scp, Fuel2))
## correctly classified
diag(confus.scp.ct)
## as proportion
sum(diag(confus.scp.ct)) / sum(confus.scp.ct)
## overall apparent error rate
1 - sum(diag(confus.scp.ct)) / sum(confus.scp.ct)
## individual error rates
1 - (diag(confus.scp.ct) / rowSums(confus.scp.ct))

## add individual error rates to the table:
class(confus.scp.ct) <- "matrix"
confus.scp.ct <- data.frame(confus.scp.ct,
                            class.error = 1 - (diag(confus.scp.ct) /
                            rowSums(confus.scp.ct)),
                            check.names = FALSE)

## CV and bagging
require(ipred)

## force predict to return class labels only
`pred.crt` <- function(object, newdata)
    predict(object, newdata = newdata, type = "class")

## 10-fold cv of classification tree for 3-Fuel SCP data
set.seed(13)
errorest(Fuel2 ~ ., data = scp[,-1], model = rpart,
         estimator = "cv", predict = pred.crt)

set.seed(1)
scp.bag <- bagging(Fuel2 ~ ., data = scp[, -1], coob = TRUE, nbag = 500)

## The OOB error rate:
scp.bag$err

## Confusion matrix
confus.scp.bag <- table(predict(scp.bag), with(scp, Fuel2))
## correctly classified
diag(confus.scp.bag)
## as proportion
sum(diag(confus.scp.bag)) / sum(confus.scp.bag)
## overall apparent error rate
1 - sum(diag(confus.scp.bag)) / sum(confus.scp.bag)
## individual error rates
1 - (diag(confus.scp.bag) / rowSums(confus.scp.bag))

## add individual error rates to the table:
class(confus.scp.bag) <- "matrix"
confus.scp.bag <- data.frame(confus.scp.bag,
                             class.error = 1 - (diag(confus.scp.bag) /
                             rowSums(confus.scp.bag)),
                             check.names = FALSE)

## Random forests for 3-fuel scp data
require(randomForest)
set.seed(1)
scpRF <- randomForest(Fuel2 ~ ., data = scp[, -1], proximity = TRUE,
                      importance = TRUE, ntree = 500)

legcols <- c("black","red","orange","blue")
##legcols <- rep("black", 4)
leglwds <- c(2,1,1,1)
##legltys <- c("solid","dashed","dotdash","dotted")
legltys <- rep("solid", 4)
out <- plot(scpRF, col = legcols, lwd = leglwds, lty = legltys)
legend("topright",
       legend = c("OOB", with(scp, levels(Fuel2))),
       col = legcols, lwd = leglwds, bty = "n", lty = legltys, seg.len = 4)

legcols <- c("black","red","orange","blue")
leglwds <- c(1,1,1,1)
legltys <- rep("solid", 4)
matplot(out, type = "l", xlab = "Number of trees", ylab = "Error rate",
        lwd = leglwds, col = legcols, lty = legltys)
legend("topright",
       legend = c("OOB", with(scp, levels(Fuel2))),
       col = legcols, lwd = leglwds, lty = legltys, bty = "n",
       seg.len = 4)

legcols <- rep("black", 4)
leglwds <- c(1,1,1,1)
legltys <- c("solid","dashed","dotdash","dotted")
matplot(out, type = "l", xlab = "Number of trees", ylab = "Error rate",
        lwd = leglwds, col = legcols, lty = legltys)
legend("topright",
       legend = c("OOB", with(scp, levels(Fuel2))),
       col = legcols, lwd = leglwds, lty = legltys, bty = "n",
       seg.len = 4)

## Confusion matrix
(confus.scp.rf <- scpRF$confusion)

## Variable importance
impRF <- importance(scpRF)

## my own plot of variable importance...
layout(matrix(1:2, ncol = 2))
op <- par(mar = c(4,4,4,1), las = 2)
dotchart(sort(impRF[,5]), main = "Mean decrease in Gini")
dotchart(sort(impRF[,4]), main = "Mean decrease in Accuracy")
par(op)
layout(1)

## my own plot of variable importance...
## class-wise
layout(matrix(1:3, ncol = 3))
op <- par(mar = c(5,4,4,1), las = 2)
xlab <- "Mean decrease in accuracy"
dotchart(sort(impRF[,1]), main = "Coal", xlim = c(0,1.2), xlab = xlab)
dotchart(sort(impRF[,2]), main = "Oil", xlim = c(0,1.2), xlab = xlab)
dotchart(sort(impRF[,3]), main = "Oil Shale", xlim = c(0,1.2), xlab = xlab)
par(op)
layout(1)
