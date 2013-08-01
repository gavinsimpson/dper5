## Custom function for colours
ggHueColours <- function(n, h = c(0, 360) + 15, l = 65, c = 100,
                         direction = 1, h.start = 0) {
    turn <- function(x, h.start, direction) {
        (x + h.start) %% 360 * direction
    }

    if ((diff(h) %% 360) < 1) {
      h[2] <- h[2] - 360 / n
    }

    hcl(h = turn(seq(h[1], h[2], length = n), h.start = h.start,
        direction = direction), c = c, l =  l)
}

## force predict to return class labels only
`pred.crt` <- function(object, newdata)
    predict(object, newdata = newdata, type = "class")


`cpPlot` <- function(RT, cex = 0.8, legend.loc = "topright") {
    tab <- RT$cptable
    xstd <- tab[, 5L]
    xerror <- tab[, 4L]
    relerror <- tab[, 3L]
    nsplit <- tab[, 2L]
    ns <- seq_along(nsplit)
    cp0 <- tab[, 1L]
    cp <- sqrt(cp0 * c(Inf, cp0[-length(cp0)]))
    ylim <- range(min(xerror - xstd), max(xerror + xstd),
                  max(relerror), min(relerror))
    plot(ns, relerror, type = "n", cex = cex, axes = FALSE,
         xlab = "", ylab = "Relative Error", ylim = ylim)
    title(xlab = "CP", line = 4)
    minpos <- which.min(xerror)
    se1 <- min(which(xerror <= (xerror[minpos] + xstd[minpos])))
    abline(h = (xerror)[minpos], lty = "dotted",
           col = "black")
    abline(v = ns[se1], lty = "dotted", col = "black")
    abline(v = ns[minpos], lty = "dotted", col = "black")
    points(ns, relerror, type = "o", cex = cex, pch = 21, bg = "white")
    points(ns, xerror, type = "o", cex = cex, pch = 21, bg = "black")
    segments(ns, xerror - xstd, ns, xerror + xstd)
    axis(2L, cex.axis = 0.8)
    axis(1L, at = ns, lab = as.character(signif(cp, 2L)), las = 2,
         cex.axis = 0.8)
    axis(3L, at = ns, lab = as.character(nsplit + 1), las = 2,
         cex.axis = 0.8)
    mtext("Size of tree", side = 3, line = 3)
    box()
    legend(legend.loc, legend = c("Apparent", expression(italic(k)~fold~CV)),
           pch = 21, pt.bg = c("white","black"), lty = "solid", bty = "n")
}

`figCaption` <- function(label, where = c("bottomright", "bottom",
                                "bottomleft","left","topleft","top",
                                "topright","right","center"),
                         inset = 0.05, cex = 1, ...) {
    inset <- rep(inset, 2)
    usr <- par("usr")
    insetx <- inset[1L] * (usr[2L] - usr[1L])
    insety <- inset[2L] * (usr[4L] - usr[3L])
    cin <- par("cin")
    Cex <- cex * par("cex")
    text.width <- max(abs(strwidth(label, units = "user", cex = cex)))
    xc <- Cex * xinch(cin[1L], warn.log = FALSE)
    yc <- Cex * yinch(cin[2L], warn.log = FALSE)
    xchar <- xc
    xextra <- 0
    yextra <- 0
    ymax <- yc * max(1, strheight(label, units = "user", cex = cex)/yc)
    ychar <- yextra + ymax
    h <- ychar + yc
    w0 <- text.width + (1 + 1) * xchar
    w <- w0 + 0.5 * xchar
    left <- switch(where,
                   bottomright = ,
                   topright = ,
                   right = usr[2L] - w - insetx,
                   bottomleft = ,
                   left = ,
                   topleft = usr[1L] + insetx,
                   bottom = ,
                   top = , center = (usr[1L] + usr[2L] - w)/2)
    top <- switch(where,
                  bottomright = ,
                  bottom = ,
                  bottomleft = usr[3L] + h + insety,
                  topleft = ,
                  top = ,
                  topright = usr[4L] - insety,
                  left = ,
                  right = ,
                  center = (usr[3L] + usr[4L] + h)/2)
    text(left, top, label, cex = cex, ...)
}

##
`coefPlot` <- function(mod, xvar = c("norm", "lambda", "dev"),
                       label = FALSE, varNames, norm,
                       ylab = "Coefficients",
                       hoff = 1, air = 1.2,
                       allDF = FALSE, ...) {
    beta <- mod$beta
    lambda <- mod$lambda
    dev <- mod$dev.ratio
    which <- nonzeroCoef(beta)
    beta <- as.matrix(beta[which, ])
    xvar <- match.arg(xvar)
    switch(xvar, norm = {
        index <- if (missing(norm)) apply(abs(beta), 2, sum) else norm
        iname <- "L1 Norm"
    }, lambda <- {
        index <- log(lambda)
        iname <- "Log Lambda"
    }, dev = {
        index <- dev
        iname <- "Fraction Deviance Explained"
    })
    dotlist <- list(...)
    type <- dotlist$type
    if (is.null(type))
        matplot(index, t(beta), xlab = iname, ylab = ylab,
                type = "l", ...)
    else matplot(index, t(beta), xlab = xlab, ylab = ylab,
                 ...)
    df <- mod$df
    if(allDF) {
        prettyVals <- approx(x = df, y = index,
                             xout = seq_len(mod$dim[1]), rule = 2)
        atdf <- prettyVals$y
        prettydf <- prettyVals$x
    } else {
        atdf <- pretty(index)
        prettydf <- trunc(approx(x = index, y = df, xout = atdf, rule = 2)$y)
    }
    axis(3, at = atdf, label = prettydf)
    mtext("DF", side = 3, line = 2.5, cex = par("cex"))
    if (label) {
        stopifnot(require(vegan))
        if(missing(varNames))
            varNames <- mod$beta@Dimnames[[1]]
        varNames <- varNames[which]
        linestack(beta[, ncol(beta)], varNames, add = TRUE,
                  at = par("usr")[2], hoff = hoff, air = air)
        rug(beta[, ncol(beta)], side = 4, ticksize = 0.02, lwd = 1)
    }
}

`cvPlot` <- function (cvobj, sign.lambda = 1, ...) {
    xlab <- "log(Lambda)"
    if (sign.lambda < 0)
        xlab <- paste("-", xlab, sep = "")
    plot.args <- list(x = sign.lambda * log(cvobj$lambda), y = cvobj$cvm,
                      ylim = range(cvobj$cvup, cvobj$cvlo), xlab = xlab,
                      ylab = cvobj$name, type = "n")
    new.args <- list(...)
    if (length(new.args))
        plot.args[names(new.args)] <- new.args
    do.call("plot", plot.args)
    glmnet:::error.bars(sign.lambda * log(cvobj$lambda), cvobj$cvup,
                        cvobj$cvlo, width = 0.01, col = "darkgrey")
    abline(v = sign.lambda * log(cvobj$lambda.min), lty = 2)
    abline(v = sign.lambda * log(cvobj$lambda.1se), lty = 2)
    points(sign.lambda * log(cvobj$lambda), cvobj$cvm, pch = 20,
        col = "black")
    axis(side = 3, at = sign.lambda * log(cvobj$lambda),
         labels = paste(cvobj$nz), tick = FALSE, line = 0)
    invisible()
}

`plot.kohchanges` <- function (x, main, keepMargins, ...)  {
    if (is.null(main))
        main <- "Training progress"
    nmaps <- ncol(x$changes)
    if (nmaps > 1) {
        if (!is.null(colnames(x$changes))) {
            varnames <- colnames(x$changes)
        }
        else {
            varnames <- paste("Matrix", 1:ncol(x$changes))
        }
    }
    if (nmaps == 2) {
        if (!keepMargins) {
            opar <- par("mar")
            on.exit(par(mar = opar))
        }
        par(mar = c(5.1, 4.1, 4.1, 4.1))
        huhn <- x$changes
        huhn[, 2] <- max(x$changes[, 1]) * huhn[, 2]/max(x$changes[, 2])
        ticks <- pretty(x$changes[, 2], length(axTicks(2)))
    }
    else {
        huhn <- x$changes
    }
    matplot(huhn, type = "l", lty = c("solid","dashed"), main = main,
            ylab = "Mean distance to closest unit",
            xlab = "Iteration", ...)
    abline(h = 0, col = "gray")
    if (nmaps == 2)
        axis(4, col.axis = 1, at = ticks *
             max(x$changes[, 1])/max(x$changes[, 2]), labels = ticks)
    if (nmaps > 1)
        legend("topright", legend = varnames, lty = c("solid","dashed"),
               col = "black", bty = "n")
    invisible()
}

## New MRT plotting function
`mytext.rpart` <- function (x, splits = TRUE, which = 4, label = "yval",
                            FUN = text, all.leaves = FALSE, pretty = NULL,
                            digits = getOption("digits") - 2, tadj = 0.65,
                            stats = TRUE, use.n = FALSE, bars = TRUE,
                            legend = FALSE, xadj = 1, yadj = 1, bord = FALSE,
                            big.pts = FALSE, ...) {
    if (!inherits(x, "rpart"))
        stop("Not legitimate rpart")
    if (!is.null(x$frame$splits))
        x <- rpconvert(x)
    frame <- x$frame
    col <- names(frame)
    method <- x$method
    ylevels <- attr(x, "ylevels")
    if (!is.null(ylevels <- attr(x, "ylevels")))
        col <- c(col, ylevels)
    if (is.na(match(label, col)))
        stop("Label must be a column label of the frame component of the tree")
    cxy <- par("cxy")
    if (!is.null(srt <- list(...)$srt) && srt == 90)
        cxy <- rev(cxy)
    xy <- mvpart:::rpartco(x)
    node <- as.numeric(row.names(x$frame))
    is.left <- (node%%2 == 0)
    node.left <- node[is.left]
    parent <- match(node.left/2, node)
    bars <- bars & is.matrix(frame$yval2)
    text.adj <- ifelse(bars, yadj * diff(range(xy$y))/12, 0)
    if (splits) {
        left.child <- match(2 * node, node)
        right.child <- match(node * 2 + 1, node)
        rows <- labels(x, pretty = pretty, digits = digits)
        if (which == 1)
            FUN(xy$x, xy$y + tadj * cxy[2], rows[left.child], ...)
        else {
            if (which == 2 | which == 4)
                FUN(xy$x, xy$y + tadj * cxy[2], rows[left.child],  pos = 2, ...)
            if (which == 3 | which == 4)
                FUN(xy$x, xy$y + tadj * cxy[2], rows[right.child], pos = 4, ...)
        }
    }
    leaves <- if (all.leaves)
        rep(TRUE, nrow(frame))
    else frame$var == "<leaf>"
    if (stats) {
        if (is.null(frame$yval2))
            stat <- x$functions$text(yval = frame$yval[leaves],
                dev = frame$dev[leaves], wt = frame$wt[leaves],
                ylevel = ylevels, digits = digits, n = frame$n[leaves],
                use.n = use.n)
        else stat <- x$functions$text(yval = frame$yval2[leaves, ],
                                      dev = frame$dev[leaves],
                                      wt = frame$wt[leaves],
                                      ylevel = ylevels, digits = digits,
                                      n = frame$n[leaves], use.n = use.n)
        FUN(xy$x[leaves], xy$y[leaves] - tadj * cxy[2] - text.adj,
            stat, adj = 0.5, xpd = NA, ...)
    }
    if (bars) {
        ## Use a grey palette
        Gray <- colorRampPalette(c("white", "grey20"))
        bar.vals <- x$functions$bar(yval2 = frame$yval2)
        mvpart:::sub.barplot(xy$x, xy$y, bar.vals, leaves, xadj = xadj,
                             yadj = yadj, bord = bord, line = TRUE,
                             ##col = c("lightblue", "blue", "darkblue"))
                             col = Gray(5))
        rx <- range(xy$x)
        ry <- range(xy$y)
        if (!is.null(ylevels))
            bar.labs <- ylevels
        else bar.labs <- dimnames(x$y)[[2]]
        if (legend & !is.null(bar.labs))
            legend(min(xy$x) - 0.1 * rx, max(xy$y) + 0.05 * ry,
                   ##bar.labs, col = c("lightblue", "blue", "darkblue"),
                   bar.labs, bty = "n",
                   border = "black", fill = Gray(5), ...)
    }
    if (big.pts)
        points(xy$x[leaves], xy$y[leaves], pch = 16, cex = 3 *
            par()$cex, col = 2:(sum(leaves) + 1))
    invisible()
}

ggHueColours <- function(n, h = c(0, 360) + 15, l = 65, c = 100,
                         direction = 1, h.start = 0) {
    turn <- function(x, h.start, direction) {
        (x + h.start) %% 360 * direction
    }

    if ((diff(h) %% 360) < 1) {
      h[2] <- h[2] - 360 / n
    }

    hcl(h = turn(seq(h[1], h[2], length = n), h.start = h.start,
        direction = direction), c = c, l =  l)
}

`myrpart.pca` <-  function (tree, pts = TRUE, plt.allx = TRUE, speclabs = TRUE,
                            specvecs = TRUE, wgt.ave = FALSE, add.tree = TRUE,
                            cv1 = 1, cv2 = 2, chulls = TRUE,
                            ...) {
    stopifnot(require(analogue))
    if (tree$method != "mrt")
        stop("Only for multivariate trees !! \n")
    if (nrow(tree$frame) < 4)
        stop("Only 2 terminal nodes -- PCA not done !! \n")
    old.par <- par(mar = c(5,4,2,2) + 0.1)
    on.exit(par(old.par))
    on.exit(layout(1), add = TRUE)
    frame <- tree$frame
    ncf <- ncol(frame)
    data <- tree$y
    ny <- ncol(data)
    treegrps <- tree$where
    specs <- dimnames(data)[[2]]
    leaves <- frame$var == "<leaf>"
    n <- length(leaves)
    ln <- sum(leaves)
    lnot <- sum(!leaves)
    key <- dimnames(frame)[[1]]
    node <- as.numeric(key)
    even.node <- node[even <- node%%2 == 0]
    num <- length(specs)
    node.means <- as.matrix(frame[, ncf])
    tnode.means <- node.means[leaves, ]
    dimnames(node.means) <- list(key, specs)
    mat <- amat <- node.means - node.means[rep(1, n), ]
    mat <- mat[leaves, ]
    temp <- mat[rep(1:ln, frame[leaves, 2]), ]
    z <- svd(temp)
    maxd <- sum(z$d > 1e-06)
    d <- diag(z$d[1:maxd])
    xall <- z$u[, 1:maxd, drop = FALSE] %*% d
    x <- amat %*% (z$v)[, 1:maxd, drop = FALSE]
    xlv <- x[leaves, ]
    if (!wgt.ave)
        y <- z$v[, 1:maxd, drop = FALSE]
    else {
        specvecs <- FALSE
        rc <- apply(tnode.means * frame$n[leaves], 2, sum)
        wgt <- diag(1/rc) %*% t(tnode.means * frame$n[leaves])
        y <- wgt %*% xlv
    }
    label <- 4:(3 + num)
    dstat <- signif(frame[leaves, "yval2"], digits = options()$digits)
    ln <- dim(dstat)[1]
    stat <- vector("character", length = ln)
    for (i in 1:ln) stat[i] <- paste(dstat[i, ], collapse = ", ")
    ymax <- max(dstat)
    ymin <- min(0, min(dstat))
    treegrps <- as.numeric(factor(treegrps))
    xx <- (scale(as.matrix(data), center = TRUE, scale = FALSE) %*%
        z$v)[, 1:maxd, drop = FALSE]
    xrb <- rbind(x, xx)
    if (plt.allx) {
        mxx <- sqrt(apply(xrb[, c(cv1, cv2)]^2, 1, sum))
    }
    else mxx <- sqrt(apply(x[, c(cv1, cv2)]^2, 1, sum))
    ## cvar <- round((100 * z$d[1:maxd]^2)/sum(z$d[1:maxd]^2), digits = 2)
    ## cvar2 <- round(diag(cor(xall, xx[order(tree$where), ]))[1:maxd], 3)
    ## dlabs <- paste("   Dim ", c(1:maxd), " ", cvar, "% : [", cvar2, "]")
    myy <- sqrt(apply(y[, c(cv1, cv2)]^2, 1, sum))
    sc <- ifelse(wgt.ave, 1, max(mxx)/max(myy))
    layout(matrix(1:2, ncol = 2))
    plot(c(sc * y[, cv1], xx[, cv1]), c(sc * y[, cv2], xx[, cv2]),
         axes = FALSE, xlab = "PC1", ylab = "PC2", type = "n", asp = 1)
    abline(h = 0, col = "darkgrey", lty = "dotted")
    abline(v = 0, col = "darkgrey", lty = "dotted")
    cxy <- par("cxy")
    sze <- par()$fin/par()$din
    adj <- ifelse(pts, cxy[2] * sze[2], 0)
    if (specvecs)
        segments(sc * y[, cv1], sc * y[, cv2], rep(0, nrow(y)),
                 rep(0, nrow(y)), col = "gray", lty = 1)
    ## mtext(dlabs[cv1], side = 1, las = 0, adj = 0, line = 0,
    ##       cex = 0.85 * par()$cex)
    ## mtext(dlabs[cv2], side = 2, las = 0, adj = 0, line = 0,
    ##       cex = 0.85 * par()$cex)
    greyPal <- colorRampPalette(c(grey(0.8), grey(0.2)))
    ##Cols <- greyPal(6)
    ##Cols <- 1:6
    Cols <- ggHueColours(6)
    if (add.tree) {
        pp <- match(c(even.node, even.node + 1), node)
        nn <- length(even.node)
        from <- pp[1:nn]
        to <- pp[(nn + 1):(2 * nn)]
        segments(x[from, cv1], x[from, cv2], x[to, cv1], x[to, cv2])
    }
    if (chulls) {
        unitg <- sort(unique(treegrps))
        for (i in 1:length(unitg)) {
            hpts <- chull(xx[unitg[i] == treegrps, c(cv1, cv2)])
            hpts <- c(hpts, hpts[1])
            lines(xx[unitg[i] == treegrps, c(cv1, cv2)][hpts, ],
                  col = Cols[i])
        }
    }
    if (plt.allx) {
        unitg <- sort(unique(treegrps))
        for (i in 1:length(unitg))
            points(xx[unitg[i] == treegrps, cv1], xx[unitg[i] == treegrps, cv2],
                   pch = c(21:24,21:22)[i], col = 1, bg = Cols[i], cex = 0.6)
    }
    if (pts) {
        lvnode <- sort(node[leaves])
        for (i in 1:length(lvnode))
            points(xlv[, cv1][lvnode[i] == lvnode],
                   xlv[, cv2][lvnode[i] == lvnode],
                   pch = c(21:24,21:22)[i],
                   cex = 0.9, col = 1, bg = Cols[i], lwd = 1)
    }
    axis(1)
    axis(2)
    axis(3, labels = FALSE)
    axis(4, labels = FALSE)
    box()
    ## species plot
    plot(c(sc * y[, cv1], xx[, cv1]), c(sc * y[, cv2], xx[, cv2]),
         axes = FALSE, xlab = "PC1", ylab = "PC2", type = "n", asp = 1)
    abline(h = 0, col = "darkgrey", lty = "dotted")
    abline(v = 0, col = "darkgrey", lty = "dotted")
    if (speclabs) {
        ##browser()
        want <- names(chooseTaxa(abernethy[, -(37:38)], n.occ = 5, max.abun = 5,
                                 type = "AND"))
        WANT <- which(specs %in% want)
        specs2 <- make.cepnames(specs)
        text(sc * y[WANT, cv1], sc * y[WANT, cv2] + 0.5 * adj * specvecs *
             (y[WANT, cv2] > 0), specs2[WANT], col = "black", cex = 0.6)
        ## text(sc * y[, cv1], sc * y[, cv2] + 0.5 * adj * specvecs *
        ##      (y[, cv2] > 0), specs[WANT], col = "black", cex = 0.7)
    }
    axis(1)
    axis(2)
    axis(3, labels = FALSE)
    axis(4, labels = FALSE)
    box()
    invisible(list(y = sc * y, xlv = xlv, xx = xx))
}
