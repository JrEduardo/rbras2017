##======================================================================
## Tranformations
lambda2mu <- function(lambda, nu) {
    phi <- log(nu)
    mu <- lambda^(1/nu) - (nu - 1)/(2 * nu)
    list("mu" = mu, "phi" = phi)
}

mu2lambda <- function(mu, phi) {
    nu <- exp(phi)
    lambda <- (mu + (nu-1)/(2*nu))^nu
    list("lambda" = lambda, "nu" = nu)
}

##======================================================================
## Compute momentes mean and variance
calc_moments <- Vectorize(FUN = function(mu, phi) {
    pars <- mu2lambda(mu, phi)
    mu <- cmpreg::calc_mean_cmp(log(pars[["lambda"]]), phi)
    va <- cmpreg::calc_var_cmp(log(pars[["lambda"]]), phi)
    ## mu <- do.call(compoisson::com.mean, pars)
    ## va <- do.call(compoisson::com.var, pars)
    c("mean" = mu, "var" = va)
}, vectorize.args = c("mu", "phi"))

##======================================================================
## Negative log-likelihood functions

## COM-Poisson Model (Original parameterization)
llcmp <- function (params, y, X, offset = NULL, sumto = NULL) {
    if (!is.null(offset)) {
        stop("This model doesn't support offset yet.")
    }
    phi <- params[1]
    nu <- exp(phi)
    Xb <- X %*% params[-1]
    i <- 0:sumto
    zs <- sapply(Xb, function(loglambda) {
        sum(exp(i * loglambda - nu * lfactorial(i)))
    })
    Z <- sum(log(zs))
    ll <- sum(y * Xb - nu * lfactorial(y)) - Z
    return(-ll)
}

## COM-Poisson Model (mean-parameterization)
llcmp2 <- function (params, y, X, offset = NULL, sumto = NULL) {
    if (!is.null(offset)) {
        stop("This model doesn't support offset yet.")
    }
    phi2 <- params[1]
    nu <- exp(phi2)
    ##-------------------------------------------
    mu <- exp(X %*% params[-1])
    Xb <- nu * log(mu + (nu - 1)/(2 * nu))
    ##-------------------------------------------
    i <- 0:sumto
    zs <- sapply(Xb, function(loglambda) {
        sum(exp(i * loglambda - nu * lfactorial(i)))
    })
    Z <- sum(log(zs))
    ll <- sum(y * Xb - nu * lfactorial(y)) - Z
    return(-ll)
}

##======================================================================
## Framework for fit count models (using bbmle::mle2 function)
fitcm <- function(formula, data, start = NULL,
                  model = c("CP", "CP2"), ...) {
    dots <- list(...)
    model <- match.arg(model)
    ##-------------------------------------------
    frame <- model.frame(formula, data)
    terms <- attr(frame, "terms")
    y <- model.response(frame)
    X <- model.matrix(terms, frame)
    off <- model.offset(frame)
    if (is.null(start)) {
        m0 <- glm.fit(x = X, y = y, offset = off, family = poisson())
        start <- c(0, m0$coefficients)

    }
    ##-------------------------------------------
    if (model == "CP") {
        names(start)[1] <- "phi"
        bbmle::parnames(llcmp) <- names(start)
        sumto <- dots$sumto
        dots$sumto <- NULL
        args <- append(
            list(minuslog = llcmp, start = start,
                 data = list(y = y, X = X, offset = off,
                             sumto = sumto, model = model),
                 vecpar = TRUE), dots)
        fit <- suppressWarnings(
            do.call(what = bbmle::mle2, args = args))
    }
    if (model == "CP2") {
        names(start)[1] <- "phi2"
        bbmle::parnames(llcmp2) <- names(start)
        sumto <- dots$sumto
        dots$sumto <- NULL
        args <- append(
            list(minuslog = llcmp2, start = start,
                 data = list(y = y, X = X, offset = off,
                             sumto = sumto, model = model),
                 vecpar = TRUE), dots)
        fit <- suppressWarnings(
            do.call(what = bbmle::mle2, args = args))
    }
    slot(fit, "call.orig") <- match.call()
    return(fit)
}

##======================================================================
## Analysis of Deviance Table to glm and mle2 objetcs
getAnova <- function(..., print = TRUE) {
    model.list <- list(...)
    test <- "Chisq"
    cls <- sapply(model.list, function(x) class(x)[1])
    if (all(cls == "glm")) {
        fam <- sapply(model.list, function(x) x$family$family)
        if (all(fam == "quasipoisson")) {
            test <- "F"
        }
    }
    nps <- sapply(model.list, function(x) attr(logLik(x), "df"))
    aic <- sapply(model.list, function(x) AIC(x))
    if (test == "Chisq") {
        lls <- sapply(model.list, function(x) logLik(x))
        cst <- 2 * diff(lls)
        pvs <- pchisq(cst, df = diff(nps), lower.tail = FALSE)
        tab <- cbind(
            "np" = nps,
            "ll" = lls,
            "AIC" = aic,
            "2*dif" = c(NA, cst),
            "dif np" = c(NA, diff(nps)),
            "Pr(>Chisq)" = c(NA, pvs))
    }
    if (test == "F") {
        ## ver pág. 187 expresssão 7.10 em Modern Applied, Ripley
        sigma <- sapply(model.list, function(x) summary(x)$dispersion)
        df.sigma <- sapply(model.list, function(x) x$df.residual)
        dev <- sapply(model.list, function(x) deviance(x))
        dfs <- c(NA, -diff(dev))
        dnp <- c(NA, diff(nps))
        F <- c(NA, -diff(dev))/(dnp * sigma)
        pvs <- pf(F, df1 = dnp, df2 = df.sigma, lower.tail = FALSE)
        tab <- cbind(
            "np" = nps,
            "dev" = dev,
            "AIC" = aic,
            "F" = F,
            "dif np" = dnp,
            "Pr(>F)" = pvs)
    }
    rownames(tab) <- 1:length(model.list)
    if (print) printCoefmat(tab, na.print = "", cs.ind = 1)
    invisible(tab)
}

## Get the estimated coefficients of the parameters
getCoefs <- function(model, digits = 3, rownames = NULL) {
    ##-------------------------------------------
    glue <- function(est, d = digits) {
        est
        ## apply(est, 1, FUN = function(x)
        ##     paste0(format(x[1], digits = d, nsmall = d), " (",
        ##            format(x[2], digits = d, nsmall = d), ")"))
    }
    ##-------------------------------------------
    cls <- class(model)[1]
    est <- switch(
        cls,
        "glm" = {
            fam <- model$family$family
            if (!fam %in% c("poisson", "quasipoisson")) {
                stop("Família de distribuições não contemplada")
            }
            switch(
                fam,
                "poisson" = {
                    glue(rbind(c(NA, NA),
                               summary(model)$coef[, c(1, 3)]))
                },
                "quasipoisson" = {
                    glue(rbind(c(summary(model)$dispersion, NA),
                          summary(model)$coef[, c(1, 3)]))
                }
            )
        },
        "mle2" = {
            glue(bbmle::summary(model)@coef[, c(1, 3)])
        }
    )
    if (!is.null(rownames)) {
        rownames(est) <- rownames
    }
    colnames(est) <- c("Est", "Est/EP")
    return(est)
}

##======================================================================
## Predict Values

## Get the interval values of linear predictors (for Delta method with
## cholesky decomposition)
cholV_eta <- function(V, X, b, qn) {
    eta <- c(X %*% b)
    if (length(qn) > 1) {
        U <- chol(V)
        stderr <- sqrt(
            apply(X %*% t(U), MARGIN = 1, FUN = function(x) sum(x^2))
        )
        eta <- sweep(x = outer(stderr, qn, FUN = "*"),
                     MARGIN = 1,
                     STATS = eta,
                     FUN = "+")
    }
    return(eta)
}

## Get the means (inverse link function)
calc_mean <- function(eta, dispersion, model = c("CP", "CP2"),
                      tol = 1e-5, offset = NULL, ...) {
    if (is.null(offset)) {
        offset <- 1L
    }
    model <- match.arg(model)
    names(eta) <- NULL
    names(dispersion) <- NULL
    pars <- data.frame(eta = eta, dispersion = dispersion)
    ##-------------------------------------------
    if (model == "CP") {
        ##-------------------------------------------
        pars <- transform(pars, lambda = exp(eta),
                          nu = exp(dispersion))
        sumto <- list(...)$sumto
        ##-------------------------------------------
        approxmu <- with(pars, lambda^(1/nu) - (nu - 1)/(2 * nu))
        sigma <- with(pars, (1/nu) * approxmu)
        sigma <- ifelse(sigma < 0, 0, sigma)
        ymax <- ceiling(max(approxmu + 5 * sqrt(sigma)))
        etamax <- max(pars$eta)
        dispersionmin <- min(pars$dispersion)
        pmax <- exp(
            -llcmp(params = c(dispersionmin, etamax), y = ymax, X = 1,
                   sumto = sumto))
        while (pmax > tol) {
            ymax <- ymax + 1
            pmax <- exp(
                -llcmp(params = c(dispersionmin, etamax),
                        y = ymax, X = 1, sumto = sumto))
        }
        y <- 1:ymax
        estmean <- mapply(
            eta = pars$eta,
            nu = pars$nu,
            MoreArgs = list(y = y),
            FUN = function(eta, nu, y) {
                i <- 0:sumto
                zs <- sapply(eta, function(eta) {
                    sum(exp(i * eta - nu * lfactorial(i)))
                })
                py <- exp(y * eta - nu * lfactorial(y) - sum(log(zs)))
                sum(y * py)
            },
            SIMPLIFY = TRUE)
    }
    if (model == "CP2") {
        estmean <- exp(eta)
    }
    names(estmean) <- NULL
    return(c(estmean))
}

## Calculate the confidence interval of predict values
predictcm <- function(object, newdata, offset = NULL, level = 0.95) {
    ##-------------------------------------------
    if (missing(newdata)) {
        X <- object@data$X
    } else {
        X <- newdata
    }
    ##-------------------------------------------
    qn <- -qnorm((1 - level[1])/2)
    qn <- qn * c(fit = 0, lwr = -1, upr = 1)
    #--------------------------------------------
    V <- vcov(object)
    Vc <- V[-1, -1] - V[-1, 1] %*% solve(V[1, 1]) %*% V[1, -1]
    ## Vc <- V[-1, -1]
    eta <- cholV_eta(Vc, X, b = coef(object)[-1], qn = qn)
    ##-------------------------------------------
    args <- list(dispersion = coef(object)[1],
                 model = object@data$model,
                 offset = offset)
    if (object@data$model == "CP") {
        args$sumto <- object@data$sumto
    }
    pred <- apply(as.matrix(eta), MARGIN = 2,
                  FUN = function(x) {
                      do.call(what = calc_mean,
                              args = append(list(eta = x), args))
                  })
    colnames(pred) <- names(qn)
    return(pred)
}

##======================================================================
## Trellis panel functions for lattice graphics (use wzRfun)
library(wzRfun)

body(panel.cbH)[[14]][[3]][[3]] <-
    quote(panel.lines(xo, lyo, col = col.line, lty = list(...)$lty))
body(panel.cbH)[[14]][[3]][[4]] <-
    quote(panel.lines(xo, uyo, col = col.line, lty = list(...)$lty))

## To profile log-likelihood
myprofile <- function(model, which) {
    sel <- c("param", "z", "focal")
    prof <- profile(model, which = which)
    as.data.frame(prof)[, sel]
}

## To plot profile in lattice
xyprofile <- function(prof, conf = c(0.9, 0.95, 0.99),
                      namestrip = NULL, subset = 4,
                      scales.x = "free", ...) {
    ##-------------------------------------------
    conf <- conf[order(conf, decreasing = TRUE)]
    da <- subset(prof, abs(z) <= subset)
    ##-------------------------------------------
    fl <- levels(da$param)
    if (!is.null(namestrip)) {
        fl <- namestrip
    }
    xyplot(abs(z) ~ focal | param,
           data = da,
           layout = c(NA, 1),
           ## ylab = expression(abs(z)~~(sqrt(~Delta~"deviance"))),
           ylab = expression(sqrt(2~(L(hat(phi))-L(phi)))),
           scales = list(x = scales.x),
           type = c("l", "g"),
           strip = strip.custom(
               factor.levels = fl),
           panel = function(x, y, subscripts, ...) {
               conf <- c(0.9, 0.95, 0.99)
               hl <- sqrt(qchisq(conf, 1))
               ##-------------------------------------------
               mle <- x[y == 0]
               xl <- x[x < mle]; yl <- y[x < mle]
               xr <- x[x > mle]; yr <- y[x > mle]
               ##-------------------------------------------
               funleft <- approxfun(x = yl, y = xl)
               funright <- approxfun(x = yr, y = xr)
               vxl <- funleft(hl)
               vxr <- funright(hl)
               vz <- c(hl, hl)
               ##-------------------------------------------
               panel.xyplot(x, y, ...)
               panel.arrows(c(vxl, vxr), 0, c(vxl, vxr), vz,
                            code = 1, length = 0.1, lty = 2,
                            col = "gray40")
               panel.segments(vxl, vz, vxr, vz, lty = 2,
                              col = "gray40")
               panel.abline(h = 0, v = mle, lty = 3)
               panel.text(x = rep(mle, 2), y = vz + 0.1,
                          labels = paste(conf*100, "%"),
                          col = "gray20")
           }, ...)
}
