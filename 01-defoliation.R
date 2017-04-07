## =====================================================================
## Analysis of Artificial Defoliation in Cotton Plants
##                                                        Eduardo Junior
##                                                    edujrrib@gmail.com
##                                                            2017-03-16
## =====================================================================

##----------------------------------------------------------------------
## Load package and functions

library(bbmle)
library(multcomp)
library(plyr)

source("config.R")
source("functions.R")

## Colors for legends
cols <- trellis.par.get("superpose.line")$col

##----------------------------------------------------------------------
## Load data
data(cottonBolls, package = "cmpreg")
str(cottonBolls)

##----------------------------------------------------------------------
## Exploratory analysis

## Scatter plot with a beewarm taste
xy1 <- xyplot(ncap ~ des | est,
              data = cottonBolls,
              layout = c(2, 3),
              as.table = TRUE,
              type = c("p", "smooth", "g"),
              xlab = "Níveis de desfolha artificial",
              ylab = "Número de capulhos produzidos",
              spread = 0.05,
              panel = panel.beeswarm)

## Sample variance vs sample mean (evidence in favor of the
## underdispersion).
mv <- aggregate(ncap ~ est + des, data = cottonBolls,
                FUN = function(x) c(mean = mean(x), var = var(x)))
xlim <- ylim <- extendrange(c(mv$ncap), f = 0.05)

xy2 <- xyplot(ncap[, "var"] ~ ncap[, "mean"],
              data = mv,
              type = c("p", "r", "g"),
              xlim = xlim,
              ylim = ylim,
              xlab = expression("Média"~"amostral"~(bar(y))),
              ylab = expression("Variância"~"amostral"~(s^2)),
              panel = function(...) {
                  panel.xyplot(...)
                  panel.abline(a = 0, b = 1, lty = 2)
              })

#+ fig.height=4.5, fig.width=9
print(xy1, split = c(1, 1, 2, 1), more = TRUE)
print(xy2, split = c(2, 1, 2, 1), more = FALSE)

##----------------------------------------------------------------------
## Fit models

mnames <- c("PO", "C1", "C2", "QP")

## Predictor, following Zeviani et al. (2014)
form1 <- ncap ~ est:(des + I(des^2))

m1PO <- glm(form1, data = cottonBolls, family = poisson)
system.time(
    m1C1 <- fitcm(form1, data = cottonBolls, model = "CP", sumto = 50)
)
system.time(
    m1C2 <- fitcm(form1, data = cottonBolls, model = "CP2", sumto = 50)
)
m1QP <- glm(form1, data = cottonBolls, family = quasipoisson)

models.ncap <- list(m1PO, m1C1, m1C2, m1QP)
names(models.ncap) <- mnames

## Numbers of calls to loglik and numerical gradient
models.ncap$C1@details$counts
models.ncap$C2@details$counts

##----------------------------------------------------------------------
## LRT and profile extra parameter

## LRT between Poisson and COM-Poisson (test: phi == 0)
getAnova(m1PO, m1C2)

profs.ncap <- lapply(list(c(m1C1, "phi"), c(m1C2, "phi2")),
                     function(x) myprofile(x[[1]], x[[2]]))
profs.ncap <- do.call("rbind", profs.ncap)

## Plot profile extra parameter
snames <- parse(text = c("'COM-Poisson'~(phi)",
                         "'COM-Poisson'[mu]~(phi)"))
xyprofile(profs.ncap, namestrip = snames,
          ylim = extendrange(c(0, 3)),
          xlab = "Dispersion/Precision parameter")

##----------------------------------------------------------------------
## Goodness of fit measures and estimate parameters

## GoF measures
measures.ncap <- sapply(models.ncap, function(x)
    c("LogLik" = logLik(x), "AIC" = AIC(x), "BIC" = BIC(x)))
measures.ncap

## Get the estimates
est <- lapply(models.ncap, FUN = function(x) getCoefs(x))
est.ncap <- do.call(cbind, est)
est.ncap

##----------------------------------------------------------------------
## Prediction

## Data for prediction
pred <- with(cottonBolls,
             expand.grid(
                 est = levels(est),
                 des = seq(min(des), max(des), l = 20)
             ))
qn <- qnorm(0.975) * c(fit = 0, lwr = -1, upr = 1)

## Design matrix for prediction
X <- model.matrix(update(form1, NULL~.), pred)

## Considering Poisson
aux <- exp(confint(
    glht(m1PO, linfct = X), calpha = univariate_calpha())$confint)
colnames(aux) <- c("fit", "lwr", "upr")
aux <- data.frame(modelo = "Poisson", aux)
predPO.ncap <- cbind(pred, aux)

## Considering COM-Poisson
aux <- predictcm(m1C1, newdata = X)
aux <- data.frame(modelo = "COM-Poisson", aux)
predC1.ncap <- cbind(pred, aux)

## Considering COM-Poisson (mean parametrization)
aux <- predictcm(m1C2, newdata = X)
aux <- data.frame(modelo = "COM-Poisson2", aux)
predC2.ncap <- cbind(pred, aux)

## Considering Quasi-Poisson
aux <- exp(confint(
    glht(m1QP, linfct = X), calpha = univariate_calpha())$confint)
colnames(aux) <- c("fit", "lwr", "upr")
aux <- data.frame(modelo = "Quasi-Poisson", aux)
predQP.ncap <- cbind(pred, aux)

## Representing the confidence intervals
pred.ncap <- rbind(predPO.ncap, predC1.ncap, predC2.ncap, predQP.ncap)

## Legend
key <- list(columns = 2,
            lines = list(col = cols[1:4], lty = rev(1:4)),
            text = list(parse(
                text = c("'Poisson'", "'COM-Poisson'",
                         "'COM-Poisson'[mu]", "'Quasi-Poisson'"))
                ))

#+ fig.height=3.8, fig.width=8
## Graph
xyplot(ncap ~ des | est,
       data = cottonBolls,
       layout = c(NA, 1),
       as.table = TRUE,
       type = c("p", "g"),
       xlab = "Níveis de desfolha artificial",
       ylab = "Número de capulhos produzidos",
       spread = 0.05,
       alpha = 0.6, key = key,
       panel = panel.beeswarm) +
    as.layer(
        xyplot(fit ~ des | est,
               auto.key = TRUE,
               data = pred.ncap,
               groups = modelo,
               type = "l",
               layout = c(NA, 1),
               as.table = TRUE,
               ly = pred.ncap$lwr, uy = pred.ncap$upr,
               cty = "bands", fill = "gray80", alpha = 0.1,
               panel = panel.superpose,
               panel.groups = panel.cbH,
               prepanel = cmpreg::prepanel.cbH,
               lty = rev(1:4))
    )

##----------------------------------------------------------------------
## Correlation between estimates
corr.ncap <- purrr::map_df(models.ncap[c("C1", "C2")],
                           function(x) cov2cor(vcov(x))[1, -1])
corr.ncap
