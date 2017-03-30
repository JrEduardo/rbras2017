## =====================================================================
## Analysis of number of grains in Soyabean under Soil Moisture and
## Potassium Doses
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
## data(nitrofen, package = "boot")
data(Paula, package = "labestData")
nitrofen <- PaulaEx4.6.20
nitrofen <- transform(nitrofen, dose = dose/100)

##----------------------------------------------------------------------
## Exploratory analysis

## Scatter plot
xy1 <- xyplot(novos ~ dose,
              data = nitrofen,
              xlab = "Nível de adubação potássica",
              ylab = "Número de grãos por parcela",
              type = c("p", "g", "smooth"),
              spread = 0.05,
              panel = panel.beeswarm)

## Sample variance vs sample mean (evidence in favor of the
## underdispersion).
## Sample variance vs sample mean (evidence in favor of the
## underdispersion).
mv <- aggregate(novos ~ dose, data = nitrofen,
                FUN = function(x) c(mean = mean(x), var = var(x)))
xlim <- ylim <- extendrange(c(mv$novos), f = 0.05)

xy2 <- xyplot(novos[, "var"] ~ novos[, "mean"],
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

print(xy1, split = c(1, 1, 2, 1), more = TRUE)
print(xy2, split = c(2, 1, 2, 1), more = FALSE)

##----------------------------------------------------------------------
## Fit models
mnames <- c("PO", "C1", "C2", "QP")

## Predictor
form <-  novos ~ dose + I(dose^2) + I(dose^3)

mPO <- glm(form, data = nitrofen, family = poisson)
mC1 <- fitcm(form, data = nitrofen, model = "CP", sumto = 100)
mC2 <- fitcm(form, data = nitrofen, model = "CP2", sumto = 100)
mQP <- glm(form, data = nitrofen, family = quasipoisson)

models.novos <- list(mPO, mC1, mC2, mQP)
names(models.novos) <- mnames

## Numbers of calls to loglik and numerical gradient
models.novos$C1@details$counts
models.novos$C2@details$counts

##----------------------------------------------------------------------
## LRT and profile extra parameter

## LRT between Poisson and COM-Poisson (test: phi == 0)
getAnova(mPO, mC2)

profs.novos <- lapply(list(c(mC1, "phi"), c(mC2, "phi2")),
                     function(x) myprofile(x[[1]], x[[2]]))
profs.novos <- do.call("rbind", profs.novos)

## Plot profile extra parameter
snames <- parse(text = c("'COM-Poisson'~(phi)",
                         "'COM-Poisson'[mu]~(phi)"))
xyprofile(profs.novos, namestrip = snames,
          ylim = extendrange(c(0, 3)),
          xlab = "Dispersion/Precision parameter",
          scales.x = NULL) +
    layer(panel.abline(v = 0, col = 2, lty = 2))

##----------------------------------------------------------------------
## Goodness of fit measures and estimate parameters

## GoF measures
measures.novos <- sapply(models.novos, function(x)
    c("LogLik" = logLik(x), "AIC" = AIC(x), "BIC" = BIC(x)))

## Get the estimates
est <- lapply(models.novos, FUN = function(x) getCoefs(x))
est.novos <- do.call(cbind, est)

##----------------------------------------------------------------------
## Prediction

## Data for prediction
pred <- with(nitrofen,
             data.frame("dose" = seq(min(dose), max(dose),
                                     length.out = 100)))
qn <- qnorm(0.975) * c(fit = 0, lwr = -1, upr = 1)

## Design matrix for prediction
X <- model.matrix(update(form, NULL~.), pred)

## Considering Poisson
aux <- exp(confint(
    glht(mPO, linfct = X), calpha = univariate_calpha())$confint)
colnames(aux) <- c("fit", "lwr", "upr")
aux <- data.frame(modelo = "Poisson", aux)
predPO.novos <- cbind(pred, aux)

## Considering COM-Poisson
aux <- predictcm(mC1, newdata = X)
aux <- data.frame(modelo = "COM-Poisson", aux)
predC1.novos <- cbind(pred, aux)

## Considering COM-Poisson (mean parametrization)
aux <- predictcm(mC2, newdata = X)
aux <- data.frame(modelo = "COM-Poisson2", aux)
predC2.novos <- cbind(pred, aux)

## Considering Quasi-Poisson
aux <- exp(confint(
    glht(mQP, linfct = X), calpha = univariate_calpha())$confint)
colnames(aux) <- c("fit", "lwr", "upr")
aux <- data.frame(modelo = "Quasi-Poisson", aux)
predQP.novos <- cbind(pred, aux)

## Representing the confidence intervals
pred.novos <- rbind(predPO.novos, predC1.novos, predC2.novos, predQP.novos)
ord <- order(pred.novos$dose, pred.novos$modelo)
pred.novos <- pred.novos[ord, ]

## Legend
key <- list(columns = 2,
            lines = list(col = cols[1:4], lty = rev(1:4), cex = 0.7),
            text = list(parse(
                text = c("'Poisson'", "'COM-Poisson'",
                         "'COM-Poisson'[mu]", "'Quasi-Poisson'"))
                )
            )

## Graph
update(xy1, layout = c(NA, 1), type = c("p", "g"),
       alpha = 0.6, key = key) +
    as.layer(
        xyplot(fit ~ dose,
               auto.key = TRUE,
               data = pred.novos,
               groups = modelo,
               type = "l",
               ly = pred.novos$lwr, uy = pred.novos$upr,
               cty = "bands", fill = "gray80", alpha = 0.1,
               panel = panel.superpose,
               panel.groups = panel.cbH,
               prepanel = prepanel.cbH,
               lty = rev(1:4))
    )

##----------------------------------------------------------------------
## Correlation between estimates
corr.novos <- purrr::map_df(models.novos[c("C1", "C2")],
                           function(x) cov2cor(vcov(x))[1, -1])
