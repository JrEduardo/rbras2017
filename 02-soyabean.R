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
data(soyaBeans, package = "cmpreg")
soyaBeans <- soyaBeans[-74, ] ## Incorrect observation
soyaBeans <- transform(soyaBeans, K = K/100)

##----------------------------------------------------------------------
## Exploratory analysis

## Scatter plot
xy1 <- xyplot(ngra ~ K | umid,
              data = soyaBeans,
              xlab = "Nível de adubação potássica",
              ylab = "Número de grãos por parcela",
              type = c("p", "g", "smooth"),
              as.table =  TRUE,
              layout = c(2, 2),
              strip = strip.custom(
                  strip.names = TRUE, var.name = "Umidade",
                  factor.levels = paste0(levels(soyaBeans$umid), "%")))

## Sample variance vs sample mean (evidence in favor of the
## underdispersion).
## Sample variance vs sample mean (evidence in favor of the
## underdispersion).
mv <- aggregate(ngra ~ K + umid, data = soyaBeans,
                FUN = function(x) c(mean = mean(x), var = var(x)))
xlim <- ylim <- extendrange(c(mv$ngra), f = 0.05)

xy2 <- xyplot(ngra[, "var"] ~ ngra[, "mean"],
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

## Predictor
form <-  ngra ~ bloc + umid * K + I(K^2)

mPO <- glm(form, data = soyaBeans, family = poisson)
mC1 <- fitcm(form, data = soyaBeans, model = "CP", sumto = 700)
mC2 <- fitcm(form, data = soyaBeans, model = "CP2", sumto = 700)
mQP <- glm(form, data = soyaBeans, family = quasipoisson)

models.ngra <- list(mPO, mC1, mC2, mQP)
names(models.ngra) <- mnames

## Numbers of calls to loglik and numerical gradient
models.ngra$C1@details$counts
models.ngra$C2@details$counts

##----------------------------------------------------------------------
## LRT and profile extra parameter

## LRT between Poisson and COM-Poisson (test: phi == 0)
getAnova(mPO, mC2)

profs.ngra <- lapply(list(c(mC1, "phi"), c(mC2, "phi2")),
                     function(x) myprofile(x[[1]], x[[2]]))
profs.ngra <- do.call("rbind", profs.ngra)

## Plot profile extra parameter
snames <- parse(text = c("'COM-Poisson'~(phi)",
                         "'COM-Poisson'[mu]~(phi)"))
xyprofile(profs.ngra, namestrip = snames,
          ylim = extendrange(c(0, 3)),
          xlab = "Dispersion/Precision parameter",
          scales.x = NULL)

##----------------------------------------------------------------------
## Goodness of fit measures and estimate parameters

## GoF measures
measures.ngra <- sapply(models.ngra, function(x)
    c("LogLik" = logLik(x), "AIC" = AIC(x), "BIC" = BIC(x)))

## Get the estimates
est <- lapply(models.ngra, FUN = function(x) getCoefs(x))
est.ngra <- do.call(cbind, est)

##----------------------------------------------------------------------
## Prediction

## Data for prediction
pred <- with(soyaBeans,
             expand.grid(
                 bloc = factor(levels(bloc)[1], levels = levels(bloc)),
                 umid = levels(umid),
                 K = seq(min(K), max(K), l = 20)
             ))
qn <- qnorm(0.975) * c(fit = 0, lwr = -1, upr = 1)

## Design matrix for prediction
X <- model.matrix(update(form, NULL~.), pred)
bl <- attr(X, "assign") == 1
X[, bl] <- X[, bl] + 1/(sum(bl) + 1)

## Considering Poisson
aux <- exp(confint(
    glht(mPO, linfct = X), calpha = univariate_calpha())$confint)
colnames(aux) <- c("fit", "lwr", "upr")
aux <- data.frame(modelo = "Poisson", aux)
predPO.ngra <- cbind(pred, aux)

## Considering COM-Poisson
aux <- predictcm(mC1, newdata = X)
aux <- data.frame(modelo = "COM-Poisson", aux)
predC1.ngra <- cbind(pred, aux)

## Considering COM-Poisson (mean parametrization)
aux <- predictcm(mC2, newdata = X)
aux <- data.frame(modelo = "COM-Poisson2", aux)
predC2.ngra <- cbind(pred, aux)

## Considering Quasi-Poisson
aux <- exp(confint(
    glht(mQP, linfct = X), calpha = univariate_calpha())$confint)
colnames(aux) <- c("fit", "lwr", "upr")
aux <- data.frame(modelo = "Quasi-Poisson", aux)
predQP.ngra <- cbind(pred, aux)

## Representing the confidence intervals
pred.ngra <- rbind(predPO.ngra, predC1.ngra, predC2.ngra, predQP.ngra)

## Legend
key <- list(columns = 2,
            lines = list(col = cols[1:4], lty = rev(1:4)),
            text = list(parse(
                text = c("'Poisson'", "'COM-Poisson'",
                         "'COM-Poisson'[mu]", "'Quasi-Poisson'"))
                )
            )

## Graph
update(xy1, layout = c(NA, 1), type = c("p", "g"),
       alpha = 0.6, key = key) +
    as.layer(
        xyplot(fit ~ K | umid,
               data = pred.ngra,
               groups = modelo,
               type = "l",
               ly = pred.ngra$lwr, uy = pred.ngra$upr,
               cty = "bands", fill = "gray80", alpha = 0.1,
               panel = panel.superpose,
               panel.groups = panel.cbH,
               prepanel = cmpreg::prepanel.cbH,
               lty = rev(1:4))
    )

##----------------------------------------------------------------------
## Correlation between estimates
corr.ngra <- purrr::map_df(models.ngra[c("C1", "C2")],
                           function(x) cov2cor(vcov(x))[1, -1])


source("functions.R")
library(bbmle)

soil <- paste0(c("37.5", "50.0", "62.5"), "%")
soyabeans <- expand.grid(K = c(0L, 30L, 60L, 120L, 180L),
                         soil = factor(x = soil, levels = soil),
                         bloc = factor(c("I", "II", "III", "IV", "V")),
                         KEEP.OUT.ATTRS = FALSE)

## Response variables
soyabeans$npod <- c(56L, 62L, 66L, 68L, 82L, 63L, 86L, 94L, 86L, 97L,
                    58L, 80L, 99L, 110L, 106L, 48L, 64L, 62L, 60L, 62L,
                    52L, 103L, 103L, 87L, 101L, 64L, 77L, 96L, 94L,
                    105L, 44L, 63L, 84L, 66L, 58L, 56L, 81L, 92L, 91L,
                    100L, 50L, 71L, 94L, 89L, 89L, 60L, 75L, 67L, 76L,
                    63L, 51L, 60L, 67L, 90L, 80L, 62L, 67L, 75L, 90L,
                    87L, 36L, 57L, 53L, 69L, 50L, 57L, 61L, 92L, 96L,
                    110L, 69L, 63L, 95L, 55L, 90L)
soyabeans$ngra <- c(136L, 159L, 156L, 171L, 190L, 140L, 193L, 200L,
                    208L, 237L, 132L, 178L, 231L, 271L, 262L, 120L,
                    146L, 159L, 152L, 159L, 124L, 239L, 248L, 190L,
                    240L, 140L, 167L, 237L, 240L, 248L, 106L, 157L,
                    204L, 169L, 147L, 133L, 199L, 211L, 223L, 261L,
                    128L, 164L, 216L, 215L, 228L, 125L, 185L, 159L,
                    176L, 154L, 133L, 149L, 170L, 215L, 186L, 147L,
                    154L, 179L, 229L, 216L, 92L, 133L, 139L, 166L, 128L,
                    131L, 150L, 218L, 236L, 261L, 150L, 156L, 215L,
                    128L, 208L)

## Remove the observation number 74, considered atypical case by
## researcher
soyabeans <- soyabeans[-74, ]
soyabeans <- transform(soyabeans, K = K/100)

predictor2 <-  ~ bloc + soil * K + I(K^2)


source("functions.R")
library(bbmle)
m2CP.ngra <- fitcm(update(predictor2, ngra ~ .), data = soyabeans,
                   model = "CP", sumto = 700)
