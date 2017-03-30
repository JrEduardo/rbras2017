## =====================================================================
## Analysis of Approximation moments
##                                                        Eduardo Junior
##                                                    edujrrib@gmail.com
##                                                            2017-03-16
## =====================================================================

## https://github.com/noamross/cmp/blob/master/R/cmp_fit.R
## https://github.com/James-Thorson/Conway-Maxwell-Poisson/tree/master/R

##----------------------------------------------------------------------
## Load package and functions

source("config.R")
source("functions.R")

## Colors for legends
cols <- trellis.par.get("superpose.line")$col

##----------------------------------------------------------------------
## COM-Poisson probability mass function (mean parametrization)
dcmp <- Vectorize(FUN = function(y, mu, phi, sumto = 100)
    exp(-llcmp2(c(phi, log(mu)), y = y, X = 1, sumto = sumto)),
    vectorize.args = c("y", "mu", "phi"))

grid <- expand.grid(mu = c(2, 8, 15), phi = log(c(0.5, 1, 2.5)))
y <- 0:30
py <- mapply(FUN = dcmp,
             mu = grid$mu,
             phi = grid$phi,
             MoreArgs = list(y = y, sumto = 100),
             SIMPLIFY = FALSE)
grid <- cbind(grid[rep(1:nrow(grid), each = length(y)), ],
              y = y,
              py = unlist(py))

## COM-Poisson p.m.f. to different combination betwenn phi and mu
leg_phi <- parse(
    text = paste("phi == \"",
                 formatC(unique(grid$phi), 1, format = "f"),
                 "\""))
barchart(py ~ y | factor(mu),
         groups = factor(phi),
         data = grid,
         horizontal = FALSE,
         layout = c(NA, 1),
         as.table = TRUE,
         axis = axis.grid,
         origin = 0,
         xlim = extendrange(y, f = 0.01),
         border = "transparent",
         scales = list(x = list(at = pretty(y))),
         ylab = expression(P(Y==y)),
         xlab = expression(y),
         par.settings = list(
             superpose.polygon = list(col = c(mycol[2:4])),
             superpose.line = list(col = c(mycol[2:4]))
         ),
         auto.key = list(
             columns = 3,
             rectangles = FALSE,
             lines = TRUE,
             text = leg_phi
         ),
         strip = strip.custom(
             strip.names = TRUE,
             var.name = expression(mu == ""),
             sep = ""))

##----------------------------------------------------------------------
## Study the approximation

##-------------------------------------------
## Mean and variance relationship
aux <- expand.grid(
    mu = seq(2, 30, length.out = 50),
    phi = seq(log(0.3), log(2.5), length.out = 50))

moments <- mapply(FUN = calc_moments,
                  mu = aux$mu,
                  phi = aux$phi,
                  SIMPLIFY = FALSE)
grid <- cbind(aux, t(do.call(cbind, moments)))
grid <- transform(grid, va = mu / exp(phi))

col <- brewer.pal(n = 8, name = "RdBu")
myreg <- colorRampPalette(colors = col)
xy1 <- xyplot(var ~ mean,
              groups = phi,
              data = grid,
              type = "l",
              lwd = 2,
              axis = axis.grid,
              xlab = expression(E(X)),
              ylab = expression(V(X)),
              xlab.top = expression(phi),
              legend = list(
                  top = list(
                      fun = draw.colorkey,
                      args = list(
                          key = list(
                              space = "top",
                              col = myreg(length(unique(grid$phi))),
                              at = unique(grid$phi),
                              draw = FALSE)))),
              par.settings = list(
                  superpose.line = list(
                      col = myreg(length(unique(grid$phi)))),
                  layout.heights = list(
                      key.axis.padding = -1)),
              panel = function(x, y, ...) {
                  panel.xyplot(x, y, ...)
                  panel.curve(1*x, min(x), max(x), lty = 2)
              })

##-------------------------------------------
## Errors in approximations for E[Y] and V[Y]
grid <- transform(grid,
                  emu = (mu - mean)^2,
                  eva = (va - var)^2)

myreg <- colorRampPalette(c("gray90",  "gray20"))(100)
xy2 <- levelplot(emu ~ phi + mu, data = grid,
                 col.regions = myreg,
                 xlab = expression(phi),
                 ylab = expression(mu),
                 colorkey = list(space = "top"),
                 panel = function(x, y, z, ...) {
                     panel.levelplot(x, y, z, ...)
                     panel.curve(10 - ( exp(x) - 1)/(2 * exp(x)),
                                 lty = 2)
                     panel.abline(v = 0, lty = 2)
                 })

xy3 <- levelplot(eva ~ phi + mu, data = grid,
                 col.regions = myreg,
                 xlab = expression(phi),
                 ylab = expression(mu),
                 colorkey = list(space = "top"),
                 panel = function(x, y, z, ...) {
                     panel.levelplot(x, y, z, ...)
                     panel.curve((10 - ( exp(x) - 1)/
                                  (2 * exp(x)))/exp(x), lty = 2)
                     panel.abline(v = 0, lty = 2)
                 })

print(xy1, split = c(1, 1, 3, 1), more = TRUE)
print(xy2, split = c(2, 1, 3, 1), more = TRUE)
print(xy3, split = c(3, 1, 3, 1), more = FALSE)
