## =====================================================================
## Make paper companion
##                                                        Eduardo Junior
##                                                    edujrrib@gmail.com
##                                                            2017-04-06
## =====================================================================

library(legdown)
setwd("..")

files <- c(
    "General functions" = "functions.R",
    "Graphics configs" = "config.R",
    "Moments Approximation" = "00-approx.R",
    "Defoliation experiment" = "01-defoliation.R",
    "Soybean experiment" = "02-soyabean.R",
    "Nitrofen experiment" = "03-nitrofen.R"
)
authors <- c(
    "[Eduardo Elias Ribeiro Junior](http://www.leg.ufpr.br/~eduardojr)",
    "[Walmes Marques Zeviani](http://www.leg.ufpr.br/~walmes)",
    "[Wagner Hugo Bonat](http://www.leg.ufpr.br/~wagner)",
    "[Clarice Garcia Borges Demétrio](http://www4.esalq.usp.br/pesquisa/node/25)"
)
title <- paste(
    "Reparametrization of COM-Poisson Regression Models with",
    "Applications in the Analysis of Experimental Count Data"
)

## Generate paper companion
make_papercomp(path = ".",
               r_files = files,
               title = title,
               authors = authors,
               date = "March 31, 2017",
               output_dir = "papercomp")

##-------------------------------------------
