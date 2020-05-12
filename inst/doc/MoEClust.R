## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(fig.width=7, fig.height=7, fig.align = 'center', fig.show='hold',
                      dev.args=list(type="cairo"), warning=FALSE, message=FALSE, 
                      progress=FALSE, collapse=TRUE, comments="#>")

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages('devtools')
#  devtools::install_github('Keefe-Murphy/MoEClust')

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages('MoEClust')

## -----------------------------------------------------------------------------
library(MoEClust)

## -----------------------------------------------------------------------------
data(CO2data)
CO2   <- CO2data$CO2
GNP   <- CO2data$GNP

## ---- results='hide'----------------------------------------------------------
m1    <- MoE_clust(CO2, G=1:3, verbose=FALSE)
m2    <- MoE_clust(CO2, G=2:3, gating= ~ GNP, verbose=FALSE)
m3    <- MoE_clust(CO2, G=1:3, expert= ~ GNP, verbose=FALSE)
m4    <- MoE_clust(CO2, G=2:3, gating= ~ GNP, expert= ~ GNP, verbose=FALSE)
m5    <- MoE_clust(CO2, G=2:3, equalPro=TRUE, verbose=FALSE)
m6    <- MoE_clust(CO2, G=2:3, expert= ~ GNP, equalPro=TRUE, verbose=FALSE)

## -----------------------------------------------------------------------------
comp  <- MoE_compare(m1, m2, m3, m4, m5, m6, optimal.only=TRUE)

## -----------------------------------------------------------------------------
(mod1 <- MoE_stepwise(CO2, GNP, verbose=FALSE))

## -----------------------------------------------------------------------------
(mod2 <- MoE_stepwise(CO2, GNP, noise=TRUE, verbose=FALSE))

## -----------------------------------------------------------------------------
(best  <- MoE_compare(mod1, mod2, comp, pick=1)$optimal)

## ---- eval=FALSE--------------------------------------------------------------
#  (summ <- summary(best, classification=TRUE, parameters=FALSE, networks=FALSE))

## ---- echo=FALSE--------------------------------------------------------------
(summ <- suppressWarnings(summary(best)))

## -----------------------------------------------------------------------------
plot(best, what="gpairs", jitter=FALSE)

## ---- echo=FALSE--------------------------------------------------------------
res    <- best
G      <- res$G
expert <- res$expert
x.name <- names(res$net.covs)
y.name <- names(res$data)
x      <- as.matrix(res$net.covs)
y      <- as.matrix(res$data)
plot(x=x, y=y, main=substitute(atop(paste('CO'[2], " Data"), paste(Mname, " model, G=", rG, ", equal mixing proportions, incl. expert network covariate: GNP")), list(Mname=res$modelName, rG=G)), ylab=expression('CO'[2]), xlab="GNP", type="n", cex.main=0.95)
x.new  <- setNames(seq(par("usr")[1], par("usr")[2], length=100), rep("GNP", 100))
y.new  <- setNames(seq(par("usr")[3], par("usr")[4], length=100), rep("CO2", 100))
grid   <- expand.grid(x.new, y.new)
getden <- function(x, y, res) {
  sig  <- res$parameters$variance$modelName
  den  <- do.call(cbind, lapply(seq_len(G), function(k) dnorm(y, predict(expert[[k]], newdata=setNames(data.frame(x), names(x)[1]), type="response"), sqrt(ifelse(sig == "V", res$parameters$variance$sigmasq[k], res$parameters$variance$sigmasq)))))
  rowSums(matrix(res$parameters$pro, nrow=length(x), ncol=G, byrow=TRUE) * den)
}
mat    <- matrix(getden(grid[,1L], grid[,2L], res), length(x.new), length(y.new))
image(x.new, y.new, mat, col=c("white", heat.colors(30L)[30L:1L]), xlab="GNP", ylab=paste('CO'[2]), add=TRUE)
box(lwd=1)
contour(x.new, y.new, mat, add=TRUE, col="grey60", labcex=0.75)
points(GNP, CO2)

## ---- fig.height=5.5, fig.width=5.5-------------------------------------------
mod <- as.Mclust(comp$optimal)
plot(mod, what="classification")
plot(mod, what="uncertainty")

## ---- eval=FALSE--------------------------------------------------------------
#  as.vector(predict(comp$optimal)$y)

## ---- echo=FALSE--------------------------------------------------------------
as.vector(suppressWarnings(predict(comp$optimal)$y))

## ---- results="hide"----------------------------------------------------------
ind     <- sample(1:nrow(CO2data), 2)
res2    <- MoE_clust(CO2data[-ind,]$CO2, G=3, expert=~GNP, 
                     equalPro=TRUE, network.data=CO2data[-ind,])

## ---- eval=FALSE--------------------------------------------------------------
#  # Using new covariates only
#  predict(res2, newdata = CO2data[ind,], use.y = FALSE)[1:3]
#  
#  # Using both new covariates & new response data
#  predict(res2, newdata = CO2data[ind,])[1:3]

## ---- echo=FALSE--------------------------------------------------------------
# Using new covariates only
suppressWarnings(predict(res2, newdata = CO2data[ind,], use.y = FALSE)[1:3])

# Using both new covariates & new response data
predict(res2, newdata = CO2data[ind,])[1:3]         

## -----------------------------------------------------------------------------
data(ais)
hema  <- ais[,3:7]

## ---- eval=FALSE--------------------------------------------------------------
#  ?MoE_control

## ---- eval=FALSE--------------------------------------------------------------
#  mod   <- MoE_clust(hema, G=1:3, expert= ~ sex, gating= ~ BMI,
#                     network.data=ais, tau0=0.1, noise.gate=FALSE)

## ---- echo=FALSE--------------------------------------------------------------
mod   <- suppressWarnings(MoE_clust(hema, G=1:3, expert= ~ sex, gating= ~ BMI, network.data=ais, tau0=0.1, noise.gate=FALSE, verbose=FALSE))

## -----------------------------------------------------------------------------
plot(mod, what="gpairs")

## -----------------------------------------------------------------------------
plot(mod, what="gpairs", response.type="density")

## -----------------------------------------------------------------------------
plot(mod, what="gpairs", response.type="uncertainty")

## -----------------------------------------------------------------------------
plot(mod, what="uncertainty", type="profile")

## ---- echo=FALSE--------------------------------------------------------------
plot(mod, what="criterion", legendArgs=list(x="right"))

## -----------------------------------------------------------------------------
plot(mod, what="gating", x.axis=ais$BMI, xlab="BMI")

## ---- eval=FALSE--------------------------------------------------------------
#  plot(mod, what="loglik")

## ---- echo=FALSE--------------------------------------------------------------
plot(mod, what="loglik")

## ---- eval=FALSE--------------------------------------------------------------
#  require("lattice")
#  z <- factor(mod$classification[mod$classification > 0],
#              labels=paste0("Cluster", seq_len(mod$G)))
#  splom(~ hema | ais$sex, groups=z)
#  splom(~ hema | z, groups=ais$sex)

## ---- echo=FALSE, fig.height=4------------------------------------------------
require("lattice", quietly=TRUE)
z <- factor(mod$classification[mod$classification > 0], labels=paste0("Cluster", seq_len(mod$G)))
splom(~ hema | ais$sex, groups=z, xlab=NULL)
splom(~ hema | z, groups=ais$sex, xlab=NULL)

