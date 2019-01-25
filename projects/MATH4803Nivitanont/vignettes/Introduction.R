## ---- eval=FALSE,collapse=TRUE-------------------------------------------
#  InformedBetaBF(alpha2=2, beta2=5, n=10, x=4)

## ---- echo=FALSE,collapse=TRUE, fig.width=8,fig.height=8-----------------
library(MATH4803Nivitanont)
InformedBetaBF2(alpha2=2, beta2=5, n=10, x=4)

## ---- eval=FALSE---------------------------------------------------------
#  data(clerical)
#  corrcheck(clerical)

## ---- echo=FALSE, fig.width=8,fig.height=8, fig.show='hold'--------------
library(MATH4803Nivitanont)
data(clerical)
corrcheck1(clerical)
corrcheck2(clerical)

## ----eval=FALSE, fig.show='hold'-----------------------------------------
#  y=clerical[,1]
#  x=as.matrix(clerical[,-1])
#  xnew=c(5000, 75, 900, 200, 650)
#  testobj=RidgeReg(x,y,xnew, DIC=TRUE)
#  

## ---- echo=FALSE, fig.show='hold', collapse=TRUE, fig.width=8,fig.height=8----
library(MATH4803Nivitanont)
data(codaSamples)
mcmc=coda::as.mcmc.list(codaSamples)
diag1(mcmc)


## ---- echo=FALSE, fig.show='hold', collapse=TRUE, fig.width=8,fig.height=8----
library(MATH4803Nivitanont)
data(codaSamples)
mcmc=coda::as.mcmc.list(codaSamples)
diag2(mcmc)

