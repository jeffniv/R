#' Simple Linear Regression with Uninformative Priors
#'
#' This function takes in a vector of predictor variables \code{x} and a vector of data \code{y} and creates a simple linear model using MCMC sampling. For the Bayesian inference, this method uses uninformative prior distributions. The point and interval estimates are taken from the posterior distribution from the samples.
#'
#' @author Jeff Nivitanont, \email{jeffniv@ou.edu}
#'
#' @param x A vector of predictor variables
#'
#' @param y A vector of data
#'
#' @param xnew A value of a predictor for making inference
#'
#' @param iter The number of iterations per chain
#'
#' @param burn The burn in period for chains
#'
#' @param thin The thinning interval for monitoring
#'
#' @param alpha For setting the \eqn{(1-\alpha)} Highest Density Interval
#'
#' @param parallel \code{TRUE/FALSE}. Option to use multiple CPUs
#'
#' @param DIC \code{TRUE/FALSE}. Option to track DIC. Cannot be used in conjunction with \code{parallel} method
#'
#' @examples
#' \dontrun{
#' data("electric")
#' y=electric[,3]
#' x=electric[,1]
#' testobj=UninformedSLR(x,
#'                       y,
#'                       xnew=1700,
#'                       iter=3334,
#'                       burn=5000,
#'                       thin=3,
#'                       alpha=0.05,
#'                       parallel=TRUE,
#'                       DIC=TRUE)
#' }
#' @export


UninformedSLR = function( x,
                          y,
                          xnew=NULL,
                          iter=3334,
                          burn=5000,
                          thin=3,
                          alpha=0.05,
                          parallel=TRUE,
                          DIC=FALSE){
  #check cores for parallel processing
  try( runjags::runjags.options( inits.warning=FALSE , rng.warning=FALSE ) )
  nChainsDefault = 3
  runjagsMethodDefault = "rjags" #default not parallel

  if(parallel){
    nCores = parallel::detectCores()
    if ( !is.finite(nCores) ) { nCores = 1 }
    if ( nCores > 4 ) {
      nChainsDefault = 4  # because JAGS has only 4 rng's.
      runjagsMethodDefault = "parallel"
    }
    if ( nCores == 4 ) {
      nChainsDefault = 3  # save 1 core for other processes.
      runjagsMethodDefault = "parallel"
    }
    if ( nCores < 4 ) {
      nChainsDefault = 3
      runjagsMethodDefault = "rjags" # NOT parallel
    }
  }
  if(DIC){
    if(parallel){
      print('Cannot run in parallel when tracking DIC. Switching to default mode.')
    }
    nChainsDefault = 3
    runjagsMethodDefault = "rjags" # NOT parallel
  }
  fileNameRoot="Jags-MATH4803Nivitanont-SLR" # For output file names.

  # Put the information into a list.
  Ntotal = length(y)
  dataList = list()
  if(!is.null(xnew)){
    dataList=list(
      y = y,
      x = x,
      xdash = xnew,
      Ntotal = Ntotal)
  } else{
    dataList=list(
      y = y,
      x = x,
      Ntotal = Ntotal)
  }
  # Define the model:
  modelString = paste0("
  # center the data
  data {
  ym <- mean(y)
  xm <- mean(x)
  for ( i in 1:Ntotal ) {
  zy[i] <- ( y[i] - ym )
  zx[i] <- ( x[i] - xm )
  }
  }
  model{
  for( i in 1 : Ntotal ) {
  zy[i] ~ dnorm(mu[i], tau)
  mu[i] <- zbeta0 + beta1 * zx[i]
  }
  zbeta0 ~ dnorm(0.0, 1.0E-6)
  beta1 ~ dnorm(0.0, 1.0E-6)
  sigma ~ dunif(0, 100)
  tau <- pow(sigma, -2)
  ",ifelse(!is.null(xnew),"#given x value
  ydash ~ dnorm(mudash, tau)
  mudash <- beta0 + beta1 * xdash",""),"
  #transform data from scaled
  beta0 <- zbeta0 + ym - beta1*xm
  }" # close quote for modelStringlose quote for modelString
  )
  writeLines( modelString , con="SLRmodel.txt" )

  variables=c('beta0', 'beta1', 'sigma' )
  if(!is.null(xnew)){
    variables=append(variables,'ydash')
  }
  if(DIC){
    variables=append(variables,'dic')
  }
  #start timer
  start.time <- proc.time()
  # Run the chains, output is runjags object:
  runJagsOut = runjags::run.jags(model = modelString,
                        data = dataList,
                        n.chains = 3,
                        monitor = variables,
                        burnin = burn,
                        adapt = 1000,
                        sample = iter,
                        thin = thin,
                        method = runjagsMethodDefault,
                        summarise = FALSE)
  #stop timer
  end.time <- proc.time()
  time.taken <- end.time - start.time
  cat(paste('Simulation took ',round(time.taken[3],2),' seconds to run.\n'))

  #convert runjags object to mcmc object
  codaSamples = coda::as.mcmc.list( runJagsOut )
  save( codaSamples , file=paste0(fileNameRoot,"Mcmc.Rdata") )

  #extract DIC from samples
  if(DIC){
    DICscore=runjags::extract.runjags(x=runJagsOut, what='dic')
    cat('DIC\n')
    print(DICscore)
  }

  #plot diagnostics
  diagGGmcmc(codaSamples)
  #point and interval estimates
  postsum=summary(codaSamples)
  coefs=postsum$statistics[,1]
  HPD=coda::HPDinterval(codaSamples, prob=1-alpha)

 
  #command line release
  cat('Point Estimates\n')
  print(coefs)
  cat('Interval Estimates\n')
  print(HPD[1])

  returnlist = list(sim=codaSamples, y=y, x=x, coefs=coefs, summary=postsum, HPD=HPD)
  if(DIC){
    returnlist = c(returnlist, DIC=DICscore)
  }
  return(returnlist)
}
