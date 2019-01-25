#' LASSO Regression using MCMC Sampling
#'
#' This function takes in a matrix of predictor variables \code{x} and a vector of data \code{y} and runs a LASSO regression. This is especially useful for dealing with multicollinearity issues with the goal of eliminating predictor variables from the model.
#'
#' @author Jeff Nivitanont, \email{jeffniv@ou.edu}
#'
#' @param x A matrix or data frame of predictor variables
#'
#' @param y A vector of data
#'
#' @param xnew A vector of predictors for making inference
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
#' @param phi The constraint value on predictor priors
#'
#' @param DIC \code{TRUE/FALSE}. Option to track DIC. Cannot be used in conjunction with \code{parallel} method
#'
#' @examples
#' \dontrun{
#' data(clerical)
#' y=clerical[,1]
#' x=as.matrix(clerical[,-1])
#' xnew=c(5000, 75, 900, 200, 650)
#' testobj=LASSOReg( x,
#'             y,
#'             xnew,
#'             iter=3334,
#'             burn=5000,
#'             thin=3,
#'             alpha=0.05,
#'             parallel=TRUE,
#'             phi=NULL,
#'             DIC=FALSE)
#'
#' ## Point Estimates
#' ##beta0       beta[1]       beta[2]       beta[3]       beta[4]       beta[5]         sigma           phi         ydash
#' ##67.053598720   0.001220227   0.108111808   0.012961810  -0.043823459   0.055245865  11.225304473   0.500000000 120.176385903
#' ##Interval Estimates
#' ##[[1]]
#' ##lower        upper
#' ##beta0   49.315400000  84.33490000
#' ##beta[1] -0.000400664   0.00300408
#' ##beta[2]  0.017713400   0.19708500
#' ##beta[3] -0.005275380   0.02977720
#' ##beta[4] -0.078546200  -0.00877705
#' ##beta[5]  0.034261800   0.07875780
#' ##sigma    9.130320000  13.81670000
#' ##phi      0.500000000   0.50000000
#' ##ydash   97.900300000 141.97700000
#' ##attr(,"Probability")
#' ##[1] 0.94991
#' }
#' @export


LASSOReg = function( x,
                     y,
                     xnew=NULL,
                     iter=3334,
                     burn=5000,
                     thin=3,
                     alpha=0.05,
                     parallel=TRUE,
                     phi=NULL,
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

  fileNameRoot="Jags-MATH4803Nivitanont-LASSO" # For output file names.

  # Put the information into a list.
  Ntotal = length(y)
  Nx = length(x[1,])
  dataList = list()
  if(!is.null(xnew)){
    dataList=list(
      y = y,
      x = x,
      xdash = xnew,
      Nx = Nx,
      Ntotal = Ntotal)
  } else{
    dataList=list(
      y = y,
      x = x,
      Nx = Nx,
      Ntotal = Ntotal)
  }
  # Define the model:
  modelString = paste0("
  # center the data
  data {
  ym <- mean(y)
  for ( i in 1:Ntotal ) {
  zy[i] <- ( y[i] - ym )
  }
  for ( j in 1:Nx ) {
  xm[j]  <- mean(x[,j])

  for ( i in 1:Ntotal ) {
  zx[i,j] <- ( x[i,j] - xm[j] )
  }
  }
  }

  model{
  for( i in 1 : Ntotal ) {
  zy[i] ~ dnorm( mu[i], tau )
  mu[i] <- zbeta0 + sum( beta[1:Nx] * zx[i,1:Nx] )
  }
  for ( j in 1 : Nx ) {
  beta[j] ~ ddexp(0.0, 1/(phi^2))     # LASSO regression priors
  }
  tau <- pow(sigma, -2)
  sigma ~ dgamma(1.0E-2, 1.0E-2)
  zbeta0 ~ ddexp(0.0, 1/(phi^2))
  ", ifelse(is.null(phi), "phi ~ dgamma(1.0E-2, 1.0E-2)",paste0('phi<-',phi)),
  ifelse(!is.null(xnew),"
  #y given x
  ydash ~ dnorm( mudash, tau)
  mudash <- beta0 + sum( beta[1:Nx] * xdash[1:Nx] )",""), "
  #transform from scaled data
  beta0 <- zbeta0 + ym - sum(beta[1:Nx] * xm[1:Nx])
  }

  " # close quote for modelStringlose quote for modelString
  )#close parenthesis for paste0
  writeLines( modelString , con="LASSOmodel.txt" )

  variables=c('beta0', paste0('beta[',1:Nx,']'), 'sigma', 'phi' )
  if(!is.null(xnew)){
    variables=append(variables,'ydash')
  }
  if(DIC){
    variables=append(variables,'dic')
  }
  #jagsModel = jags.model( file="TEMPass4model.txt" , data=dataList , n.chains=3 , n.adapt=500 )
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

  #plot GGMCMC diagnostics
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
