#' General Linear Model MCMC Sampling
#'
#' This function takes in a matrix of predictor variables \code{x} and a vector of data \code{y} and runs a general linear model using uninformative priors. It will require many iterations for the autocorrelation to shrink if the predictors have high multicollinearity. Before using this function, it is suggested to use the \code{corrcheck()} function to determine whether it is more appropriate to use the \code{RidgeReg()} function or the \code{LASSOReg()} functions.
#'
#' @author Jeff Nivitanont, \email{jeffniv@ou.edu}
#'
#' @param x A matrix or data frame of predictor variables
#'
#' @param y A vector of data
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
#' data('leukemia')
#' y=leukemia[,1]
#' x=as.matrix(leukemia[,-1])
#' testobj=   UninformedGLM( x,
#'                           y,
#'                           iter=3334,
#'                           burn=5000,
#'                           thin=3,
#'                           alpha=0.05,
#'                           parallel=TRUE,
#'                           DIC=FALSE)
#' }
#' @export


UninformedGLM = function( x,
                          y,
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
  fileNameRoot="Jags-MATH4803Nivitanont-GLM" # For output file names.

  # Put the information into a list.
  Ntotal = length(y)
  Nx = length(x[1])
  dataList=list(
    y = y,
    x = x,
    Nx = Nx,
    Ntotal = Ntotal)

  # Define the model:
  modelString = "
    # center the data
    data {
      for ( j in 1:Nx ) {
        xm[j]  <- mean(x[,j])
        for ( i in 1:Ntotal ) {
          zx[i,j] <- ( x[i,j] - xm[j] )
        }
      }
    }

    model{
      for( i in 1 : Ntotal ) {
        y[i] ~ dbern( P[i] )
        logit(P[i]) <- zbeta0 + sum( beta[1:Nx] * zx[i,1:Nx] )
      }
      for ( j in 1 : Nx ) {
        beta[j] ~ dnorm(0.0, 1.0E-6)     # Uninformative regression priors
      }
      tau <- pow(sigma, -2)
      sigma ~ dunif(0, 100)
      zbeta0 ~ dnorm(0.0, 1.0E-6)

      #transform from scaled data
      beta0 <- zbeta0 - sum(beta[1:Nx] * xm[1:Nx])
    }

    " # close quote for modelString, paste0
  writeLines( modelString , con="GLMModel.txt" )

  variables=c('beta0', paste0('beta[',1:Nx,']'), 'sigma', 'P' )
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

  #plot GGMCMC diagnostics
  diagGGmcmc(codaSamples[,1:(Nx+2)])
  
  #caterpillar plot
  ggm=ggmcmc::ggs(codaSamples[,c(-(1:(Nx+2)))])
  dev.new(noRStudioGD = TRUE)
  catplot=ggmcmc::ggs_caterpillar(ggm)
  show(catplot)
  ggplot2::ggsave('caterpillar.png')
  
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

