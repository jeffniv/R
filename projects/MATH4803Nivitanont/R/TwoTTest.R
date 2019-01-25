#' Two-Sample Bayesian T-Test
#'
#' This function takes in two vectors of data \code{y1,y2} and performs a t-test using MCMC sampling. For the Bayesian inference, this method uses uninformative prior distributions. The point and interval estimates are taken from the posterior distribution from the samples. The function also generates an interactive 3D contour plot of the two sample densities.
#'
#' @author Jeff Nivitanont, \email{jeffniv@ou.edu}
#'
#' @param y1,y2 A vector of data
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
#' @examples
#' \dontrun{
#' data("reading")
#' y1=reading$Grade[reading$Method == 'c']
#' y2=reading$Grade[reading$Method == 't']
#' testobj=two.t.test(y1,
#'                    y2,
#'                    iter=3334,
#'                    burn=5000,
#'                    thin=3,
#'                    alpha=0.05,
#'                    parallel=TRUE)
#' }
#' @export

two.t.test = function(y1,
                      y2,
                      iter=3334,
                      burn=5000,
                      thin=3,
                      alpha=0.05,
                      parallel=FALSE){
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
  #temp file save name
  fileNameRoot = 'Jags-MATH4803Nivitanont-two.t.test'
  #JAGS model
  modelString ="
  model{
  for( i in 1 : n1 ) {
  y1[i] ~ dt(mu1, tau1, nu)
  }
  for( j in 1 : n2 ) {
  y2[j] ~ dt(mu2, tau2, nu)
  }
  mu1 ~ dnorm(0.0, 1.0E-6)
  mu2 ~ dnorm(0.0, 1.0E-6)
  nuMinus1 ~ dexp(1/29)
  sigma1 ~ dunif(0, 100)
  sigma2 ~ dunif(0, 100)
  muDiff <- mu1 - mu2
  nu <- nuMinus1 + 1
  sigDiff <- sigma1 - sigma2
  tau1 <- pow(sigma1,  -2)
  tau2 <- pow(sigma2,  -2)
  }"

  #data list for JAGS
  dataList = list(
    y1=y1,
    y2=y2,
    n1=length(y1),
    n2=length(y2)
  )

  #write model to txt
  writeLines( modelString , con="two.t.testmodel.txt" )

  #observed params
  variables=c('mu1', 'mu2', 'muDiff', 'sigma1', 'sigma2', 'sigDiff', 'nu')

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

  #plot diagnostics
  diagGGmcmc(codaSamples)

  #convert mcmc to matrix for plotting
  xx=as.matrix(codaSamples[,1])
  yy=as.matrix(codaSamples[,2])
  #kernel density 2d
  sim.kde = MASS::kde2d(xx, yy, n=300)
  dev.new(noRStudioGD = TRUE)
  filled.contour(sim.kde, xlab = expression(mu[1]), ylab = expression(mu[2]), main=expression(paste(mu[1],' vs ', mu[2])))
  savePlot(filename='kerneldensitycompare.png',type = 'png')
  rgl::persp3d(sim.kde, col='lightblue', box=FALSE, xlab='mu1', ylab='mu2', zlab='')
  rgl::snapshot3d("kerneldensity3d.png", fmt='png')

  #compile y values for return
  df=data.frame(matrix(nc=2, nr=length(y1)))
  df[,1] = y1
  df[,2] = y2
  names(df) = c('y1', 'y2')

  #point and interval estimates
  postsum=summary(codaSamples)
  coefs=postsum$statistics[,1]
  HPD=coda::HPDinterval(codaSamples, prob=1-alpha)
  #command line release
  cat('Point Estimates\n')
  print(coefs)
  cat('Interval Estimates\n')
  print(HPD[1])

  returnlist = list(sim=codaSamples, coefs=coefs, summary=postsum, HPD=HPD, y=df)

  return(returnlist)
}
