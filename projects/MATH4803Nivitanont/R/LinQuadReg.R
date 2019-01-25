#' Linear vs Quadratic Regression Comparison for One Predictor Variable
#'
#' This function takes in a vector of predictor variables \code{x} and a vector of data \code{y} and creates a quadratic model and a simple linear model using MCMC sampling. For the Bayesian inference, this method uses uninformative prior distributions. The point and interval estimates are taken from the posterior distribution from the samples. The models are then plotted over the data and DIC is calculated to determine goodness of fit
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
#' @examples
#' \dontrun{
#' data("electric")
#' y=electric[,3]
#' x=electric[,1]
#' LinQuadReg(  x,
#'              y,
#'              xnew=NULL,
#'              iter=3334,
#'              burn=5000,
#'              thin=3,
#'              alpha=0.05)
#' }
#' @export



LinQuadReg=function(x,
           y,
           xnew=NULL,
           iter=3334,
           burn=5000,
           thin=3,
           alpha=0.05){

  try( runjags::runjags.options( inits.warning=FALSE , rng.warning=FALSE ) )
  nChainsDefault = 3
  runjagsMethodDefault = "rjags" #default not parallel for DIC computation


  fileNameRoot="Jags-MATH4803Nivitanont-LinQuadReg" # For output file names.

  # Put the information into a list.
  Ntotal = length(y)
  dataList = list(
    y = y ,
    x = x,
    Ntotal = Ntotal
  )
  dataList_l = list(
    y = y,
    x = x,
    Ntotal = Ntotal
  )
  if(!is.null(xnew)){
    dataList=c(dataList, xdash=xnew)
  }
  # Define the model:
  modelString = paste0("
  # center the data
  data {
  ym <- mean(y)
  xm <- mean(x)
  for ( i in 1:Ntotal ) {
  zy[i] <- ( y[i] - ym )
  zy_l[i] <- zy[i]
  zx[i] <- ( x[i] - xm)
  }
  }
  model{
  for( i in 1 : Ntotal ) {
  zy[i] ~ dnorm(mu[i], tau)
  mu[i] <- zbeta0 + zbeta1 * zx[i] + beta2 * zx[i]^2

  zy_l[i] ~ dnorm(mu_l[i], 1/(sigma_l^2))
  mu_l[i] <- zbeta0_l + beta1_l * zx[i]
  }

  zbeta0 ~ dnorm(0.0, 1.0E-6)
  zbeta1 ~ dnorm(0.0, 1.0E-6)
  beta2 ~ dnorm(0.0, 1.0E-6)
  sigma ~ dunif(0, 100)
  tau <- pow(sigma,  -2)

  zbeta0_l ~ dnorm(0.0, 1.0E-6)
  beta1_l ~ dnorm(0.0, 1.0E-6)
  sigma_l ~ dunif(0, 100)


  #transform data from scaled
  beta0 <- zbeta0 + ym - zbeta1*xm + beta2*xm^2
  beta1 <- zbeta1 - beta2*2*xm

  beta0_l <- zbeta0_l + ym - beta1_l*xm
  ",ifelse(!is.null(xnew),"
  #given x value
  ydash ~ dnorm(mudash, tau)
  mudash <- beta0 + beta1 * xdash",""),"

  }" # close quote for modelStringlose quote for modelString
  )

  modelString_l = "
  # center the data
  data {
  ym <- mean(y)
  xm <- mean(x)
  for ( i in 1:Ntotal ) {
  zy_l[i] <- ( y[i] - ym )
  zx[i] <- ( x[i] - xm)
  }
  }
  model{
  for( i in 1 : Ntotal ) {
  zy_l[i] ~ dnorm(mu_l[i], 1/(sigma_l^2))
  mu_l[i] <- zbeta0_l + beta1_l * zx[i]
  }
  zbeta0_l ~ dnorm(0.0, 1.0E-6)
  beta1_l ~ dnorm(0.0, 1.0E-6)
  sigma_l ~ dunif(0, 100)

  #transform data from scaled
  beta0_l <- zbeta0_l + ym - beta1_l*xm
  }" # close quote for modelStringlose quote for modelString

  writeLines( modelString , con="TEMPass4modelQuad.txt" )

  #parameters to monitor
  parameters = c('beta0', 'beta1', 'beta2','sigma', 'dic')
  if(!is.null(xnew)){
    parameters=append(parameters,'ydash')
  }
  parameters_l = c('beta0_l', 'beta1_l', 'sigma_l','dic')

  #start timer
  start.time <- proc.time()
  # Run the chains, output is runjags object:
  runJagsOut = runjags::run.jags(model = modelString,
                        data = dataList,
                        n.chains = 3,
                        monitor = parameters,
                        burnin = burn,
                        adapt = 1000,
                        sample = iter,
                        thin = thin,
                        method = runjagsMethodDefault,
                        summarise = FALSE)
  runJagsOut_l = runjags::run.jags(model = modelString_l,
                          data = dataList_l,
                          n.chains = 3,
                          monitor = parameters_l,
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
  #convert runjags object to mcmc object
  codaSamples_l = coda::as.mcmc.list( runJagsOut_l )
  save( codaSamples_l , file=paste0(fileNameRoot,"_l_Mcmc.Rdata") )

  DICscore=runjags::extract.runjags(x=runJagsOut, what='dic')
  DIC_lscore=runjags::extract.runjags(x=runJagsOut_l, what='dic')

  #plot diagnostics
  HDI=1-alpha
  diagGGmcmc(codaSamples)
  diagGGmcmc(codaSamples_l)

  #combine sim data
  simdf=c(codaSamples, codaSamples_l)
  #point and interval estimates
  postsum=summary(codaSamples)
  coefs=postsum$statistics[,1]
  HPD=coda::HPDinterval(codaSamples, prob=HDI)

  postsum_l=summary(codaSamples_l)
  coefs_l=postsum_l$statistics[,1]
  HPD_l=coda::HPDinterval(codaSamples_l, prob=HDI)

  dic=round(do.call(sum,DICscore[1])+do.call(sum,DICscore[2]),2)
  dic_l=round(do.call(sum,DIC_lscore[1])+do.call(sum,DIC_lscore[2]),2)
  #plot models for comparison
  xmin=min(x);  xmax=max(x)
  ymin=min(y);  ymax=max(y)
  xvals=seq(xmin, xmax, length=100)
  quad=function(x){
    coefs[1] + coefs[2]*x + coefs[3]*x^2
  }
  linear=function(x){
    coefs_l[1] + coefs_l[2]*x
  }
  ggdf <- data.frame(
    x = rep(xvals,2),
    y = c(quad(xvals), linear(xvals)),
    Model = rep(c("Quadratic","Linear"),
                each = 100)
  )
  ggdf2 = data.frame(
    x=x,
    y=y,
    Model= rep('Original', length(x)))
  dev.new(noRStudioGD = TRUE)
  p <- ggplot2::ggplot(ggdf, ggplot2::aes(x=x, y=y, Model=Model))
  p = p + ggplot2::geom_line(ggplot2::aes(colour=Model)) + ggplot2::geom_point(data = ggdf2, colour='black' )
  p = p + ggplot2::labs(title=paste("Linear vs Quadratic Fit  (alpha=",alpha,')'))
  p = p + ggplot2::annotate('text',
                   label=paste('Quadratic:',round(coefs[1],4), '+', round(coefs[2],4),'*x + ', round(coefs[3],4),'*x^2\n DIC:', round(dic, 4)),
                   x = mean(x),
                   y=ymax*1.05)
  p = p + ggplot2::annotate('text',label=paste('Linear:',round(coefs_l[1],4), '+', round(coefs_l[2],4),'*x\n DIC:', round(dic_l, 4)),
                   x = mean(x),
                   y=ymin*1.05)
  show(p)
  ggplot2::ggsave(filename=paste(fileNameRoot,'-LinVsQuad.png'), device = 'png', dpi=200)
  #command line release
  cat('Point Estimates\n')
  cat('Quadratic:\n')
  print(coefs)
  cat('Linear:\n')
  print(coefs_l)
  cat('Interval Estimates\n')
  cat('Quadratic:\n')
  print(HPD[1])
  cat('Linear:\n')
  print(HPD_l[1])
  cat('DIC\n')
  cat('Quadratic:\n')
  print(DICscore)
  cat('Linear:\n')
  print(DIC_lscore)

  returnedlist=list(sim=codaSamples, y=y, x=x,
                    coefs=append(coefs,coefs_l),
                    summary=append(postsum,postsum_l),
                    HPD=append(HPD,HPD_l),
                    DIC=c(DICscore, DIC_lscore))
  return(returnedlist)
}
