LASSOReg2 = function( x,
                     y,
                     xnew=NULL,
                     iter=3334,
                     burn=5000,
                     thin=3,
                     alpha=0.05,
                     phi=NULL){
  #check cores for parallel processing
  try( runjags::runjags.options( inits.warning=FALSE , rng.warning=FALSE ) )
  nChainsDefault = 3
  runjagsMethodDefault = "rjags" #default not parallel
  

  fileNameRoot="Jags-MATH4803Nivitanont-LASSOsteps" # For output file names.
  
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
  
  variables=c('beta0', paste0('beta[',1:Nx,']'), 'sigma', 'phi', 'dic' )
  if(!is.null(xnew)){
    variables=append(variables,'ydash')
  }
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
  
  #convert runjags object to mcmc object
  codaSamples = coda::as.mcmc.list( runJagsOut )
  
  #extract DIC from samples
  DICscore=runjags::extract.runjags(x=runJagsOut, what='dic')
  dic=round(do.call(sum,DICscore[1])+do.call(sum,DICscore[2]),2)

  #point and interval estimates
  postsum=summary(codaSamples)
  coefs=postsum$statistics[,1]
  coefs=append(coefs,dic=dic)
  return(coefs)
}


StepLASSO = function( x,
                      y,
                      xnew=NULL,
                      iter=3334,
                      burn=5000,
                      thin=3,
                      alpha=0.05,
                      phi=NULL,
                      DIC=FALSE){
  
  obj = LASSOReg2( x, y, xnew, iter=3334, burn=5000, thin=3, alpha=0.05, phi=NULL)
  
  phi=
  while(phi < 0.1)
  
  return(returnlist)
}


data(clerical)
y=clerical[,1]
x=as.matrix(clerical[,-1])
xnew=c(5000, 75, 900, 200, 650)
