#' @export


InformedBetaBF2 = function(alpha2,beta2, n, x){
  #function to calculate bayes factor vs uniformed prior
  bayesfactor=function(aa,bb){
    (beta(alpha1, beta1)*beta(x+aa, n-x+bb))/(beta(x+alpha1, n-x+beta1)*beta(aa, bb))
  }
  
  alpha1=1
  beta1=1
  #get bayes factor
  BF=bayesfactor(aa=alpha2, bb=beta2)
  #release to command line
  cat('Bayes Factor (informed vs uniformed)\n',BF)
  
  #plot possible values
  alpha2=0:50+alpha2
  beta2=0:50+beta2
  obj=outer(alpha2, beta2, FUN=bayesfactor)
  filled.contour(x=alpha2, y=beta2,z=obj, color.palette = heat.colors, main = expression(paste("Contour Plot of ",alpha,"vs ",beta) ), xlab=expression(alpha), ylab=expression(beta), key.title = title(main='BF'))
  
  return(BF)
}