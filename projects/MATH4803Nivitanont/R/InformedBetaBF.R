#' Bayes Factor for an Informed Beta Prior Distribution
#'
#' This function takes parameters for an informed Beta prior distribution and compares it to the uninformed prior \code{beta(1,1)} and releases a command line output of the Bayes Factor.
#'
#' @author Jeff Nivitanont, \email{jeffniv@ou.edu}
#'
#' @param alpha2 The alpha parameter for the informed Beta prior distribution
#'
#' @param beta2 The beta parameter for the informed Beta prior distribution
#'
#' @param n The total number of trials
#'
#' @param x The total number of successes
#'
#' @return Returns \eqn{P(M_2)/P(M_1)}, where \eqn{M_1} is the uninformed prior.
#' @examples
#' InformedBetaBF(alpha2=2, beta2=5, n=10, x=4)
#' ##Bayes Factor (informed vs uniformed)
#' ##1.442308
#' ##contour plots not shown
#' @export

InformedBetaBF = function(alpha2,beta2, n, x){
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
  dev.new(noRStudioGD = TRUE)
  filled.contour(x=alpha2, y=beta2,z=obj, color.palette = heat.colors, main = expression(paste("Contour Plot of ",alpha,"vs ",beta) ), xlab=expression(alpha), ylab=expression(beta), key.title = title(main='BF'))

  return(BF)
}
