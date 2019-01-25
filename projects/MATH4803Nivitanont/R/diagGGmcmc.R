#' Plot MCMC Diagnostics
#'
#' This function takes an \code{mcmc.list} object and outputs diagnostic plots to include trace plots, autocorrelation plots, running mean plots, and posterior density plots. It is especially useful to examine if MCMC chains are well-behaved.
#'
#' @author Jeff Nivitanont, \email{jeffniv@ou.edu}
#'
#' @param mcmc an \code{mcmc.list} object
#'
#' @return Returns a summary of the \code{mcmc.list} object
#'
#' @examples
#' \dontrun{
#' data(codaSamples)
#' mcmc=coda::as.mcmc.list(codaSamples)
#' diagGGmcmc(codaSamples, filetype='png', dpi=200)
#' }
#' @export


diagGGmcmc = function(mcmc,filetype='png', dpi=200){
  #create a ggmcmc object
  ggs1=ggmcmc::ggs(mcmc)
  # Examine the chains:
  # Convergence diagnostics:
  f1=ggmcmc::ggs_traceplot(ggs1) + ggthemes::theme_hc()
  f2=ggmcmc::ggs_autocorrelation(ggs1) + ggthemes::theme_hc()
  f3=ggmcmc::ggs_running(ggs1) + ggthemes::theme_hc()
  #posterior density
  f4=ggmcmc::ggs_density(ggs1) + ggthemes::theme_hc()
  dev.new(noRStudioGD = TRUE)
  ggr1=gridExtra::grid.arrange(f1,f2, nrow=1, ncol=2 )
  ggplot2::ggsave(plot = ggr1, filename = paste0('MCMCDiagnosticsPg1.',filetype), dpi = dpi)
  dev.new(noRStudioGD = TRUE)
  ggr2=gridExtra::grid.arrange(f3,f4, nrow=1, ncol=2 )
  ggplot2::ggsave(plot = ggr2, filename = paste0('MCMCDiagnosticsPg2.',filetype), dpi = dpi)
  summary.mcmc=summary(mcmc)
  print(summary.mcmc)
  return(summary.mcmc)
}
