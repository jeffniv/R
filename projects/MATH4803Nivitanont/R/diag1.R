#' @export
#' 
#' 

diag1=function(mcmc){
  ggs1=ggmcmc::ggs(mcmc)
  # Examine the chains:
  # Convergence diagnostics:
  f1=ggmcmc::ggs_traceplot(ggs1) + ggthemes::theme_hc()
  f2=ggmcmc::ggs_autocorrelation(ggs1) + ggthemes::theme_hc()
  
  ggr1=gridExtra::grid.arrange(f1,f2, nrow=1, ncol=2 )
}
