
#' @export
#' 

diag2=function(mcmc){
  ggs1=ggmcmc::ggs(mcmc)
  f3=ggmcmc::ggs_running(ggs1) + ggthemes::theme_hc()
  f4=ggmcmc::ggs_density(ggs1) + ggthemes::theme_hc()
  ggr2=gridExtra::grid.arrange(f3,f4, nrow=1, ncol=2 )
  summary.mcmc=summary(mcmc)
  return(summary.mcmc)
}
