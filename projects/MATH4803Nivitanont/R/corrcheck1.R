#' @export

#for vignette plot
corrcheck1=function(df, filetype='png', dpi=200){
  corrgram::corrgram(df, upper.panel=NULL)
}