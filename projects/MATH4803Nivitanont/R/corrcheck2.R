#' @export

#for vignette plot
corrcheck2 = function(df, filetype='png', dpi=200){
  s20x::pairs20x(df)
  corrmat=cor(df)
  print(corrmat)
}