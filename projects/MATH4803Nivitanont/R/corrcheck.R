#' Correlation check
#'
#' This function takes in a dataframe object and creates a correlation heatmap, a correlation description matrix, and releases the correlation matrix to the command line. It is useful for screening data sets for possible multicollinearity issues.
#'
#' @author Jeff Nivitanont, \email{jeffniv@ou.edu}
#'
#' @param df A data frame that ontains the data for analysis.
#'
#' @return Returns a correlation matrix.
#'
#' @examples
#' \dontrun{
#' data(clerical)
#' corrcheck(clerical)
#' ##               Hours         Mail       Gifts     Charge     Returns      Checks
#' ##Hours    1.000000000 -0.007650103  0.29281923 0.46151908  0.08479822  0.58731901
#' ##Mail    -0.007650103  1.000000000  0.01128202 0.05480359 -0.04311752 -0.27658574
#' ##Gifts    0.292819235  0.011282017  1.00000000 0.24521511  0.03686148 -0.01588972
#' ##Charge   0.461519082  0.054803588  0.24521511 1.00000000  0.47780716  0.50899367
#' ##Returns  0.084798217 -0.043117518  0.03686148 0.47780716  1.00000000  0.44280516
#' ##Checks   0.587319010 -0.276585736 -0.01588972 0.50899367  0.44280516  1.00000000
#' ##Hours         Mail       Gifts     Charge     Returns      Checks
#' ##Hours    1.000000000 -0.007650103  0.29281923 0.46151908  0.08479822  0.58731901
#' ##Mail    -0.007650103  1.000000000  0.01128202 0.05480359 -0.04311752 -0.27658574
#' ##Gifts    0.292819235  0.011282017  1.00000000 0.24521511  0.03686148 -0.01588972
#' ##Charge   0.461519082  0.054803588  0.24521511 1.00000000  0.47780716  0.50899367
#' ##Returns  0.084798217 -0.043117518  0.03686148 0.47780716  1.00000000  0.44280516
#' ##Checks   0.587319010 -0.276585736 -0.01588972 0.50899367  0.44280516  1.00000000
#' }
#' @export



corrcheck = function(df, filetype='png', dpi=200){
  dev.new(noRStudioGD = TRUE)
  corrgram::corrgram(df, upper.panel=NULL)
  dev.new(noRStudioGD = TRUE)
  s20x::pairs20x(df)
  corrmat=cor(df)
  print(corrmat)
  return(corrmatrix=corrmat)
}
