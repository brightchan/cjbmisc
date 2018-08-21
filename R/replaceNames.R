#' Replace the names 
#' 
#' Replace the names from vector or dataframe
#' @param data a vector or a dataframe
#' @param c names before change (vector names or column names)
#' @param n names after change
#' @return the dataframe or vector which using new names (not affect the data)
#' @examples 
#' v <- setNames(c("a","b","c"),c("a","b","c"))
#' replaceNames(v,c("b","c"),c("d","f"))
#' @rdname replaceNames
#' @export
#
replaceNames <- function(data,c,n){
  if(is.vector(data)) names(data)[na.omit(match(c,names(data)))] <- n[c%in%names(data)]
  if(is.data.frame(data)|is.matrix(data)) colnames(data)[na.omit(match(c,colnames(data)))] <- n[c%in%colnames(data)]
  return(data)
}