#' Discretize a continuous vector
#'
#' A group of functions to divide continous numeric vectors into discret parts.
#' Similar to \code{\link[stats]{quantile}} then \code{\link[base]{cut}} but easier to use.
#' \itemize{
#'   \item discretize: cut a vector by a function or by certain cutoff value using cut.\cr
#'     Return a factor by replacing the original ones with the discrete value regions (<=xxx, xxx-xxx, >xxx).\cr
#'     \emph{NOTE: NAs are automatically removed}
#'   \item trisect: quantile a vector into three equal parts. To be used with discretize: \code{discretize(x,trisect)}
#' }
#' @param x numeric vector
#' @param cutoff numeric value (single or vector) as cutting point(s) (Don't include min or maximum),
#'  or a function generating one cutting points (median, mean, etc) or multiple cutting points which include (quantile).
#' @param label change the default output region labels to customized labels
#' @param include.lowest see \code{\link[base]{cut}} parameter \emph{include.lowest}
#' @param summary when TRUE, message a summary of the count of each discrete items.
#' @examples
#' discretize(rnorm(20),trisect,summary=TRUE)
#' discretize(rnorm(20),mean,label=c("low","high"),summary=TRUE)
#' @name discretizes
NULL
#' @rdname discretizes
#' @export
#cut continuous vector to discrete with gt/lt number label.
discretize <- function(x,cutoff,label=NULL,include.lowest=TRUE,summary=FALSE){
  if (is.function(cutoff)) {
    if(length(cutoff(na.omit(x)))==1) {
      cutoff <- c(min(x,na.rm=TRUE),cutoff(na.omit(x)),max(x,na.rm=TRUE))
    } else cutoff <- cutoff(na.omit(x))
  } else cutoff <- c(min(x,na.rm=TRUE),cutoff,max(x,na.rm=TRUE))
  #add label for the discretized vector
  cutlabel <- signif(cutoff,3)
    if(is.null(label)){
      label <- c(paste0("<=",cutlabel[2]))
      if(length(cutoff)>3) for(i in 2:(length(cutoff)-2)) label <- c(label,paste0(cutlabel[i],"-",cutlabel[i+1]))
      label <- c(label,paste0(">",cutlabel[length(cutoff)-1]))
    }
  out <- cut(x,cutoff,label,include.lowest=include.lowest)
  if(summary) print(table(out))
  return(out)
}
#' @rdname discretizes
#' @export
trisect <- function(x,na.rm=FALSE){
  quantile(x,c(0,1/3,2/3,1),na.rm=na.rm)
}
