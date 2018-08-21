#' Subset rows from a dataframe by various criteria
#' 
#' A group of functions to subset rows from a dataframe 
#' \itemize{
#'   \item topvar: subset top p% rows by variance (or other methods like mean/median) 
#'   \item topvar2: subset rows with at least thres variance.
#'   \item topvar3: subset top n rows by variance.
#'   \item \emph{NOTE: for topvars, rows with zero sum are automatically removed}
#'   \item subset_nonzero: subset rows with at least p% non-zero numbers. (For filtering non-expressed genes)
#'   \item subset_nonzero2: subset rows with at least n non-zero numbers
#'   \item subset_rowname: Subset one dataframe according to another, by rownames. Return df has same order as the reference df
#' }
#' @param df numeric dataframe
#' @param v method to calculate the variance (\code{\link{cv}}, \code{\link[stats]{mad}}, \code{\link[stats]{mean}},etc)
#' @param p percentage, 0-100
#' @param thres a threshold of variance for the cutoff
#' @param n integer
#' @return dataframe
#' @name subsetRows
NULL
#' @rdname subsetRows
#' @export
# Pick rows with top variance by methods specified in v and perc. in p.
topvar <- function(df,p,v=cv){
  df <- df[rowSums(df)>0,] # DELETE Rows with zero sum. Quantile would vary
  criterion=apply(df,1,v)
  thres <- quantile(criterion,(1-p/100))
  subdata <- df[criterion>=thres,]
  return(subdata)
}
#' @rdname subsetRows
#' @export
# pick rows with at least thres variance.
topvar2 <- function(df,thres,v=cv){
  df <- df[rowSums(df)>0,] # DELETE Rows with zero sum
  criterion=apply(df,1,v)
  subdata <- df[criterion>=thres,]
  return(subdata)
}
#' @rdname subsetRows
#' @export
# pick rows with top n variance.
topvar3 <- function(df,n,v=cv){
  df <- df[rowSums(df)>0,] # DELETE Rows with zero sum
  subdata <- df[order(apply(df,1,v),decreasing = T)[1:n],]
  return(subdata)
}
#' @rdname subsetRows
#' @export
# Pick rows with at least p% non-zero numbers
subset_nonzero <- function(df,p){
  nonz <- rowSums(df!=0)
  thres <- ncol(df)*p/100
  df <- df[nonz>=thres,]
}
#' @rdname subsetRows
#' @export
# Pick rows with at least n non-zero numbers
subset_nonzero2 <- function(df,n){
  nonz <- rowSums(df!=0)
  df <- df[nonz>=n,]
}
#' @rdname subsetRows
#' @export
# Get subset of x according to y by rownames and order by rownames of y. Input data.frame
subset_rowname <- function(x,y,order=T) { 
  xrown <- row.names(x)
  yrown <- row.names(y)
  subx <- x[(xrown %in% yrown),]
  message("Rows:", as.character(nrow(x))," -> ",as.character(nrow(subx)))
  if(order) subx <- subx[yrown[yrown %in% xrown],]
  return(subx)
}
