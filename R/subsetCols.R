#' Subset columns from a dataframe by various criteria
#' 
#' A group of functions to subset columns from a dataframe 
#' \itemize{
#'   \item del_dupcol: remove duplicated columns by name.
#'   \item del_alldupcol: remove columns that are not unique by name
#'   \item subset_colname: Subset one dataframe according to another, by colnames. Return df has same order as the reference df
#'   \item subcol_tcgaTnN: Subset tumor/normal from tcgacount
#'   \item subset_tcgasam: Subset one dataframe according to another, by colnames of patient name. Will truncate long names to the patient barcode for you.
#' }
#' @param df,x,y dataframe
#' @param tum TRUE for tumor and FALSE for normal
#' @return dataframe
#' @name subsetCols
NULL
#' @rdname subsetCols
#' @export
# remove duplicated columns from DF
del_dupcol <- function(df){
  coln <- colnames(df)
  df <- df[,coln[!duplicated(coln)]]
}
#' @rdname subsetCols
#' @export
# remove ALL duplicated columns from DF
del_alldupcol <- function(df){
  coln <- colnames(df)
  df <- df[,!(coln %in% coln[duplicated(coln)])]
}
#' @rdname subsetCols
#' @export
# Get subset of x according to y by colnames and order by colnames of y. Input data.frame
subset_colname <- function(x,y,order=T) { 
  xcoln <- colnames(x)
  ycoln <- colnames(y)
  subx <- x[,(xcoln %in% ycoln)]
  message("Cols:", as.character(ncol(x))," -> ",as.character(ncol(subx)))
  if(order) subx <- subx[,ycoln[ycoln %in% xcoln]]
  return(subx)
}
#' @rdname subsetCols
#' @export
# Subset tumor/normal from tcgacount
subcol_tcgaTnN <- function (df,tum=T) { #tum=T for tumor and tum=F for normal
  sname <- sapply(strsplit(colnames(df),"[.]|-"),"[[",4)
  sname <- as.numeric(gsub("\\D","",sname))
  if (tum){
    df <- df[,(sname<10)]
    return(df)
  }
  else {
    df  <- df[,(sname>9)]
    return(df)
  }
}
#' @rdname subsetCols
#' @export
# Subset x according to y, by patient name in colname only
subset_tcgasam <- function(x,y) { 
  colnames(x) <- substr(colnames(x),1,12) #truncate to patient name
  x <- del_dupcol(x)
  colnames(y) <- substr(colnames(y),1,12) 
  y <- del_dupcol(y)
  x <- subset_colname(x,y)
  return(x) 
}
