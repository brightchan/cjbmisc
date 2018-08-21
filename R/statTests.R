#' Simplified stat tests
#'
#' A group of functions to quickly run certain statistical tests from a dataframe.
#' \itemize{
#'   \item fish.t: Fisher's exact test with columns in dataframe.
#'   \item chisq.t Chi Square Test with columns in dataframe
#'   \item one.v.all.t: One versus all test of contingency table, return a named p value vector
#' }
#' @param df numeric dataframe
#' @param v1,v2 character of the column name in the df
#' @param alt passing to alternative in \code{\link[stats]{fisher.test}}, should be one of "two.sided", "greater" or "less". You can specify just the initial letter. Only used in the 2 by 2 case.
#' @param ws passing to workspace size in \code{\link[stats]{fisher.test}},  integer specifying the size of the workspace used in the network algorithm. In units of 4 bytes.
#' @param test the test to be used (fisher.test, etc). the ... is passed to the test function
#' @seealso  \code{\link[stats]{fisher.test}},\code{\link[stats]{chisq.test}},
#' @name stattests
NULL
#' @rdname stattests
#' @export
# Fisher's exact test with columns in dataframe
fish.t <- function(df,v1,v2,alt="two.sided",ws=2e7){#v1,v2:name of col
  f <- as.formula(paste0("~",v1,"+",v2))
  cont.table <- xtabs(f,data=df)
  fisher.test(cont.table,alternative=alt,workspace=ws)
}
#' @rdname stattests
#' @export
# Chi Square Test with columns in dataframe
chisq.t <- function(df,v1,v2){#v1,v2:name of col
  f <- as.formula(paste0("~",v1,"+",v2))
  cont.table <- xtabs(f,data=df)
  chisq.test(cont.table)
}
#' @rdname stattests
#' @export
# One versus all test of contingency table
one.v.all.t <- function(df,test=fisher.test,...){
  sapply(setNames(colnames(df),colnames(df)),
         function(s){
    tmp <- data.frame(df[,s],rowSums(df)-df[,s])
    pvalue <- test(tmp,...)$p.value
    if (pvalue<0.01) message(s," p=",pvalue,"**") 
    else if (pvalue<0.05) message(s," p=",pvalue,"*")
    else message(s," p=",pvalue)
  })
}
