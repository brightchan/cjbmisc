#' Returning only the p values 
#' 
#' A group of functions to return only the p values from various of tests.
#' \itemize{
#'   \item discretize: cut a vector by a function or by certain cutoff value using cut.\cr
#'     Return a factor by replacing the original ones with the discrete value regions (<=xxx, xxx-xxx, >xxx).\cr
#'     \emph{NOTE: NAs are automatically removed}
#'   \item trisect: quantile a vector into three equal parts. To be used with discretize: \code{discretize(x,trisect)}
#' }
#' @param df numeric dataframe
#' @param v1,v2 character of the column name in the df
#' @param alt passing to alternative in \code{\link[stats]{fisher.test}}, should be one of "two.sided", "greater" or "less". You can specify just the initial letter. Only used in the 2 by 2 case.
#' @param ws passing to workspace size in \code{\link[stats]{fisher.test}},  integer specifying the size of the workspace used in the network algorithm. In units of 4 bytes.
#' @seealso \code{\link{chisq.t}}, \code{\link{fish.t}}, \code{\link[stats]{fisher.test}},\code{\link[stats]{chisq.test}},
#' \code{\link[stats]{aov}},\code{\link[stats]{kruskal.test}},\code{\link[stats]{lm}} 
#' @name getPvalue
NULL
#' @rdname getPvalue
#' @export
# Only return p values:
p_fish.chi.t <- function(df,v1,v2,alt="two.sided",ws=2e6){
  tryCatch(fish.t(df,v1,v2,alt,ws)$p.value,
           error=function(e){
             warning(e)
             warning("Using ChiSQ test instead.")
             chisq.t(df,v1,v2)$p.value})
}
#' @rdname getPvalue
#' @export
p_aov.t <- function(df,v1,v2){
  if (plyr::is.discrete(df[,v1]))
    f <- as.formula(paste0(v2,"~",v1)) else
      f <- as.formula(paste0(v1,"~",v2))
    aov <- aov(f, df)
    summary(aov)[[1]][["Pr(>F)"]][1]
}
#' @rdname getPvalue
#' @export
p_krus <- function(df,v1,v2){
  if (plyr::is.discrete(df[,v1]))
    f <- as.formula(paste0(v2,"~",v1)) else
      f <- as.formula(paste0(v1,"~",v2))
    kruskal.test(f, df)$p.value
}
#' @rdname getPvalue
#' @export
p_lm <- function(df,v1,v2){
  x <- df[,v1]
  y <- df[,v2]
  pvalue <- coef(summary(lm(y~x)))[2,4]
}
