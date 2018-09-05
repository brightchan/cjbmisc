#' Returning only the p values
#'
#' A group of functions to return only the p values from various of tests.
#' \itemize{
#'   \item p_fish.chi.t: Try first use \code{\link[cjbmisc]{fish.t}} to get p value,\cr
#'   if failed (usually due to too many catergories), then use \code{\link[cjbmisc]{chisq.t}} with popup warning.
#'   \item p_aov.t: p-value from \code{\link[stats]{aov}}
#'   \item p_krus: p-value from \code{\link[stats]{kruskal.test}}
#'   \item p_lm: p-value of the coefficient from the univariate \code{\link[stats]{lm}}
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
