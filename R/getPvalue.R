#' Returning only the p values
#'
#' A group of functions to return only the p values from various of tests.
#' \itemize{
#'   \item p_fish.chi.t: Try first use \code{\link[cjbmisc]{fish.t}} to get p value,\cr
#'   if failed (usually due to too many catergories), then use \code{\link[cjbmisc]{chisq.t}} with popup warning. \cr
#'   If only Fisher's exact test or ChiSQ is wanted, specify it in "p.test" parameter
#'   \item p_aov.t: p-value from \code{\link[stats]{aov}}
#'   \item p_ContDisc: p-value from \code{\link[stats]{t.test}},\code{\link[stats]{wilcox.test}},\code{\link[stats]{kruskal.test}}
#'   \item p_ContCont: p-value from \code{\link[stats]{cor.test}} or from\code{\link[cjbmisc]{p_lm}}
#'   \item p_lm: p-value of the coefficient from the univariate \code{\link[stats]{lm}}
#'   \item p_adjust_mat: adjust a matrix of p-values
#'   \item p_xVsAll: generate a vector of pvalues by contrasting x vs all the rest of the variables in a dataframe.
#'   \item p_dfAll: generate a df of pvalues by contrasting all x.coln vs all y.coln. parameters are the same as in p_xVsAll
#'   \item p_uniCox: generate a list of univariate Cox models and a vector of the minimal pvalue from each model
#'   \item p_feat_subtype: p-value of a set of features according to a subtype in a manner similar to \code{\link[cjbmisc]{plot_feat_subtype}}
#' }
#' @param df numeric dataframe
#' @param v1,v2 character of the column name in the df
#' @param method p_ContDisc: the test method to-be-used, can be t.test, kruskal.test or wilcox.test
#' @param alt passing to alternative in \code{\link[stats]{fisher.test}}, should be one of "two.sided", "greater" or "less". You can specify just the initial letter. Only used in the 2 by 2 case.
#' @param p.test Should be "fisher" or "chi", to use only that method to get p values.
#' @param ws passing to workspace size in \code{\link[stats]{fisher.test}},  integer specifying the size of the workspace used in the network algorithm. In units of 4 bytes.
#' @param df dataframe with rows of samples and columns of features, one column should contain the subtype info.
#' @param subtype column name of the subtype
#' @param feat character vector of the features to be calculated
#' @param cont.test the significance test to be used for continuous variables. Use the name of the r basic tests.
#' @param disc.test the significance test to be used for discrete variables. Should be "fisher","chi" or"both"
#' @param ... for p_feat_subtype, pass to \code{\link[cjbmisc]{p_fish.chi.t}}
#' @seealso \code{\link{chisq.t}}, \code{\link{fish.t}}, \code{\link[stats]{fisher.test}},\code{\link[stats]{chisq.test}},
#' \code{\link[stats]{aov}},\code{\link[stats]{kruskal.test}},\code{\link[stats]{lm}}
#' @name getPvalue
NULL
#' @rdname getPvalue
#' @export
# Only return p values:
p_fish.chi.t <- function(df,v1,v2,alt="two.sided",p.test="both",ws=2e6){
  df <- as.data.frame(df)
  if(p.test=="both"){
    tryCatch(fish.t(df,v1,v2,alt,ws)$p.value,
             error=function(e){
               warning(e)
               warning("Using ChiSQ test instead.")
               chisq.t(df,v1,v2)$p.value})
  }else if (p.test=="fisher"){
    fish.t(df,v1,v2,alt,ws)$p.value
  }else if (p.test=="chi"){
    chisq.t(df,v1,v2)$p.value
  }else stop('p.test must be "fisher" or "chi"')
}
#' @rdname getPvalue
#' @export
p_aov.t <- function(df,v1,v2){
  df <- as.data.frame(df)
  if (plyr::is.discrete(df[,v1,drop=T]))
    f <- as.formula(paste0(v2,"~",v1)) else
      f <- as.formula(paste0(v1,"~",v2))
    aov <- aov(f, df)
    summary(aov)[[1]][["Pr(>F)"]][1]
}
#' @rdname getPvalue
#' @export
p_ContCont <- function(df,v1,v2,method="spearman"){
  # method to choose:
  # spearman,pearson,kendall,lm
  df <- as.data.frame(df)

  if (method=="lm"){
    return(p_lm(df,v1,v2))
  }else{
    return(cor.test(df[,v1],df[,v2],method = method)$p.value)
  }
}
#' @rdname getPvalue
#' @export
p_ContDisc <- function(df,v1,v2,method="kruskal.test"){
  # remove NAs to find out how many groups
  df <- df[!is.na(df[[v1]]) & !is.na(df[[v2]]), ]

  if (plyr::is.discrete(df[[v1]])) {
    df[[v1]] <- as.factor(df[[v1]])
    if (length(unique(na.omit(df[[v1]]))) == 1) {
      warning(paste0("Only one group found for ", v1,
                     ", returning NA for ", method))
      return(NA)
    }
    f <- as.formula(paste0(v2, "~", v1))
  }
  else {
    df[[v2]] <- as.factor(df[[v2]])
    if (length(unique(na.omit(df[[v2]]))) == 1) {
      warning(paste0("Only one group found for ", v1,
                     ", returning NA for ", method))
      return(NA)
    }
    f <- as.formula(paste0(v1, "~", v2))
  }
  get(method)(f, df)$p.value
}
#' @rdname getPvalue
#' @export
p_lm <- function(df,v1,v2){
  x <- df[[v1]]
  y <- df[[v2]]
  pvalue <- coef(summary(lm(y~x)))[2,4]
}
#' @rdname getPvalue
#' @export
p_xVsAll <- function(df,x.coln,y.coln=NULL,
                     num.num.test="spearman",
                     cat.num.test="kruskal.test",
                     cat.cat.test="both"){

  df <- as.data.frame(df)
  # use all columns if y.coln not provided
  if(is.null(y.coln)) y.coln <- setdiff(colnames(df),x.coln)

  df <- df[,c(x.coln,y.coln)]

  feat.cate <- colnames(df)[sapply(df,function(x)!is.numeric(x))]

  # when x is categorical:

  if(x.coln %in% feat.cate){
    vec.p <- sapply(y.coln,function(y){
      if(y %in% feat.cate){
        return(p_fish.chi.t(df,x.coln,y,p.test=cat.cat.test))
      }else{
        return(p_ContDisc(df,x.coln,y,cat.num.test))
      }
    })

  }else{
    vec.p <- sapply(y.coln,function(y){
      if(y %in% feat.cate){
        return(p_ContDisc(df,x.coln,y,cat.num.test))
      }else{
        return(p_ContCont(df,x.coln,y,num.num.test))
      }
    })
  }

  vec.p <- setNames(vec.p,y.coln)

  return(vec.p)
}

#' @rdname getPvalue
#' @export
p_dfAll <- function(df,x.coln,y.coln=NULL,...){
  if(is.null(y.coln)) y.coln <- setdiff(colnames(df),x.coln)

  lst.pval <- lapply(setNames(x.coln,x.coln),function(x)
    p_xVsAll(df,x,y.coln,...) %>% as.data.frame())


  df.pval <- lapply(names(lst.pval), function(n) {
    colnames(lst.pval[[n]]) <- n
    lst.pval[[n]]$Features <- row.names(lst.pval[[n]])
    return(lst.pval[[n]])
  }) %>% Reduce(full_join, .) %>% select(Features, everything())

  return(df.pval)
}

#' @rdname getPvalue
#' @export
p_uniCox <- function(survdf,
                     feat,
                     surv.time = "PFS.days",
                     surv.status = "PFS.status",
                     signif.cutoff=0.05,
                     keep.all.in.barplot=F,
                     plot.surv=T,
                     survp.xlim=c(0,1400),
                     survp.timebreak= 365){
  feat.missing <- setdiff(feat,colnames(survdf))
  if(length(feat.missing)!=0) message(paste0(feat.missing,collapse = ", "),
                                      " are not in the provided dataframe!")
  feat <- intersect(feat,colnames(survdf))

  # fit each feat into a Cox model and return the minimal pvalue.
  lst.fit <- lapply(feat %>% setNames(.,.),
                    function(x){
                      message(x," ",appendLF =F)
                      if(length(unique(na.omit(pull(survdf,!!x))))<2) return(NA)
                      fit <- coxph(as.formula(paste0("Surv(",surv.time,",",surv.status,")~`",x,"`")),
                                   data = survdf)
                    })
  pval <- sapply(lst.fit,function(x)if(!is.na(x))min(summary(x)$coefficients[,5]) else return(NA))
  names(pval) <- gsub("\\.pvalue","",names(pval))

  # Plot the pvalues

  if(keep.all.in.barplot){
    par(mar=c(max(nchar(feat))/3+0.5,5,1,1))
    barplot(-log(sort(pval)),
            las=2,cex.names =0.8,ylab="-log(Pvalue)",
            main="Pvalues of Univariate Cox")
    abline(h=-log(signif.cutoff),lty=2)
  }else if(sum(pval<signif.cutoff,na.rm=T)>0){
    par(mar=c(1+(max(nchar(feat))/2.8),5,1,1))
    barplot(-log(sort(pval[pval<signif.cutoff])),
            las=2,cex.names =0.8,ylab="-log(Pvalue)",
            main="Pvalues of Univariate Cox")
    abline(h=-log(signif.cutoff),lty=2)

  }

  ## plot the Kaplan-Meier plots
  if(plot.surv & (sum(pval<signif.cutoff,na.rm=T)>0) ){
    lst.plot <- lapply(names(sort(pval[pval<signif.cutoff])) %>%
                         setNames(.,.),function(x){
                           p <- cjbmisc::ggsurvpWrap(survdf %>% as.data.frame(),
                                                     x,surv.time,surv.status,
                                                     risk.table = T,xlim = survp.xlim,
                                                     break.time.by = survp.timebreak)$plot+
                             theme(legend.key.size=unit(0.2,"cm"))
                           return(p)
                         }
    )

    p <-ggpubr::ggarrange(plotlist = lst.plot)
    plot(p)

    return(list(fit=lst.fit,pval=pval,plot=p))
  }else

    return(list(fit=lst.fit,pval=pval))

}


#' @rdname getPvalue
#' @export
p_adjust_mat <- function(pvaldf, p.adjust.method = "BH") {

  pvec <- as.vector(pvaldf)
  nai <- is.na(pvec)
  qvec <- rep(NA, length(pvec))
  qvec[!nai] <- p.adjust(pvec[!nai], method = p.adjust.method)
  qmat <- matrix(qvec, nrow = nrow(pvaldf))
  dimnames(qmat) <- dimnames(pvaldf)
  qmat

}
