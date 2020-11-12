#' Plots to explore the correlations between features
#'
#' Plot the correlations focusing on a variable x vs all the rest
#' of the variables.

#'
#' @param plotdf dataframe with rows of samples and columns of features.
#' @param x.coln column name of the x axis of the plot
#' @param y.coln character vector of the column names of features to be plotted as y axis
#' @param cat.num.test the significance test to be used for categorical vs numerical variables. Use the name of the r basic tests (Default "kruskal.test").
#' @param cat.cat.test the significance test to be used for categorical vs categorical variables. Should be "fisher","chi" or "both"(Default)
#' @param num.num.test the significance test to be used for numerical vs numerical variables. Should be "spearman"(Default), "pearson", "kendall", or "lm"(using the pvalue of the independent variable in lm).
#' @param plot.nrow,plot.ncol The number of rows and columns in the combined plot
#' @param plot.signif.only whether to plot only the significant items
#' @param plot.it Whether to plot it out (T/F)
#' @param outpdir If not NULL, save all plots and pvalues (as table) to the outpdir
#' @param ... pass to \code{\link[ggstatsplot]{ggbarstats}} , \code{\link[ggstatsplot]{ggbetweenstats}} ,\code{\link[ggstatsplot]{ggscatterstats}}
#' @return For plot_corr, List of returns from plot_corr_numeric or plot_corr_categorical. For plot_corr_categorical and plot_corr_numeric: List of two: "plot" of ggarrange object which arrange all plot into one, and "pvalues" of a named vector.
#'
#' @name plot_corr
#' @import ggstatsplot
#' @import ggpubr
#' @import dplyr
#'
#' @export
plot_corr_categorical <- function(plotdf,x.coln,y.coln=NULL,
                                  cat.num.test="kruskal.test",
                                  cat.cat.test="both",
                                  plot.it=F,plot.nrow=NULL,plot.ncol=NULL,
                                  signif.cutoff=0.05,
                                  plot.signif.only=F,seed=999,...){

  set.seed(seed)

  # use all columns if y.coln not provided
  if(is.null(y.coln)) y.coln <- setdiff(colnames(plotdf),x.coln)
  # exclude the x.coln in y.coln
  y.coln <- setdiff(y.coln,x.coln)

  feat.cate <- colnames(plotdf)[sapply(plotdf,function(x)!is.numeric(x))]

  # generate all the plots
  message("Plotting: ",appendLF=F)

  ### combine_plots and purrr

  plot.list <- lapply(y.coln,function(x){
    message(x,appendLF=F)
    # Categorical
    if(x %in% feat.cate){
      message(":Cate ",appendLF=F)
      p <- ggbarstats(plotdf[!is.na(plotdf[,x.coln]),],!!x,!!x.coln,
                      title=x,conf.level=(1-signif.cutoff),
                      proportion.test = F,
                      ...)
    } else {
      # Numeric
      message(":Num ",appendLF=F)
      # TODO: add ggwithinstats support for groups with paired samples
      p <- ggbetweenstats(plotdf[!is.na(plotdf[,x.coln]),],!!x.coln,!!x,title=x,
                          type = "np",# usually we should use non-parametric tests
                          pairwise.comparisons = TRUE,conf.level=(1-signif.cutoff),...)
      return(p) }
  })
  plot.list <- setNames(plot.list,y.coln)

  # Arrange the seqeunce of plots by pvalues
  pvalue <- p_xVsAll(plotdf,x.coln,y.coln,cat.num.test=cat.num.test,
                     cat.cat.test=cat.cat.test)
  plot.list <- plot.list[order(pvalue)]


  # Plot only siginificant plots
  if(plot.signif.only){
    if(sum(pvalue<signif.cutoff)>0){
      plot.list <- plot.list[(pvalue<signif.cutoff)[order(pvalue)]]
    } else{
      warning("No significant correlation found.")
      return(list(plot=NULL,pval=pvalue))
    }
  }


  if(is.null(plot.nrow)|is.null(plot.ncol)) outp <- ggpubr::ggarrange(plotlist = plot.list)
  else  outp <- ggpubr::ggarrange(plotlist = plot.list,
                                       nrow = plot.nrow,ncol = plot.ncol,
                                       align = "hv")

  if(plot.it) {
    plot(outp)
  }

  return(list(plot=outp,pval=pvalue))
}

#' @rdname plot_corr
#' @export
plot_corr_numeric <- function(plotdf,x.coln,y.coln,
                              cat.num.test="kruskal.test",
                              num.num.test="spearman",
                              plot.it=F,plot.nrow=NULL,plot.ncol=NULL,
                              signif.cutoff=0.05,
                              plot.signif.only=F,seed=999,...){

  set.seed(seed)

  # use all columns if y.coln not provided
  if(is.null(y.coln)) y.coln <- setdiff(colnames(plotdf),x.coln)
  # exclude the x.coln in y.coln
  y.coln <- setdiff(y.coln,x.coln)


  feat.cate <- colnames(plotdf)[sapply(plotdf,function(x)!is.numeric(x))]

  # generate all the plots
  message("Plotting: ",appendLF=F)
  plot.list <- lapply(y.coln,function(y){
    message(y,appendLF=F)
    # Categorical
    if(y %in% feat.cate){
      message(":Cate ",appendLF=F)
      p <- ggbetweenstats(plotdf[!is.na(plotdf[,y]),],!!y,!!x.coln,title=y,
                          type = "np",# usually we should use non-parametric tests
                          pairwise.comparisons = TRUE,conf.level=(1-signif.cutoff),...)
    } else {
      # Numeric
      message(":Num ",appendLF=F)
      # TODO: add ggwithinstats support for groups with paired samples
      p <- ggscatterstats(plotdf[!is.na(plotdf[,x.coln]),],!!x.coln,!!y,title=y,
                          marginal.type="densigram",...)
      return(p) }
  })
  plot.list <- setNames(plot.list,y.coln)

  # Arrange the seqeunce of plots by pvalues
  pvalue <- p_xVsAll(plotdf,x.coln,y.coln,cat.num.test=cat.num.test,
                     num.num.test=num.num.test)
  plot.list <- plot.list[order(pvalue)]


  # Plot only siginificant plots
  if(plot.signif.only){
    if(sum(pvalue<signif.cutoff)>0){
      plot.list <- plot.list[(pvalue<signif.cutoff)[order(pvalue)]]
    } else{
      warning("No significant correlation found.")
      return(list(plot=NULL,pval=pvalue))
    }
  }


  if(is.null(plot.nrow)|is.null(plot.ncol)) outp <- ggpubr::ggarrange(plotlist = plot.list)
  else  outp <- ggpubr::ggarrange(plotlist = plot.list,
                                       nrow = plot.nrow,ncol = plot.ncol,
                                       align = "hv")

  if(plot.it) {
    plot(outp)
  }

  return(list(plot=outp,pval=pvalue))

}

#' @rdname plot_corr
#' @export
plot_corr <- function(plotdf,x.coln,y.coln=NULL,
                      cat.num.test="kruskal.test",
                      cat.cat.test="both",
                      num.num.test="spearman",
                      plot.it=F,plot.nrow=NULL,plot.ncol=NULL,
                      signif.cutoff=0.05,
                      plot.signif.only=F,seed=999,
                      outpdir=NULL,...){

  # use all columns if y.coln not provided
  if(is.null(y.coln)) y.coln <- setdiff(colnames(plotdf),x.coln)

  feat.cate <- colnames(plotdf)[sapply(plotdf,function(x)!is.numeric(x))]

  lst.out <- lapply(setNames(x.coln,x.coln),
         function(x){
           if(x %in% feat.cate)
             plot_corr_categorical(
               plotdf,x,y.coln=setdiff(y.coln,x),
               cat.num.test="kruskal.test",
               cat.cat.test="both",
               plot.it=plot.it,plot.nrow=plot.nrow,plot.ncol=plot.ncol,
               signif.cutoff=signif.cutoff,
               plot.signif.only=plot.signif.only,seed=seed,...
             )
           else
             plot_corr_numeric(
             plotdf,x,y.coln=setdiff(y.coln,x),
             cat.num.test="kruskal.test",
             num.num.test="spearman",
             plot.it=plot.it,plot.nrow=plot.nrow,plot.ncol=plot.ncol,
             signif.cutoff=signif.cutoff,
             plot.signif.only=plot.signif.only,seed=seed,...
           )
           })

  # plot all plots and tables and put inside outpdir
  if(!is.null(outpdir)){

    # define a basic single plot width and height
    w.plot <- 3.5
    h.plot <- 4

    for(n in names(lst.out)){
      if(!is.null(lst.out[[n]]$plot)){

        # get the num of row and col to decide the plot size
        # following the "arrangeGrob" function used by ggarrange
        n.plots <- length(lst.out[[n]]$plot)
        if (is.null(plot.nrow) && !is.null(plot.ncol)) {
          plot.nrow <- ceiling(n.plots/plot.ncol)
        }else if (is.null(plot.ncol) && !is.null(plot.nrow)) {
          plot.ncol <- ceiling(n.plots/plot.nrow)
        }else if (is.null(plot.nrow) && is.null(plot.ncol)){
          nm <- grDevices::n2mfrow(n.plots)
          plot.nrow = nm[1]
          plot.ncol = nm[2]
        }

        ggsave(paste0(outpdir,"\\",n,".pdf"),lst.out[[n]]$plot,
               width=w.plot*plot.nrow,height = h.plot*plot.ncol)
      }

    }

    # save the pvalues as a table
    lst.pval <- sapply(lst.out,"[[","pval") %>% lapply(as.data.frame)
    df.pval <- lapply(names(lst.pval), function(n){
      colnames(lst.pval[[n]]) <- n
      lst.pval[[n]]$Features <- row.names(lst.pval[[n]])
      return(lst.pval[[n]])
    }) %>% Reduce(full_join,.)

    write.table(df.pval,paste0(outpdir,"\\pvalues.tsv"),sep = "\t",row.names = F)
  }

  return(outp=lst.out)

}
