#' Plots to explore the correlations between features
#'
#' Plot the correlations focusing on a variable x vs all the rest
#' of the variables.
#' The workflow is: 1. remove small groups if "min.group.size" is defined;
#' 2. calculate the p values for all pairs of variables
#' 3. select the ones that pass pvalue threshold for plotting.
#'   Pvalues are by default non-parametric. Can choose if p.adjust should be used.
#' 4. calculate padj and save the plots and pvalue tables.

#'
#' @param plotdf dataframe with rows of samples and columns of features.
#' @param x.coln column name of the x axis of the plot
#' @param y.coln character vector of the column names of features to be plotted as y axis
#' @param cat.num.test the significance test to be used for categorical vs numerical variables. Use the name of the r basic tests (Default "kruskal.test").
#' @param cat.cat.test the significance test to be used for categorical vs categorical variables. Should be "fisher","chi" or "both"(Default)
#' @param num.num.test the significance test to be used for numerical vs numerical variables. Should be "spearman"(Default), "pearson", "kendall", or "lm"(using the pvalue of the independent variable in lm).
#' @param plot.nrow,plot.ncol The number of rows and columns in the combined plot
#' @param padj.method,padj.method.each,padj.method.all pvalue adjustment method. Should follow \code{\link[ggstatsplot]{ggbetweenstats}}.
#' In plot_corr_one, if padj.method is specified, pvalue adjustment will be done and only those pass the padj threshold will be plotted.
#' In plot_corr, padj.method.each will be passed to plot_corr_one while padj.method.all will be used to genearate the final padj table, adjusting for all pvalues.
#' @param plot.stattest Pass to the "type" parameter in \code{\link[ggstatsplot]{ggbarstats}} , \code{\link[ggstatsplot]{ggbetweenstats}} ,\code{\link[ggstatsplot]{ggscatterstats}}, defining the stats test to be used. Default "np" is non-parametric.
#' @param plot.signif.only whether to plot only the significant items
#' @param plot.max maximum how many plots to be plotted. If set to NULL then plot all.
#' @param min.group.size.x for categorical x, remove groups that are smaller than this number
#' @param min.group.size.y for categorical y, remove groups that are smaller than this number
#' @param plot.it Whether to plot it out (T/F)
#' @param outpdir If not NULL, save all plots and pvalues (as table) to the outpdir
#' @param plot.w,plot.h Width and height of each individual plot
#' @param fn.suffix filename suffix
#' @param ... pass to \code{\link[ggstatsplot]{ggbarstats}} , \code{\link[ggstatsplot]{ggbetweenstats}} ,\code{\link[ggstatsplot]{ggscatterstats}}
#' @return For plot_corr, List of returns from plot_corr_one. For plot_corr_one: List of two: "plot" of ggarrange object which arrange all plot into one, and "pvalues" of a named vector.
#'
#' @name plot_corr
#' @import ggstatsplot
#' @import ggpubr
#' @import dplyr
#'
#' @export
plot_corr_one <- function(plotdf,x.coln,y.coln,
                          cat.num.test="kruskal.test",
                          num.num.test="spearman",
                          plot.it=F,plot.nrow=NULL,plot.ncol=NULL,
                          signif.cutoff=0.05,
                          plot.stattest="np",
                          p.adj.method=NULL,
                          plot.signif.only=F,
                          plot.max=40,
                          min.group.size.x=3,
                          min.group.size.y=3,
                          seed=999,
                          outpdir=NULL,
                          plot.w=7,plot.h=7.5,
                          fn.suffix="",...){


  set.seed(seed)

  # use all columns if y.coln not provided
  if(is.null(y.coln)) y.coln <- setdiff(colnames(plotdf),x.coln)
  # exclude the x.coln in y.coln
  y.coln <- setdiff(y.coln,x.coln)


  # check which are categorical columns
  feat.cate <- colnames(plotdf)[sapply(plotdf,function(x)!is.numeric(x))]


  # remove small groups from categorical features by setting them to NA
  # Keep a Note in the output file suffix
  if(!is.null(min.group.size.x)&x.coln%in%feat.cate){
    fn.suffix <- paste0(fn.suffix,"_grp.x.gt",min.group.size.x)
    plotdf <- plotdf %>%
      mutate_at(vars(!!x.coln),~ifelse(table(.)[.]<min.group.size.x,NA,.))
  }
  if(!is.null(min.group.size.y)&any(y.coln%in%feat.cate)){
    fn.suffix <- paste0(fn.suffix,"_grp.y.gt",min.group.size.y)
    plotdf <- plotdf %>%
      mutate_at(vars(one_of(intersect(y.coln,feat.cate))),
                ~ifelse(table(.)[.]<min.group.size.y,NA,.))
  }


  # Prompt WARNING and stop if x.coln only has one group
  if(x.coln%in%feat.cate&(length(unique(na.omit(as.vector(plotdf[[x.coln]]))))<2)){
    warning("Catergorical variable ", x.coln, " only have one effective group. Stop plotting.")
    return(list(plot=NULL,pval=NULL))
  }

  # remove categorical y.coln if only one group remains
  if(any(y.coln%in%feat.cate)){
    count.feat <- plotdf %>% select(one_of(intersect(y.coln,feat.cate))) %>%
      apply(2,function(x)length(unique(na.omit(x))))
    y.coln <- setdiff(y.coln,names(count.feat[count.feat<2]))
  }
  # Prompt WARNING and stop if no more y.coln left
  if(length(y.coln)==0){
    warning("No more effective y column left. Stop plotting.")
    return(list(plot=NULL,pval=NULL))
  }



  message("*****Plotting: ",x.coln," vs ",appendLF=F)

  ### calculate the pvalues
  pvalue.raw <- p_xVsAll(plotdf,x.coln,y.coln,cat.num.test=cat.num.test,
                     num.num.test=num.num.test)

  if(!is.null(p.adj.method)) pvalue <- p.adjust(pvalue.raw,p.adj.method)
  else pvalue <- pvalue.raw


  ### generate all the plots


  if(plot.signif.only){
    # Plot only siginificant plots
    # by ignoring all non-significant ones
    if(sum(pvalue<signif.cutoff,na.rm=T)>0){
      y.coln.plot <- na.omit(names(pvalue[pvalue<signif.cutoff]))
    } else{
      warning("No significant correlation found.")
      return(list(plot=NULL,pval=pvalue))
    }
  }


  # Arrange the seqeunce of plots by pvalues
  y.coln.plot <- y.coln.plot[order(pvalue[y.coln.plot])]

  # remove NA values
  y.coln.plot <- na.omit(y.coln.plot)

  # plot only the max number of plots if number of significant plot is more than plot.max
  if(!is.null(plot.max)&(length(y.coln.plot)>plot.max)) y.coln.plot <- y.coln.plot[1:plot.max]


  # use different plots for different combinations of catergorical/numerical values
  if(x.coln %in% feat.cate){
    plot.list <- lapply(y.coln.plot,function(x){
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
                            type = plot.stattest,# usually we should use non-parametric tests
                            pairwise.comparisons = TRUE,conf.level=(1-signif.cutoff),...)
        return(p) }
    })
  } else{
    plot.list <- lapply(y.coln.plot,function(y){
      message(y,appendLF=F)
      # Categorical
      if(y %in% feat.cate){
        message(":Cate ",appendLF=F)
        p <- ggbetweenstats(plotdf[!is.na(plotdf[,y]),],!!y,!!x.coln,title=y,
                            type = plot.stattest,# usually we should use non-parametric tests
                            pairwise.comparisons = TRUE,conf.level=(1-signif.cutoff),...)
      } else {
        # Numeric
        message(":Num ",appendLF=F)
        # TODO: add ggwithinstats support for groups with paired samples
        p <- ggscatterstats(plotdf[!is.na(plotdf[,x.coln]),],!!x.coln,!!y,title=y,
                            marginal.type="densigram",type =plot.stattest, ...)
        return(p) }
    })
  }


  plot.list <- setNames(plot.list,y.coln.plot)




  ## generate the ggarrange plot for the output

  # get the num of row and col to decide the plot size
  # following the "arrangeGrob" function used by ggarrange
  # but set a maximum size of 5*4

  plot.nrow.max <- 5
  plot.ncol.max <- 4

  n.plots <- length(plot.list)

  if (is.null(plot.nrow) && !is.null(plot.ncol)) {
    plot.nrow <- ceiling(n.plots/plot.ncol)
    if(plot.nrow>plot.nrow.max) plot.nrow=plot.nrow.max
  }else if (is.null(plot.ncol) && !is.null(plot.nrow)) {
    plot.ncol <- ceiling(n.plots/plot.nrow)
    if(plot.ncol>plot.ncol.max) plot.ncol=plot.ncol.max
  }else if (is.null(plot.nrow) && is.null(plot.ncol)){
    if(n.plots>plot.nrow.max*plot.ncol.max){
      plot.nrow=plot.nrow.max
      plot.ncol=plot.ncol.max
    }else{
      nm <- grDevices::n2mfrow(n.plots)
      plot.nrow = nm[1]
      plot.ncol = nm[2]
    }
  }

  # Arrange the plots into one pdf with ggarrange
  outp <- ggpubr::ggarrange(plotlist = plot.list,
                            nrow = plot.nrow,ncol = plot.ncol,
                            align = "hv")


  # save figure to outpdir
  if(!is.null(outpdir)){
    dir.create(outpdir,showWarnings=F,recursive=T)
    pdf(paste0(outpdir,"/",x.coln,fn.suffix,".pdf"),
        width=plot.w*plot.ncol,height = plot.h*plot.nrow)
    print(outp)
    dev.off()
  }

  # plot to screen
  if(plot.it) {
    print(outp)
  }

  return(list(plot=outp,pval=pvalue.raw))

}



#' @rdname plot_corr
#' @export
plot_corr <- function(plotdf,x.coln,y.coln=NULL,
                      cat.num.test="kruskal.test",
                      cat.cat.test="both",
                      num.num.test="spearman",
                      plot.it=F,plot.nrow=NULL,plot.ncol=NULL,
                      signif.cutoff=0.05,
                      plot.stattest="np",
                      plot.signif.only=F,
                      p.adj.method.each=NULL,
                      p.adj.method.all="bonferroni",
                      plot.max=50,
                      seed=999,
                      outpdir=NULL,
                      plot.w=7,plot.h=7.5,
                      min.group.size.x=3,
                      min.group.size.y=3,
                      fn.suffix="",...){

  # use all columns if y.coln not provided
  if(is.null(y.coln)) y.coln <- setdiff(colnames(plotdf),x.coln)

  feat.cate <- colnames(plotdf)[sapply(plotdf,function(x)!is.numeric(x))]

  lst.out <- lapply(setNames(x.coln,x.coln),
         function(x){
           plot_corr_one(
             plotdf,x,y.coln=setdiff(y.coln,x),
             cat.num.test="kruskal.test",
             num.num.test="spearman",
             plot.it=plot.it,plot.nrow=plot.nrow,plot.ncol=plot.ncol,
             signif.cutoff=signif.cutoff,plot.stattest=plot.stattest,
             p.adj.method=p.adj.method.each,
             plot.signif.only=plot.signif.only,plot.max=plot.max,seed=seed,
             outpdir=outpdir,min.group.size.x=min.group.size.x,
             min.group.size.y=min.group.size.y,fn.suffix=fn.suffix,...
           )
           })

  # save pvalue table to outpdir
  if(!is.null(outpdir)){

    # for(n in names(lst.out)){
    #   if(!is.null(lst.out[[n]]$plot)){
    #
    #     # get the num of row and col to decide the plot size
    #     # following the "arrangeGrob" function used by ggarrange
    #     n.plots <- length(lst.out[[n]]$plot)
    #     if (is.null(plot.nrow) && !is.null(plot.ncol)) {
    #       plot.nrow <- ceiling(n.plots/plot.ncol)
    #     }else if (is.null(plot.ncol) && !is.null(plot.nrow)) {
    #       plot.ncol <- ceiling(n.plots/plot.nrow)
    #     }else if (is.null(plot.nrow) && is.null(plot.ncol)){
    #       nm <- grDevices::n2mfrow(n.plots)
    #       plot.nrow = nm[1]
    #       plot.ncol = nm[2]
    #     }
    #
    #     ggsave(paste0(outpdir,"/",n,".pdf"),lst.out[[n]]$plot,
    #            width=plot.w*plot.ncol,height = plot.h*plot.nrow)
    #   }
    #
    # }

    # save the pvalues as a table
    lst.pval <- lapply(lst.out,"[[","pval") %>% lapply(as.data.frame)
    df.pval <- lapply(names(lst.pval), function(n){
      colnames(lst.pval[[n]]) <- n
      lst.pval[[n]]$Features <- row.names(lst.pval[[n]])
      return(lst.pval[[n]])
    }) %>% Reduce(full_join,.) %>%
      select(Features,everything())

    write.table(df.pval,paste0(outpdir,"/pvalues.tsv"),sep = "\t",row.names = F)

    write.table(cbind(df.pval[,1],p_adjust_mat(df.pval[,-1],p.adj.method.all)),
                paste0(outpdir,"/padj.",p.adj.method.all,".tsv"),sep = "\t",row.names = F)
  }

  return(outp=lst.out)

}
