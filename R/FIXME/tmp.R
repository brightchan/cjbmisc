plot_feat_subtype <- function(plotdf,subtype,feat.plot,cont.plot="boxplot",
                              cont.test="wilcox.test",disc.test="both",
                              plot.it=T,plot.nrow=NULL,plot.ncol=NULL,
                              p.label=NULL,p.symnum.args=list(),
                              cont.show.mean=FALSE,show.count=TRUE,
                              cont.col=comp_hm_colist_full$disc,disc.col=comp_hm_colist_full,
                              disc.width=0.8,cont.width=0.9,
                              disc.border=TRUE,cont.border=TRUE,
                              anno.textsize=4,cont.anno.vjust=0.2,disc.anno.vjust=0.04,
                              aspect.ratio=1.4,legend.position="bottom",
                              highlight.signif=TRUE,highlight.signif.col="red",signif.cutoff=0.05,
                              plot.signif.only=F,...){

  feat.disc <- colnames(plotdf)[sapply(plotdf,function(x)!is.numeric(x))]
  feat.cont <- colnames(plotdf)[sapply(plotdf,function(x)is.numeric(x))]


  message("Plotting: ",appendLF=F)
  plot.list <- lapply(feat.plot,function(x){
    message(x,appendLF=F)
    # Discrete
    if(x %in% feat.disc){
      message(":disc ",appendLF=F)
      if (x %in% names(disc.col)) {
        col=disc.col[[x]]
      } else{
        if (length(unique(na.omit(plotdf[,x])))==2) col=comp_hm_colist_full$bool else col=comp_hm_colist_full$disc
      }
      p <- perc_barplot(plotdf,subtype,x,title=x,show.count=show.count,width=disc.width,
                        col=col,xlab="",ylab="",anno.textsize=anno.textsize,
                        anno.vjust=disc.anno.vjust,aspect.ratio=aspect.ratio,border=disc.border,
                        highlight.signif=highlight.signif,highlight.signif.col=highlight.signif.col,signif.cutoff=signif.cutoff,...)
    } else {
      # Continuous
      message(":cont ",appendLF=F)
      p <- violin(plotdf[!is.na(plotdf[,subtype]),],subtype,x,
                  title=x,col=cont.col,width=cont.width,anno.textsize=anno.textsize,
                  show.mean=cont.show.mean,show.count=show.count,
                  plot.type=cont.plot,border=cont.border,aspect.ratio=aspect.ratio,
                  p.label=p.label,p.symnum.args=p.symnum.args,
                  highlight.signif=highlight.signif,highlight.signif.col=highlight.signif.col,signif.cutoff=signif.cutoff,...)+
        theme(axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.title=element_blank())
    }
    p <- p+theme(legend.position=legend.position,
                 ...)
    return(p)
  })
  plot.list <- setNames(plot.list,feat.plot)

  # Arrange by pvalues
  pvlaue <- p_feat_subtype(plotdf,subtype,feat.plot,cont.test=cont.test,disc.test=disc.test)
  plot.list <- plot.list[order(pvalue)]

  if(plot.signif.only){
    plot.list <- plot.list[(pvlaue<signif.cutoff)[order(pvalue)]]
  }else{
    plot.list <- plot.list
  }

  if(plot.it) {
    if(is.null(plot.nrow)|is.null(plot.ncol)) ggpubr::ggarrange(plotlist = plot.list)
    else  ggpubr::ggarrange(plotlist = plot.list,nrow = plot.nrow,ncol = plot.ncol)}

  return(plot.list)
}
