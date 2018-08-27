#' Plot multiple features of subtypes
#'
#' Plot percentage barplot or violin/boxplot for discrete or continuous variables of each subtype.
#' Return a list of ggplots that can be arranged into figure panel with ggarrange.
#'
#' @param plotdf dataframe with rows of samples and columns of features, one column should contain the subtype info.
#' @param subtype column name of the subtype
#' @param feat.plot character vector of the features to be plotted
#' @param cont.plot character of either "boxplot" or "violin"
#' @param cont.test the test to be used for the significant test for continuous variables
#' @param cont.show.mean boolean of whether to show the mean as a green dot for continuous variables
#' @param show.count whether to show the total count of each catergory
#' @param disc.col list of named vectors for colors of discrete variables
#' @param cont.col vector of colors for continuous variables
#' @param cont.anno.vjust,disc.anno.vjust numeric to tweak the position of the pvalue annotation on the plot
#' @param label,symnum.args change the label style ("p.signif" or "p.format") of the significance. See \code{\link[ggpubr]{stat_compare_means}}
#' @param anno.textsize change text size of the significance values
#' @name plot_feat_subtype
#' @import ggplot2
#' @importFrom ggpubr stat_compare_means
#' @importFrom EnvStats stat_n_text
#' @importFrom gginnards which_layers
#' @export
plot_feat_subtype <- function(plotdf,subtype,feat.plot,cont.plot="boxplot",
                              cont.test="wilcox.test",label=NULL,symnum.args=list(),
                              cont.show.mean=FALSE,show.count=TRUE,
                              cont.col=comp_hm_colist_full$disc,disc.col=comp_hm_colist_full,
                              disc.width=0.8,cont.width=0.9,
                              anno.textsize=4,cont.anno.vjust=0.2,disc.anno.vjust=0.04,
                              aspect.ratio=1.4){
  feat.disc <- colnames(plotdf)[sapply(plotdf,function(x)!is.numeric(x))]
  feat.cont <- colnames(plotdf)[sapply(plotdf,function(x)is.numeric(x))]

  para.plot <- lapply(feat.cont,function(x){ #parameter for continuous plots
    dat <- plotdf[,x]
    spr <- max(dat,na.rm=TRUE)-min(dat,na.rm=TRUE) #spread
    adj <- spr*cont.anno.vjust # unit of adjustment as anno.vjust of the spread
    if(show.count){
      ylim <- c(min(dat,na.rm=TRUE)-1.5*adj,max(dat,na.rm=TRUE)+2*adj)
    } else ylim <- c(min(dat,na.rm=TRUE),max(dat,na.rm=TRUE)+2*adj)
    return(ylim) #the bottom, the top of the plot, and the adj value
  }) %>% setNames(feat.cont)

  plot.list <- lapply(feat.plot,function(x){
    if(x %in% feat.disc){
      if (x %in% names(disc.col)) {
        col=comp_hm_colist_full[[x]]
      } else{
        if (length(unique(na.omit(plotdf[,x])))==2) col=comp_hm_colist_full$bool else col=comp_hm_colist_full$disc
      }
      p <- perc_barplot(plotdf,subtype,x,title=x,show.count=show.count,width=disc.width,
                        col=col,xlab="",ylab="",anno.textsize=anno.textsize,
                        anno.vjust=disc.anno.vjust,aspect.ratio=aspect.ratio)
    } else {
      ylim=para.plot[[x]][1:2]
      adj=para.plot[[x]][3]
      if (cont.plot=="boxplot"){
        p <- ggboxplot(plotdf,subtype,x,fill=subtype,palette=cont.col,width=cont.width,
                       add = "jitter",add.params=list(size=0.2,color="grey30",alpha=0.6),
                       outlier.shape=NA)
      } else if (cont.plot=="violin"){
        p <- ggviolin(plotdf,subtype,x,fill=subtype,palette=cont.col,width=cont.width,
                      add = "boxplot",add.params = list(fill = "grey90",width=0.15*cont.width))
      }
      # Add theme features
      p <- p+
        ggtitle(x)+ ylim(ylim) +
        theme_minimal()+theme(axis.title.x=element_blank(),
                              axis.title.y=element_blank(),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              #axis.text = element_text(size=15,color="black"),
                              panel.border = element_rect(linetype = "solid",fill=NA),
                              #plot.title=element_text(face="bold",hjust=0.5,size=18),
                              aspect.ratio=aspect.ratio)
      # Add count
      if(show.count) p <- p+stat_n_text(y.pos=ylim[1]+0.15*adj)
      # add significance
      if (length(na.omit(unique(plotdf[,subtype])))==2){
        p <- p+
          stat_compare_means(comparisons=list(c(1,2)),method=cont.test,
                             symnum.args=symnum.args,label=label,
                             label.y=c(ylim[2]-0.5*adj))
        p$layers[[which_layers(p, "GeomSignif")]]$aes_params$textsize <- anno.textsize
      } else if (length(na.omit(unique(plotdf[,subtype])))==3){
        p <- p+
          stat_compare_means(comparisons=list(c(1,2),c(2,3),c(1,3)),method=cont.test,
                             symnum.args=symnum.args,label=label,
                             label.y=c(ylim[2]-1.4*adj,ylim[2]-1.4*adj,ylim[2]-0.45*adj))
        p$layers[[which_layers(p, "GeomSignif")]]$aes_params$textsize <- anno.textsize
      } else p <- p+stat_compare_means(method="kruskal.test",
                                       label.y=c(ylim[2]-1.3*adj),
                                       size=anno.textsize)+
        ylim(c(ylim[1],ylim[2]-adj))
      if (cont.show.mean) p <- p+stat_summary(fun.y=mean, geom="point", shape=19, size=2,color="green")
    }
    return(p)
  })
}
