#' Plot pairwise correlations of multiple features
#'
#' Plot a big complex correlation map of every features in a dataframe or specific features (feat.plot).
#' Usage: ggpairs_custom(df,plot.it=TRUE) (all features). Or p<-ggpairs_custom(df,c("Ploidy","Purity","EGFR","Smoker"))\cr
#' then ggsave("plot.pdf",p$plot)
#'
#' @param ggdf dataframe with rows of samples and columns of features
#' @param feat.plot character vector of the features to be plotted. If NULL then all columns will be plotted
#' @param colist a named list of named vector of colors to be used for each feature. Missing ones are automatically handled.
#' @param sig.col vector of 3 colors for nonsignif, lower and higher bound of the significance.
#' @param signif.cutoff cutoff point for significant p-value. A colored (sig.col[3]) frame will be shown to significant plots.
#' @param plot.it if FALSE, no plot printed, only return the object. Use ggsave(plot=p$plot) to plot.
#' @param ... pass to \code{\link[GGally]{ggpairs}}
#' @return a list of two: $plot is "gg" "ggmatrix" plotable and ggsavable object; $p.value is the table of correlation pvalues
#' @name ggpairs_custom
#' @import ggplot2
#' @import RColorBrewer
#' @import GGally
#' @export
ggpairs_custom <- function(ggdf,feat.plot=NULL,colist=comp_hm_colist_full,
                           sig.col=c("white","thistle1","orchid1"),plot.it=FALSE,
                           signif.cutoff=0.05,...){
  ####Define sub-plots.
  ####This is also the storage of functions can be used by ggpairs
  { ### For P-values
    # Continuous
    ggplot_corr<- function(data, mapping, method="pearson",col=sig.col){
      x <- data[,quo_name(mapping$x)]
      y <- data[,quo_name(mapping$y)]
      cor <- cor(x, y, method = method,use="complete.obs")
      fill=colorRampPalette(c(col))(10)[scales::rescale(cor,c(1,10),c(0,1))]
      xrange = c(0, 1)
      yrange = c(0, 1)
      ggplot() + xlim(xrange) + ylim(yrange)+
        geom_rect(aes(xmin = xrange[1], xmax = xrange[2], ymin = yrange[1], ymax = yrange[2]),
                  fill = fill)+
        annotate("text", x = mean(xrange), y = mean(yrange), label = paste0("Corr=",round(cor, digits = 4)))+
        theme_minimal()
    }

    ggplot_lmPvalue <- function(data, mapping, col=sig.col){
      x <- data[,quo_name(mapping$x)]
      y <- data[,quo_name(mapping$y)]
      pvalue <- coef(summary(lm(y~x)))[2,4]
      # add pvalue to corrTable
      corrTable[quo_name(mapping$x),quo_name(mapping$y)] <<- pvalue
      corrTable[quo_name(mapping$y),quo_name(mapping$x)] <<- pvalue
      fill=ifelse(pvalue>0.01,col[1],
                  colorRampPalette(col[2:3])(100)[scales::rescale(-log(pvalue),c(1,100),c(0,200))])
      #fill=ifelse(pvalue<0.01,col[3],ifelse(pvalue<0.001,col[2],col[1]))
      xrange = c(0, 1)
      yrange = c(0, 1)
      ggplot() + xlim(xrange) + ylim(yrange)+
        geom_rect(aes(xmin = xrange[1], xmax = xrange[2], ymin = yrange[1], ymax = yrange[2]),
                  fill = fill)+
        annotate("text", x = mean(xrange), y = mean(yrange), label = paste0(round(pvalue, digits = 4)))+
        theme(axis.line=element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              legend.position="none",
              panel.background=element_blank(),
              panel.border=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              plot.background=element_blank())
    }
    # Contin + Discret
    ggplot_anova<- function(data, mapping,col=sig.col){
      if (plyr::is.discrete(data[,quo_name(mapping$y)])){
        f <- as.formula(paste0(quo_name(mapping$y),"~",quo_name(mapping$x)))
      } else f <- as.formula(paste0(quo_name(mapping$x),"~",quo_name(mapping$y)))
        aov <- aov(f, data)
        pvalue <- summary(aov)[[1]][["Pr(>F)"]][1]
        # add pvalue to corrTable
        corrTable[quo_name(mapping$x),quo_name(mapping$y)] <<- pvalue
        corrTable[quo_name(mapping$y),quo_name(mapping$x)] <<- pvalue
        #fill=colorRampPalette(c(col))(10)[scales::rescale(aov.p,c(10,1),c(0,1))]
        fill=ifelse(pvalue<0.01,col[3],ifelse(pvalue<0.05,col[2],col[1]))
        xrange = c(0, 1)
        yrange = c(0, 1)
        ggplot() + xlim(xrange) + ylim(yrange)+
          geom_rect(aes(xmin = xrange[1], xmax = xrange[2], ymin = yrange[1], ymax = yrange[2]),
                    fill = fill)+
          annotate("text", x = mean(xrange), y = mean(yrange), label = paste0(round(pvalue, digits = 2)))+
          theme(axis.line=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks=element_blank(),
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                legend.position="none",
                panel.background=element_blank(),
                panel.border=element_blank(),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                plot.background=element_blank())
    }

    ggplot_krus <- function(data, mapping,col=sig.col){
      if (plyr::is.discrete(data[,quo_name(mapping$y)])){
        f <- as.formula(paste0(quo_name(mapping$x),"~",quo_name(mapping$y)))
      } else f <- as.formula(paste0(quo_name(mapping$y),"~",quo_name(mapping$x)))
        pvalue <- kruskal.test(f,data)$p.value
        # add pvalue to corrTable
        corrTable[quo_name(mapping$x),quo_name(mapping$y)] <<- pvalue
        corrTable[quo_name(mapping$y),quo_name(mapping$x)] <<- pvalue
        #fill=colorRampPalette(c(col))(10)[scales::rescale(aov.p,c(10,1),c(0,1))]
        fill=ifelse(pvalue<0.01,col[3],ifelse(pvalue<0.05,col[2],col[1]))
        xrange = c(0, 1)
        yrange = c(0, 1)
        ggplot() + xlim(xrange) + ylim(yrange)+
          geom_rect(aes(xmin = xrange[1], xmax = xrange[2], ymin = yrange[1], ymax = yrange[2]),
                    fill = fill)+
          annotate("text", x = mean(xrange), y = mean(yrange), label = paste0(round(pvalue, digits = 2)))+
          theme(axis.line=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks=element_blank(),
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                legend.position="none",
                panel.background=element_blank(),
                panel.border=element_blank(),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                plot.background=element_blank())
    }

    # Discret
    ggplot_FET<- function(data, mapping,col=sig.col){
      options(workspace=2e9)
      pvalue <- p_fish.chi.t(data,quo_name(mapping$x),quo_name(mapping$y))
      # add pvalue to corrTable
      corrTable[quo_name(mapping$x),quo_name(mapping$y)] <<- pvalue
      corrTable[quo_name(mapping$y),quo_name(mapping$x)] <<- pvalue
      #fill=colorRampPalette(c(col))(10)[scales::rescale(pvalue,c(10,1),c(0,1))]
      fill=ifelse(pvalue<0.01,col[3],ifelse(pvalue<0.05,col[2],col[1]))
      xrange = c(0, 1)
      yrange = c(0, 1)
      ggplot() + xlim(xrange) + ylim(yrange)+
        geom_rect(aes(xmin = xrange[1], xmax = xrange[2], ymin = yrange[1], ymax = yrange[2]),
                  fill = fill)+
        annotate("text", x = mean(xrange), y = mean(yrange), label = paste0(round(pvalue, digits = 2)))+
        theme(axis.line=element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              legend.position="none",
              panel.background=element_blank(),
              panel.border=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              plot.background=element_blank())
    }

    ### For actual plots
    ggplot_percBar <- function(data, mapping, col=colist){
      x <- quo_name(mapping$x)
      y <- quo_name(mapping$y)
      # Assign color
      if(y%in%names(col)) {
        col <- col[[y]]
        } else col <- comp_hm_colist_full$disc
      # df <- as.matrix(prop.table(xtabs(as.formula(paste0("~",x,"+",y)),data=data),1))
      # df <- data.frame(df)
      # #pvalue <- fish.t(x,y,data)$p.value
      # ggplot(df,aes_string(x,"Freq",fill=y))+
      #   geom_bar(stat="identity",colour="black")+
      #   scale_fill_manual(values=col)+
      #   theme_minimal()
      perc_barplot(data,x=x,y=y,aspect.ratio=1,col=col,anno.textsize=3,signif.cutoff=signif.cutoff,border=F)
    }

    ggplot_box <- function(data, mapping,col=colist){
      if (plyr::is.discrete(data[,quo_name(mapping$x)])){
        x <- quo_name(mapping$x)
        y <- quo_name(mapping$y)
      } else {
        x <- quo_name(mapping$y)
        y <- quo_name(mapping$x)
      }
      if(x%in%names(col)) {
        col <- col[[x]]
      } else col <- comp_hm_colist_full$disc
      boxjitter(data,x,y,aspect.ratio=1,col=col,border=F,plot.it=FALSE,
                anno.textsize=3,signif.cutoff=signif.cutoff)
    }

    ggplot_violin <- function(data, mapping,col=colist){
      if (plyr::is.discrete(data[,quo_name(mapping$x)])){
        x <- quo_name(mapping$x)
        y <- quo_name(mapping$y)
      } else {
        x <- quo_name(mapping$y)
        y <- quo_name(mapping$x)
      }
      if(x%in%names(col)) {
        col <- col[[x]]
      } else col <- comp_hm_colist_full$disc
      violin(data,x,y,aspect.ratio=1,col=col,border=F,plot.it=FALSE,
             anno.textsize=3,signif.cutoff=signif.cutoff)
    }


    ggplot_scatter <- function(data, mapping,col=comp_hm_colist_full$disc){
      ggplot(data,mapping)+
        geom_point(alpha = 1/2,color="green")+
        scale_fill_manual(values=col)+
        geom_smooth(method=lm)+
        theme_minimal()
    }

    ggplot_lm <- function(data, mapping,col=comp_hm_colist_full$disc){
      x <- quo_name(mapping$x)
      y <- quo_name(mapping$y)
      p <- plot_lm(data,x,y,plot.it=FALSE,signif.cutoff=signif.cutoff)
      lab <- gsub(" ","\n",p$labels$caption)
      p+annotate("text",label=lab, x=Inf, y = Inf,vjust=1.2,hjust=1.2,size=2.5)
    }

    ggplot_headerblank <- function(data, mapping,col=NULL){
      ggplot() + xlim(xrange) + ylim(yrange)+
        geom_rect(aes(xmin = xrange[1], xmax = xrange[2], ymin = yrange[1], ymax = yrange[2]),
                  fill = fill)+
        annotate("text", x = mean(xrange), y = mean(yrange), label = mapping$x)+
        theme(axis.line=element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              legend.position="none",
              panel.background=element_blank(),
              panel.border=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              plot.background=element_blank())
    }
  }

  if(is.null(feat.plot)) feat.plot <- colnames(ggdf)

  corrTable <- matrix(1,length(feat.plot),length(feat.plot),dimnames=list(feat.plot,feat.plot))

  p <- ggpairs(ggdf[,feat.plot],axisLabels="internal",
          upper=list(continuous =ggplot_lmPvalue,
                     combo = ggplot_krus,
                     discrete = ggplot_FET),
          lower=list(continuous = ggplot_lm,
                     combo=ggplot_violin,
                     discrete = ggplot_percBar),...)

  if(plot.it) print(p)
  return(list(plot=p,p.value=corrTable))
}
