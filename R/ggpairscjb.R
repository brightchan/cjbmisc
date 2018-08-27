#' Plot pairwise correlations of multiple features
#'
#' Only a beta version
#'
#' @param ggdf dataframe with rows of samples and columns of features
#' @param feat.plot character vector of the features to be plotted. If NULL then all columns will be plotted
#' @param ... pass to \code{\link[GGally]{ggpairs}}
#' @return a "gg" "ggmatrix" plotable and ggsavable object
#' @name ggpairs_custom
#' @import ggplot2
#' @import RColorBrewer
#' @import GGally
#' @export
ggpairs_custom <- function(ggdf,feat.plot=NULL,...){
  ####Define sub-plots
  { ### For P-values
    # Continuous
    ggplot_corr<- function(data, mapping, method="pearson",col=c("white","thistle1","red")){
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

    ggplot_lmPvalue <- function(data, mapping, col=c("white","thistle1","orchid1")){
      x <- data[,quo_name(mapping$x)]
      y <- data[,quo_name(mapping$y)]
      pvalue <- coef(summary(lm(y~x)))[2,4]
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
    ggplot_anova<- function(data, mapping,col=c("white","thistle1","orchid1")){
      if (plyr::is.discrete(eval(mapping$x,data))){
        f <- as.formula(paste0(quo_name(mapping$y),"~",quo_name(mapping$x)))
      } else f <- as.formula(paste0(quo_name(mapping$x),"~",quo_name(mapping$y)))
        aov <- aov(f, data)
        pvalue <- summary(aov)[[1]][["Pr(>F)"]][1]
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

    ggplot_krus <- function(data, mapping,col=c("white","thistle1","orchid1")){
      if (plyr::is.discrete(eval(mapping$y,data))){
        f <- as.formula(paste0(quo_name(mapping$x),"~",quo_name(mapping$y)))
      } else f <- as.formula(paste0(quo_name(mapping$y),"~",quo_name(mapping$x)))
        pvalue <- kruskal.test(f,data)$p.value
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
    ggplot_FET<- function(data, mapping,col=c("white","thistle1","orchid1")){
      options(workspace=2e9)
      pvalue <- p_fish.chi.t(data,quo_name(mapping$x),quo_name(mapping$y))
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
    ggplot_percBar <- function(data, mapping,col=comp_hm_colist_full$disc){
      x <- quo_name(mapping$x)
      y <- quo_name(mapping$y)
      df <- as.matrix(prop.table(xtabs(as.formula(paste0("~",x,"+",y)),data=data),1))
      df <- data.frame(df)
      #pvalue <- fish.t(x,y,data)$p.value
      ggplot(df,aes_string(x,"Freq",fill=y))+
        geom_bar(stat="identity",colour="black")+
        scale_fill_manual(values=col)+
        theme_minimal()
    }

    ggplot_box <- function(data, mapping,col=comp_hm_colist_full$disc){
      if (plyr::is.discrete(data[,quo_name(mapping$x)])){
        x <- quo_name(mapping$x)
        y <- quo_name(mapping$y)
      } else {
        x <- quo_name(mapping$y)
        y <- quo_name(mapping$x)
      }
      if(is.null(col)) col <- comp_hm_colist_full$disc
      ggplot(data[complete.cases(data[,as.character(c(x,y))]),],
             aes_string(x,y,fill=x))+
        geom_boxplot()+
        scale_fill_manual(values=col)+
        theme_minimal()
    }

    ggplot_scatter <- function(data, mapping,col=comp_hm_colist_full$disc){
      ggplot(data,mapping)+
        geom_point(alpha = 1/2,color="green")+
        scale_fill_manual(values=col)+
        geom_smooth(method=lm)+
        theme_minimal()
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

  ggpairs(ggdf[,feat.plot],axisLabels="internal",
          upper=list(continuous =ggplot_lmPvalue,
                     combo = ggplot_krus,
                     discrete = ggplot_FET),
          lower=list(continuous = ggplot_scatter,
                     combo=ggplot_box,
                     discrete = ggplot_percBar),...)
}
