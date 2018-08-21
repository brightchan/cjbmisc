#' Miscellaneous plotting functions
#'
#' A group of functions for miscellaneous plots.
#' \itemize{
#'   \item perc_barplot: Percentage barplot for two discrete values, showing the proportion of y in each x. fish or chisq test p value shown.
#'   \item violin: violin plot for one discrete and one continuous values. pvalue shown.
#'   \item plot_lm: dot plot with lm regression line for two continuous values. Also show R2, p value of slope
#'   \item violin_mul: plot multiple violin plot, but only the significant ones.
#' }
#' @param ggdf,df dataframe with rows of samples and columns of features
#' @param x,y character of the name of columns to be plotted
#' @param f formula, (y~x)
#' @param title plot title
#' @param col named vector of color (as long as the values are contained in the names), or unnamed vector, then color would be assigned based on the order of the levels.
#' @param xlab,ylab axis label
#' @param width perc_barplot & violin: with of the bar or violin
#' @param aspect.ratio perc_barplot & violin:aspect ratio of the final plot
#' @param show.count perc_barplot & violin:whether to show the total count of each catergory
#' @param show.mean violin:whether to show the mean with a green dot
#' @param anno.hjust,anno.textsize violin:change the position and text size of the significance values
#' @param pch plot_lm:dot style
#' @param plot.it if FALSE, no plot printed, only return the plot object
#' @param ... in perc_barplot and violin: pass to \code{\link[ggplot2]{theme}};\cr in plot_lm: pass to \code{\link[graphics]{plot}}
#' @return
#' @name miscPlots
NULL
#' @rdname miscPlots
#' @export
#' @import ggplot2
#' @importFrom scales percent
#Percentage Barplot
perc_barplot <- function(ggdf,x,y,title=NULL,col=colpal,
                         xlab=x,ylab=y,width=0.9,aspect.ratio=1.4,
                         show.count=TRUE,anno.vjust=0.04,...){
  ggdf[,x] <- droplevels(as.factor(ggdf[,x]))
  #remove NA
  ggdf <- ggdf[complete.cases(ggdf[,y]),]
  # Set names to color according to level of y
  if(is.null(names(col))|!all(levels(as.factor(ggdf[,y]))%in%names(col))){
    #message("Assigning color based on level order.")
    lel <- levels(as.factor(ggdf[,y]))
    col <- setNames(col[1:length(lel)],lel)
  }
  xcount <- table(ggdf[,x]) %>% as.data.frame() %>% mutate(Freq=paste0("n=",Freq),!!y:=NA)
  pvalue <- p_fish.chi.t(ggdf,x,y)
  df <- as.matrix(prop.table(xtabs(as.formula(paste0("~",x,"+",y)),data=ggdf),1))
  df <- data.frame(df)
  p <- ggplot(df,aes_string(x,"Freq",fill=y))+
    geom_bar(stat="identity",colour="black",width=width)+
    scale_y_continuous(labels=scales::percent,limits=c(0,1+anno.vjust),breaks=seq(0,1,0.25))+
    ggtitle(title)+
    labs(x=xlab,y=ylab)+
    annotation_custom(grid::textGrob(paste0("p=",signif(pvalue,2))),
                      xmin = -Inf, xmax = Inf,ymin=1+1.5*anno.vjust,ymax=1+anno.vjust)+
    #guides(fill=F)+
    guides(fill=guide_legend(title=y))+
    scale_fill_manual(values=col)+
    theme_minimal()+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(linetype = "solid",fill=NA),
      #plot.title=element_text(face="bold",hjust=0.5,size=18),
      aspect.ratio=aspect.ratio,...)
  if(show.count) p <- p+geom_text(data=xcount,mapping=aes(x=Var1,y=-anno.vjust,vjust=0.5,label=Freq))+
    scale_y_continuous(labels=scales::percent,limits=c(-anno.vjust,1+anno.vjust),breaks=seq(0,1,0.25))
  return(p)
}
#' @rdname miscPlots
#' @export
#' @import ggplot2
#' @import ggpubr
#' @importFrom ggsignif geom_signif
#' @importFrom EnvStats stat_n_text
#Violin Plot with significance level for 2/3 variables
violin <- function(ggdf,x,y,test="wilcox.test",title=NULL,col=NULL,xlab=x,ylab=y,
                   width=0.9,aspect.ratio=1.4,
                   show.mean=TRUE,show.count=TRUE,
                   anno.vjust=0.1,anno.textsize=4,plot.it=TRUE,...){
  dat <- ggdf[,y]
  spr <- max(dat,na.rm=TRUE)-min(dat,na.rm=TRUE) #spread
  adj <- spr*anno.vjust # unit of adjustment as anno.vjust of the spread
  ylim <- c(min(dat,na.rm=TRUE)-1.5*adj,max(dat,na.rm=TRUE)+2*adj)
  if(is.null(col)){
    p <- ggplot(subset(ggdf,!is.na(get(x))),aes(factor(get(x)),get(y)))+
      geom_violin(width=width)
  } else {
    p <- ggplot(subset(ggdf,!is.na(get(x))),aes(factor(get(x)),get(y)))+
      geom_violin(aes(fill=get(x)),width=width)+
      scale_fill_manual(values=col)+
      guides(fill=guide_legend(title=x))
  }
  p <- p+
    geom_boxplot(width=0.15*width)+
    ggtitle(title)+
    ylab(ylab)+
    xlab(xlab)+
    theme_minimal()+
    theme(aspect.ratio=aspect.ratio,...)
  #comp <- combn(unique(na.omit(ggdf[,x])),2,as.character,simplify=F)
  if (length(na.omit(unique(ggdf[,x])))==2){
    p <- p+ylim(ylim)+
      geom_signif(comparisons=list(c(1,2)),test=test,textsize=anno.textsize,
                  y_position=c(ylim[2]-0.5*adj))
  } else if (length(na.omit(unique(ggdf[,x])))==3){
    p <- p+ylim(ylim)+
      geom_signif(comparisons=list(c(1,2),c(2,3),c(1,3)),test=test,textsize=anno.textsize,
                         y_position=c(ylim[2]-1.4*adj,ylim[2]-1.4*adj,ylim[2]-0.45*adj))
  } else p <- p+stat_compare_means(method=test)
  if (show.count) p <- p+stat_n_text()
  if (show.mean) p <- p+stat_summary(fun.y=mean, geom="point", shape=19, size=2,color="green")
  if(plot.it) print(p)
  return(p)
}
#' @rdname miscPlots
#' @export
#dot plot with lm regression line and R2, p value of slope
plot_lm <- function(df,f,title="",col="#00DD006F",pch=16,...){
  plot(f,df,main=main,xlab=all.names(f)[3],ylab=all.names(f)[2],col=col,pch=pch,...)
  abline(fit <- lm(f,df))
  legend("topright", bty="n", cex=0.7,
         legend=paste("R2 = ", format(summary(fit)$adj.r.squared, digits=2),
                      "\nSlope = ", format(summary(fit)$coefficients[2,1], digits=2),
                      "\nSlope q = ",format(summary(fit)$coefficients[2,4], digits=2)))
}
#' @rdname miscPlots
#' FIXME
#' #export
violin_mul <- function(ggdf,x,y,facet,test="t.test",title=NULL,cutoff=0.05,...){
  library(EnvStats)
  library(ggplot2)
  library(ggsignif)
  if(length(unique(na.omit(ggdf[,x])))<4){
    comp <- combn(unique(na.omit(ggdf[,x])),2,as.character,simplify=F)
    tobeplot <- c() #the list of facets to be plotted
    for(f in unique(ggdf[,facet])){
      for(i in comp){
        pvalue <- tryCatch(
          get(test)(ggdf[ggdf[,facet]==f&ggdf[,x]==i[1],y],
                    ggdf[ggdf[,facet]==f&ggdf[,x]==i[2],y])$p.value,
          error=function(e) return(NA)
        )
        if(is.na(pvalue)) pvalue <- 1 # All pvalue of NA are assigned 1 to avoid error
        if(pvalue<cutoff) tobeplot <- c(tobeplot,f)
      }
    }
    f <- as.formula(paste0("~",facet))
    subset(ggdf,!is.na(get(x))) %>%
      filter(get(facet)%in%tobeplot) %>%
      ggplot(aes(factor(get(x)),get(y)))+
      facet_wrap(f,...)+
      geom_violin()+
      stat_n_text()+
      geom_boxplot(width=0.12)+
      stat_summary(fun.y=mean, geom="point", shape=19, size=2,color="green")+
      ggtitle(title)+
      ylab(y)+
      xlab(x)+
      theme_minimal()+
      theme(text = element_text(size=14))+
      geom_signif(comparisons = comp,test=test,textsize=4)
  }
  else{
    ggplot(subset(ggdf,!is.na(get(x))),aes(factor(get(x)),get(y)))+
      geom_violin()+
      stat_n_text()+
      geom_boxplot(width=0.12)+
      stat_summary(fun.y=mean, geom="point", shape=19, size=2,color="green")+
      ggtitle(title)+
      ylab(y)+
      xlab(x)+
      theme_minimal()+
      theme(text = element_text(size=14))
  }
}
