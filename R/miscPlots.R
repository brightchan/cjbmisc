#' Miscellaneous plotting functions
#'
#' A group of functions for miscellaneous plots. These are the foundation of the plot_feat_subtype and ggpairs_custom\cr
#' You may adjust via trial&error for the anno.hjust and anno.textsize parameter if the p-value text size and position are not good enough.
#' \itemize{
#'   \item perc_barplot: Percentage barplot for two discrete values, showing the proportion of y in each x. fish or chisq test p value shown.
#'   \item violin: violin plot for one discrete and one continuous values. pvalue shown.
#'   \item boxjitter: same as violin plot, but show the data with a box plot and jitter points.
#'   \item plot_lm: dot plot with lm regression line for two continuous values. Also show the correlation, R2, p value of slope
#'   \item violin_mul: plot multiple violin plot, but only the significant ones.
#' }
#' @param ggdf,df dataframe with rows of samples and columns of features
#' @param x,y character of the name of columns to be plotted
#' @param title plot title
#' @param col For perc_barplot & violin: named vector of color (as long as the values are contained in the names), or unnamed vector, then color would be assigned based on the order of the levels.\cr For plot_lm: one color with alpha.
#' @param xlab,ylab axis label
#' @param width perc_barplot & violin: with of the bar or violin
#' @param aspect.ratio perc_barplot & violin:aspect ratio of the final plot
#' @param show.count perc_barplot & violin:whether to show the total count of each catergory
#' @param show.mean violin:whether to show the mean with a green dot
#' @param show.anno plot_lm: whether to show the R2, slope and slope p value
#' @param anno.hjust,anno.textsize violin:change the position and text size of the significance/count values
#' @param border perc_barplot,violin: whether add a solid border to the plot
#' @param box.notch violin/boxjitter: whether show notch style box plot
#' @param highlight.signif,highlight.signif.col,signif.cutoff Whether to highlight (highlight.signif) a significant plot with a color border (highlight.signif.col) according to the p-value cutoff (signif.cutoff).\cr
#'           pvalue calculationg method: perc_barplot-\code{\link{p_fish.chi.t}}; violin-\code{\link{p_krus}}; plot_lm-from the univariate \code{\link[stats]{lm}}
#' @param p.label,p.symnum.args change the label style ("p.signif" or "p.format") of the significance. See label and symnum.args from \code{\link[ggpubr]{stat_compare_means}}
#' @param pch plot_lm:dot style
#' @param plot.it if FALSE, no plot printed, only return the plot object
#' @param cor.method plot_lm:the correlation method used for the annotation, pass to \code{\link[stats]{cor}}
#' @param ... in perc_barplot, violin, plot_lm: pass to \code{\link[ggplot2]{theme}}
#' @return A ggplot object
#' @name miscPlots
NULL
#' @rdname miscPlots
#' @export
#' @import ggplot2
#' @import stats
#' @importFrom scales percent
#Percentage Barplot
perc_barplot <- function(ggdf,x,y,title=NULL,col=comp_hm_colist_full$disc,
                         xlab=x,ylab=y,width=0.9,aspect.ratio=1.4,
                         p.label="p.format",p.symnum.args=list(),
                         show.count=TRUE,anno.textsize=4,anno.vjust=0.04,border=TRUE,plot.it=FALSE,
                         highlight.signif=TRUE,highlight.signif.col="orchid1",signif.cutoff=0.05,...){
  ggdf[,x] <- droplevels(as.factor(ggdf[,x]))
  #remove NA
  ggdf <- ggdf[complete.cases(ggdf[,y]),]
  #get the levels of y
  lel <- levels(as.factor(ggdf[,y]))
  # Set names to color according to level of y
  if(is.null(names(col))|!all(levels(as.factor(ggdf[,y]))%in%names(col))){
    #message("Assigning color based on level order.")
    col <- setNames(col[1:length(lel)],lel)
  }
  xcount <- table(ggdf[,x]) %>% as.data.frame() %>% mutate(Freq=paste0("n=",Freq),!!y:=NA)
  # If y only has one level, return p-value=1
  if(length(lel)==1){
    pvalue <- 1
    message("Only ONE catergory for y axis. Returning p-value=1.")
    } else pvalue <- p_fish.chi.t(ggdf,x,y)

  # Change p value format
  if(p.label=="p.signif") {
    if(length(p.symnum.args)==0) {
      pvalue <- symnum(pvalue,cutpoints=c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                       symbols=c("****", "***", "**", "*", "ns"))
    } else pvalue <- symnum(pvalue,cutpoints=p.symnum.args[[1]],
                            symbols=p.symnum.args[[2]])
  } else {pvalue <- signif(pvalue,2)}

  df <- as.matrix(prop.table(xtabs(as.formula(paste0("~",x,"+",y)),data=ggdf),1))
  df <- data.frame(df)
  p <- ggplot(df,aes_string(x,"Freq",fill=y))+
    geom_bar(stat="identity",colour="black",width=width)+
    scale_y_continuous(labels=scales::percent,limits=c(0,1+anno.vjust),breaks=seq(0,1,0.25))+
    ggtitle(title)+
    labs(x=xlab,y=ylab)+
    annotation_custom(grid::textGrob(paste0("p=",pvalue),gp = grid::gpar(fontsize = 3*anno.textsize)),
                      xmin = -Inf, xmax = Inf,ymin=1+1.5*anno.vjust,ymax=1+anno.vjust)+
    #guides(fill=F)+
    guides(fill=guide_legend(title=y))+
    scale_fill_manual(values=col)+
    theme_minimal()+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      #plot.title=element_text(face="bold",hjust=0.5,size=18),
      aspect.ratio=aspect.ratio,...)

  # Add a frame
  if (border) p <- p+theme(panel.border = element_rect(linetype = "solid",fill=NA))

  # Highlight significant plots with a color frame
  if (highlight.signif&pvalue<signif.cutoff) p <- p +
    theme(panel.border=element_rect(linetype = "solid",fill=NA,colour=highlight.signif.col))

  if(show.count) suppressMessages(p <- p+geom_text(data=xcount,mapping=aes(x=Var1,y=-anno.vjust,vjust=0.5,label=Freq),size=anno.textsize)+
    scale_y_continuous(labels=scales::percent,limits=c(-anno.vjust,1+anno.vjust),breaks=seq(0,1,0.25)))
  if(plot.it) print(p)
  return(p)
}
#' @rdname miscPlots
#' @export
#' @import ggplot2
#' @import ggpubr
#' @importFrom ggpubr stat_compare_means
#' @importFrom EnvStats stat_n_text
#' @importFrom gginnards which_layers
#Violin Plot with significance level for 2/3 variables
violin <- function(ggdf,x,y,test="wilcox.test",
                   title=NULL,col=NULL,xlab=x,ylab=y,
                   width=0.9,aspect.ratio=1.4,
                   show.mean=TRUE,show.count=TRUE,
                   p.label=NULL,p.symnum.args=list(),
                   anno.vjust=0.1,anno.textsize=4,
                   border=TRUE,plot.it=FALSE,
                   highlight.signif=TRUE,highlight.signif.col="orchid1",signif.cutoff=0.05,
                   plot.type="violin",box.notch=FALSE,...){
  dat <- ggdf[,y]
  spr <- max(dat,na.rm=TRUE)-min(dat,na.rm=TRUE) #spread
  adj <- spr*anno.vjust # unit of adjustment as anno.vjust of the spread
  # Adjust annotation position
  if(show.count){
    ylim <- c(min(dat,na.rm=TRUE)-1.5*adj,max(dat,na.rm=TRUE)+2*adj)
  } else ylim <- c(min(dat,na.rm=TRUE),max(dat,na.rm=TRUE)+2*adj)

  if(plot.type=="violin"){
    #Assign color
    if(is.null(col)){
      p <- ggplot(subset(ggdf,!is.na(get(x))),aes(factor(get(x)),get(y)))+
        geom_violin(width=width)+guides(fill="none")
    } else {
      p <- ggplot(subset(ggdf,!is.na(get(x))),aes(factor(get(x)),get(y)))+
        geom_violin(aes(fill=get(x)),width=width)+
        scale_fill_manual(values=col)+
        guides(fill=guide_legend(title=x))
    }
    p <- p+geom_boxplot(width=0.13*width)
    #comp <- combn(unique(na.omit(ggdf[,x])),2,as.character,simplify=F)
  } else if(plot.type=="boxplot"){
    width <- width*0.85
    if(is.null(col)){
      p <- ggplot(subset(ggdf,!is.na(get(x))),aes(factor(get(x)),get(y)))+
        geom_boxplot(width=width,notch = box.notch,outlier.shape = NA)+
        guides(fill="none")
    } else {
      p <- ggplot(subset(ggdf,!is.na(get(x))),aes(factor(get(x)),get(y)))+
        geom_boxplot(aes(fill=get(x)),width=width,notch = box.notch,outlier.shape = NA)+
        scale_fill_manual(values=col)+
        guides(fill=guide_legend(title=x))
    }
    p <- p+geom_jitter(width=0.35*width,alpha=0.6,shape=19)
  } else {stop("Wrong plot.type param: can only be 'violin' or 'boxplot'")}

  # Add title and labels
  p <- p+
    ggtitle(title)+
    ylab(ylab)+
    xlab(xlab)+
    theme_minimal()+
    theme(aspect.ratio=aspect.ratio,...)

  # Add significant annotations
  if (length(na.omit(unique(ggdf[,x])))==2){
    p <- p+ylim(ylim)+
      stat_compare_means(comparisons=list(c(1,2)),method=test,
                         p.symnum.args=p.symnum.args,label=p.label,
                         label.y=c(ylim[2]-0.5*adj))
    #geom_signif(comparisons=list(c(1,2)),test=test,textsize=anno.textsize,
    #            y_position=c(ylim[2]-0.5*adj))
    p$layers[[which_layers(p, "GeomSignif")]]$aes_params$textsize <- anno.textsize
  } else if (length(na.omit(unique(ggdf[,x])))==3){
    p <- p+ylim(ylim)+
      stat_compare_means(comparisons=list(c(1,2),c(2,3),c(1,3)),method=test,
                         p.symnum.args=p.symnum.args,label=p.label,
                         label.y=c(ylim[2]-1.4*adj,ylim[2]-1.4*adj,ylim[2]-0.45*adj))
    p$layers[[which_layers(p, "GeomSignif")]]$aes_params$textsize <- anno.textsize
  } else p <- p+stat_compare_means(method="kruskal.test",
                                   label.y=c(ylim[2]-1.3*adj),
                                   size=anno.textsize)+
    ylim(c(ylim[1],ylim[2]-adj))

  # Add a frame
  if (border) p <- p+theme(panel.border = element_rect(linetype = "solid",fill=NA))

  # Highlight significant plots with a color frame
  ggdf[,x] <- as.factor(ggdf[,x])

  if (highlight.signif){
    if(length(na.omit(unique(ggdf[,x])))==2) pvalue <- p_ContDisc(ggdf,x,y,method=test) else pvalue <- p_ContDisc(ggdf,x,y,method="kruskal.test")
    if(pvalue<signif.cutoff) p <- p +
    theme(panel.border=element_rect(linetype = "solid",fill=NA,colour=highlight.signif.col))
    }

  if (show.count) p <- p+stat_n_text(size=anno.textsize)
  if (show.mean) p <- p+stat_summary(fun.y=mean, geom="point", shape=19, size=1.5,color="#22DD2290")
  if(plot.it) print(p)
  return(p)
}
#' @rdname miscPlots
#' @export
boxjitter <- function(...){
  violin(...,plot.type="boxplot")
}


#' @rdname miscPlots
#' @export
#dot plot with lm regression line and R2, p value of slope
# plot_lm <- function(df,f,title=NULL,col="#00DD006F",pch=16,...){
#   plot(f,df,main=title,xlab=all.names(f)[3],ylab=all.names(f)[2],col=col,pch=pch,...)
#   abline(fit <- lm(f,df))
#   legend("topright", bty="n", cex=0.7,
#          legend=paste("R2 = ", format(summary(fit)$adj.r.squared, digits=2),
#                       "\nSlope = ", format(summary(fit)$coefficients[2,1], digits=2),
#                       "\nSlope q = ",format(summary(fit)$coefficients[2,4], digits=2)))
# }

plot_lm <- function(ggdf,x,y,title=NULL,col="#00DD006F",show.anno=TRUE,cor.method="spearman",
                    plot.it=TRUE,highlight.signif=TRUE,
                    highlight.signif.col="orchid1",signif.cutoff=0.05,...){
  fit <- lm(as.formula(paste0(y,"~",x)),ggdf)
  if(show.anno){
    anno <- paste0("R2=", format(summary(fit)$adj.r.squared, digits=2),
                  " Slope=", format(summary(fit)$coefficients[2,1], digits=2),
                  " Slope_p=",format(summary(fit)$coefficients[2,4], digits=2))
    anno.cor <- paste0(cor.method,".cor=\n",signif(cor(ggdf[,x],ggdf[,y],method=cor.method,use="com"),2))
  } else {anno <- NULL;anno.cor <- ""}
  p <- ggplot(ggdf,aes_string(x,y))+
    geom_point(alpha = 1/2,color="green")+
    scale_fill_manual(values=col)+
    labs(title=title,caption=anno)+
    geom_smooth(method=lm,fill="grey80")+
    annotate("text",x=-Inf,y=Inf,label=anno.cor,hjust=-0.1,vjust=1.1,col="#00000093")+
    theme_minimal()+
    theme(...)

  # Highlight significant plots with a color frame
  pvalue <- summary(fit)$coefficients[2,4]
  if (highlight.signif&pvalue<signif.cutoff) p <- p +
    theme(panel.border=element_rect(linetype = "solid",fill=NA,colour=highlight.signif.col))

  if(plot.it) print(p)
  return(p)
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
