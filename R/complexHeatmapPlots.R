#' Complex heatmap plots
#'
#' A group of functions for complexheatmap plots using simplied single functions.
#' \itemize{
#'   \item subtypeChaHM_col: Plot complexheatmap of certain characteristics (Age,Stage,etc), sorting by subtype, using specified color.
#' }
#' @param df dataframe
#' @param subtype character of the name of column in the df, containing subtype info
#' @param cha vector of characters to be plotted
#' @param topbar character of the name of a numeric column with which a bar will be plotted on the top
#' @param order_by character of the name of column to order the heatmap
#' @param pic function of plotting device, eg png, pdf, etc.
#' @param outp character of output folder and file name
#' @param w,h width and height of the output plot
#' @param colist the comprehensive color list.\cr By default it is comp_hm_colist_full loaded in this package
#' @return a plot saved in desinated path
#' @seealso \code{\link[ComplexHeatmap]{Heatmap}}
#' @name compHeatmap
NULL
#' @rdname compHeatmap
#' @export
#' @import ComplexHeatmap
#' @importFrom circlize colorRamp2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grid gpar
#
#input dataframe, with the subtype, and characters to be shown as heatmap for each subtype.
#plot using pic function, outp file and size w,h
# if topbar is a column, sort by that column and add top annotation bar
subtypeChaHM_col <- function(df,subtype,cha,topbar,order_by=topbar,pic,outp,w,h,colist=comp_hm_colist_full){
  lel <- levels(as.factor(df[,subtype]))
  colsubtype <- setNames(RColorBrewer::brewer.pal(8,"YlGnBu")[c(2,4,6,8)],lel)
  df <- df[order(df[,order_by]),]
  colR <- circlize::colorRamp2(breaks = c(min(df[,topbar],na.rm=T),
                                median(df[,topbar],na.rm=T),
                                max(df[,topbar],na.rm=T)),
                     colors = RColorBrewer::brewer.pal(8,"YlGnBu")[c(2,5,8)])
  # Assign color
  for (i in cha){
    if (is.numeric(df[,i])){
      #for continuous, take max,med and min for color assignment
      ma <- max(df[,i],na.rm=T)
      mi <- min(df[,i],na.rm=T)
      me <- median(df[,i],na.rm=T)
      if (i %in% names(colist)) colist[[i]] <- circlize::colorRamp2(c(mi,me,ma),colist[[i]])
      else colist[[i]] <- circlize::colorRamp2(c(mi,me,ma),colist[["cont"]])
      next
    }
    if (!i %in% names(colist)){
      lel <- levels(as.factor(df[,i]))
      if (length(lel)==2) colist[[i]] <- setNames(colist[["bool"]],lel)
      else colist[[i]] <- setNames(colist[["disc"]][1:length(lel)],lel)
    }
  }
  lel <- levels(as.factor(df[,subtype]))
  colsubtype <- setNames(RColorBrewer::brewer.pal(8,"YlGnBu")[c(2,4,6,8)],lel)
  colist <- colist[cha]

  hatop <- HeatmapAnnotation(
    barplot = anno_barplot(df[,topbar],gp = gpar(fill = colR(df[,topbar]),lty=0),
                           axis = TRUE,border=F),
    annotation_height = unit(10, "mm"), #adjust the height
    gp=gpar(lty = "solid", lwd = 1,col="white"))
  ha <- HeatmapAnnotation(df[,cha],
                          na_col = "white", col=colist,
                          gap = unit(c(rep(1,length(cha)-1)), "mm"), #adjust the gap
                          annotation_height = unit(rep(4,length(cha)), "mm"), #adjust the height
                          gp=gpar(lty = "solid", lwd = 1,col="white"))
  hm <- Heatmap(t(as.matrix(df[,subtype])),
                top_annotation=hatop,
                bottom_annotation=ha,
                col=colsubtype,
                cluster_rows = FALSE,
                show_column_names = F,
                show_row_names = F)
  pic(outp,w,h)
  draw(hm,gap = unit(5, "mm"),
       padding = unit(c(25, 45, 25, 5), "mm"), #space around plot
       show_annotation_legend = T,
       heatmap_legend_side="right",
       annotation_legend_side="top")
  for(an in colnames(df[,cha])){
    decorate_annotation(an,grid.text(an, unit(-5,"mm"), just = "right"))
  }
  decorate_annotation("barplot",grid.text(topbar, unit(-5,"mm"), just = "right"))
  dev.off()
  invisible(NULL)
}


