#' Venn diagram plotting functions
#'
#' Plot the venn diagram.
#' @param l list of vectors to be drawn. for named list the names would be the label
#' @param n alternative label of the names
#' @param title main title of the diagram
#' @param filename if not NULL, will be saved as a file (tiff); if NULL, the plot will be printed to the default device
#' @param fill character vector giving the colour of each circle's area
#' @param dist numeric vector giving the distance (in npc units) of each category name from the edge of the circle (can be negative)
#' @param mar numeric giving the amount of whitespace around the diagram in grid units
#' @param ... pass to \code{\link[VennDiagram]{venn.diagram}}
#' @name venn
#' @export
#' @import VennDiagram
#' @importFrom grid grid.draw
#Venn Diagram
venn <- function(l,n=NULL,title=NULL,filename=NULL,fill=NULL,dist=0.08,mar=0.1,...){
  if (!is.null(n)) l <- setNames(l,n)
  if (!is.null(filename)){
    venn.diagram(l,main=title,filename = filename,
                 fill=fill,cat.dist=dist, margin=mar,...)
  } else {
    par(mar=rep(0,4))
    plot.new()
    grid.draw(venn.diagram(l,main=title,filename = filename,
                           fill=fill,cat.dist=dist, margin=mar,...))
  }
}
