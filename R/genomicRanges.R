#' Genomic Ranges Functions
#'
#' A group of functions for genomic ranges
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
#'
olpByPerc <- function(gr1,gr2,thres){ # find overlapping regions by thres
  hits <- findOverlaps(gr1,gr2)
  olp <- pintersect(gr1[queryHits(hits)],
                    gr2[subjectHits(hits)])
  plop <- width(olp)/min(width(gr1[queryHits(hits)]),width(gr2[subjectHits(hits)]))
  hits <- hits[plop>thres]
}
