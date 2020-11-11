#' Color related functions and variables
#'
#' A group of color related functions and variables.\cr
#' \itemize{
#'   \item darken,lighten: darken or lighten colors. Credit to \url{https://gist.github.com/Jfortin1}
#'   \item add.alpha: add alpha to color vector. Credit to \url{https://magesblog.com/post/2013-04-30-how-to-change-alpha-value-of-colours-in/}
#'   \item assignCHMcolor: produce a color list for complexheatmap from predefined base colors.
#'   }
#' @section Predefined Colors:
#' These variables will be loaded to the global environment and might overwrite existing colors.
#' \itemize{
#'   \item colpal: vector of 6 for groups
#'   \item colpal2: c("deepskyblue3","ivory","coral") #for heatmap
#'   \item comp_hm_colist_full: for complex heatmap related color
#'   \item colonco: for mutations on oncoprint plot
#'}
#' @param color a character or vector of color
#' @param factor a numeric. To what extend you want to change the color
#' @param alpha a numeric. The alpha level.
#' @param df dataframe to be plotted
#' @param colist list of named color vector. By default is comp_hm_colist_full
#' @param columns if only selected columns to be used for generating the color, supply a vector of column names here.
#' @return \itemize{
#'  \item darken,lighten: a character or vector of the changed color
#'  \item assignCHMcolor: a list of vectors (discrete colors) or functions (for continuous colors)
#'  }
#' @name cjbColor
NULL
#' @rdname cjbColor
#' @export
darken <- function(color, factor=1.4){
  sapply(color,function(color){
    col <- col2rgb(color)
    col <- col/factor
    col <- rgb(t(col), maxColorValue=255)
    col})
}
#' @rdname cjbColor
#' @export
lighten <- function(color, factor=1.4){
  sapply(color,function(color){
    col <- col2rgb(color)
    col <- col*factor
    col[col>255] <- 255
    col <- rgb(t(col), maxColorValue=255)
    col})
}
#' @rdname cjbColor
#' @export
add.alpha <- function(color, alpha=0.7){
  if(missing(color))
    stop("Please provide a vector of colours.")
  apply(sapply(color, col2rgb)/255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha=alpha))
}
#' @rdname cjbColor
#' @export
assignCHMcolor <- function(df,colist=comp_hm_colist_full,columns=NULL){
  # . . columns: the columns to be plotted.
  library(circlize)
  if(is.null(columns)) columns <- colnames(df)
  for (i in columns){
    if (is.numeric(df[,i])){
      #for continuous, take max,med and min for color assignment
      ma <- max(df[,i],na.rm=T)
      mi <- min(df[,i],na.rm=T)
      me <- median(df[,i],na.rm=T)
      if (i %in% names(colist)) colist[[i]] <- colorRamp2(c(mi,me,ma),colist[[i]])
      else colist[[i]] <- colorRamp2(c(mi,me,ma),colist[["cont"]])
      next
    }
    if (!i %in% names(colist)){
      lel <- levels(as.factor(df[,i]))
      if (length(lel)==2) colist[[i]] <- setNames(colist[["bool"]],lel)
      else colist[[i]] <- setNames(colist[["disc"]][1:length(lel)],lel)
    }
  }
  return(colist[columns])
}




#Define default color palette and shape styles ####

#' @export
colpal=setNames(c("#111111CF","#fd0d0dCF","#0fe00fEF","#090ee0CF","#4d1b7bCF","#f0d817CF"),as.character(1:6))
#' @export
colpal2=c("deepskyblue3","ivory","coral") #for heatmap
#shapes=c(15,17,18,19)
#' @importFrom RColorBrewer brewer.pal
#' @export
# complex heatmap color list
comp_hm_colist_full <- list(
  Cohort=setNames(c("lemonchiffon","lightcyan"),c("Asian","Caucasian")),
  Smoker=c("No"="#E1FFE0","Yes"="lightsteelblue4"),
  Gender=c("Female"="mistyrose","Male"="skyblue3"),
  GenderSmoker=c("mistyrose","mistyrose3","skyblue","skyblue3"),
  Stage=setNames(RColorBrewer::brewer.pal(4,"PuRd"),c("I","II","III","IV")),
  Grade=setNames(c("gray",circlize::colorRamp2(c(0,0.5,1),RColorBrewer::brewer.pal(3, "PRGn"))(seq(0,1,0.25))),
                 c("Not specified","Poorly","Moderately to Poorly",
                   "Moderately","Well to Moderately","Well")),
  Histology=setNames(RColorBrewer::brewer.pal(9,"Set3"),
                     c("Acinar","Invasive mucinous","Lepidic","Micropapillary",
                       "Minimally invasive non-mucinous","Mixed","NOS","Papillary","Solid")),
  Signature.Groups=setNames(c("#2D8CC4","#FF4C49","#CBEF53"),c("Smoking","Aging","APOBEC")),
  RNA.Subtype=setNames(c("#FF615E","#F7A034","#CADB48","#208CBF","#00B8CF","#00E0B7"),
    c("TRU-a", "TRU", "TRU-b", "PI", "nonTRU", "PP")),
  GD=c("No"="#E1FFE0","Yes"="lightsteelblue4"),
  Age=RColorBrewer::brewer.pal(3, "YlGnBu"),
  Purity=RColorBrewer::brewer.pal(3, "Blues"),
  GII=RColorBrewer::brewer.pal(3, "BuPu"),
  TMB=RColorBrewer::brewer.pal(3,"GnBu"),
  Ploidy=RColorBrewer::brewer.pal(7, "OrRd")[c(1,3,5)],
  pLM=RColorBrewer::brewer.pal(6, "YlGnBu")[c(2,4,6)],
  Num.Clone=RColorBrewer::brewer.pal(6, "BuPu")[c(2,4,6)],
  Num.Driver=RColorBrewer::brewer.pal(6,"RdPu")[c(1,3,5)],
  bool=c("whitesmoke","lightsteelblue4"),
  cont=RColorBrewer::brewer.pal(6,"PuBuGn")[c(2,4,6)],
  disc=c(RColorBrewer::brewer.pal(8,"Pastel1"),"#111111CF","#fd0d0dCF","#0fe00fEF","#090ee0CF","#4d1b7bCF","#f0d817CF")
)
#' @importFrom RColorBrewer brewer.pal
#' @export
# oncoprint color list
colonco=RColorBrewer::brewer.pal(8,"Set1")
colonco = c("Missense_Mutation" = colonco[2], "Nonsense_Mutation" = colonco[4],
        "Nonstop_Mutation"=colonco[6],"Stop_Codon_Del"=colonco[7],
        "Splice_Site"=colonco[1], "Frame_Shift_Indel"=colonco[5],
        "In_Frame_Indel"=colonco[3],"Silent"="#989898","Translation_Start_Site"=colonco[8])
#' @export
# alter_fun for oncoprint plot
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#E8E8E8", col = NA))
  },
  Silent = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = colonco["Silent"], col = NA))
  },
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = colonco["Missense_Mutation"], col = NA))
  },
  Nonsense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = colonco["Nonsense_Mutation"], col = NA))
  },
  Splice_Site = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = colonco["Splice_Site"], col = NA))
  },
  Frame_Shift_Indel = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = colonco["Frame_Shift_Indel"], col = NA))
  },
  In_Frame_Indel = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = colonco["In_Frame_Indel"], col = NA))
  }
)

