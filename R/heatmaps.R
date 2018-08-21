#' Heatmaps
#' A group of functions for heatmaps
#' FIXME
#' 
#' 
#Heatmap of subset genes with subtypes. Input all in data.frame with match row/col names. col/rowdata use samples in row/col, filled with color name
htmap <- function(countdata,coldata,title,#rowdata=NULL,
                  breaks=seq(-1.8,1.8,length.out=201)^3+0.5) { # Adjust this color breaks for better look 
  countdata <-as.matrix(countdata[,order(coldata[,1])])
  coldata <- coldata[order(coldata[,1]),]
  #csc <- c("deepskyblue3","springgreen","coral")[as.factor(sort(coldata[,1]))]
  my_palette <- colorRampPalette(colpal2)(200)
  library(heatmap3)
  #print(rowdata)
  heatmap3(countdata,Colv=NA,#RowSideColors=rowdata,
           ColSideColors=coldata,breaks=breaks,
           col=my_palette,
           labCol = NA, labRow = NA,main=title)
}

#Heatmap using showAnn, which take only numeric and factor
htmap2 <- function(countdata,coldata,title,#rowdata=NULL,
                   breaks=seq(-1.8,1.8,length.out=201)^3+0.5) { # Adjust this color breaks for better look 
  countdata <-as.matrix(countdata[,order(coldata[,1])])
  coldata <- coldata[order(coldata[,1]),]
  #csc <- c("deepskyblue3","springgreen","coral")[as.factor(sort(coldata[,1]))]
  my_palette <- colorRampPalette(colpal2)(200)
  library(heatmap3)
  #print(rowdata)
  heatmap3(countdata,Colv=NA,#RowSideColors=rowdata,
           ColSideAnn=coldata,ColSideFun= function(x) showAnn(x),
           breaks=breaks,ColSideWidth=1.5,
           col=my_palette,labCol = NA, labRow = NA,main=title)
}


htmap3 <- function(countdata,coldata,title,#rowdata=NULL,
                   breaks=seq(-1.8,1.8,length.out=201)^3+0.5) { # Adjust this color breaks for better look 
  countdata <-as.matrix(countdata[,order(coldata[,1])])
  coldata <- coldata[order(coldata[,1]),]
  #csc <- c("deepskyblue3","springgreen","coral")[as.factor(sort(coldata[,1]))]
  my_palette <- colorRampPalette(c("deepskyblue3","ivory","coral"))(200)
  library(heatmap3)
  #print(rowdata)
  heatmap3(countdata,Colv=F,#RowSideColors=rowdata,
           ColSideColors=coldata,breaks=breaks,
           col=my_palette,hclustfun=hclust,
           labCol = NA, labRow = NA,main=title)
}
