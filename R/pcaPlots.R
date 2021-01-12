#' PCA Plots
#'
#' A group of functions for PCA plots.
#' \itemize{
#'   \item pca: PCA plot with one or two group info.
#'   \item pca_multi: Multiple PCA plots for top PCs.
#'   \item pca_bg: PCA plot with color background shading for each subtype
#' }
#' @param df numeric dataframe with rows of samples and columns of features
#' @param A,B,subtype character vector of the subtypes of samples
#' @param col color to be used for each subtype.
#' @param plot.it if FALSE, no plot printed, only return the plot object
#' @param n for pca_multi, number of PCs to be plotted
#' @param mar for pca_bg, margin between the points to the frame.
#' @param "+" use + to add ggplot layers for more customization
#' @return
#' Print the plot and return the plot
#' \itemize{
#'   \item pca: a ggplot object
#'   \item pca_multi: a plot object
#' }
#' @seealso \code{\link[ggplot2]{autoplot}} and \code{\link{prcomp}}
#' @name PCAPlots
NULL
#' @rdname PCAPlots
#' @export
#' @import ggplot2
#' @import ggfortify
# PCA plot. Input data.frame of countdata, vector of subtype and diff as in coldata for the htmap

pca <- function(df,A,B=NULL,title="PCA",
                col=c("#111111CF", "#fd0d0dCF", "#0fe00fEF", "#090ee0CF",
                      "#4d1b7bCF", "#f0d817CF"),
                nameA="subtypeA",nameB="subtypeB",
                point.size=4,shapes=c(16:18,15),
                label.points=FALSE,lab=NA,
                plot.it=TRUE,return.data=FALSE){
  pr <- prcomp(df)
  perc <- round((pr$sdev)^2 / sum(pr$sdev^2)*100,digits=2) #Percentage of each PC
  
  if(is.null(B)){
    anno <- data.frame(pr$x[,c(1:2)],A,lab=lab)
    colnames(anno)[3] <- nameA
    p <- ggplot(data=anno,aes_string(x="PC1",y="PC2",
                                     colour=nameA,label="lab"))+
      geom_point(size=point.size,shapes=shapes[1])
  }
  else{
    anno <- data.frame(pr$x[,c(1:2)],A,B,lab=lab)
    colnames(anno)[3:4] <- c(nameA,nameB)
    p <- ggplot(data=anno,aes_string(x="PC1",y="PC2",
                                     colour=nameA,label="lab"))+
      geom_point(colour=nameA,shape=nameB,
                 size=point.size)+
      scale_shape_manual(values=shapes)
  }
  
  p <- p+
    ggtitle(title)+
    xlab(paste0("PC1 (",perc[1],"%)"))+
    ylab(paste0("PC2 (",perc[2],"%)"))+
    scale_colour_manual(values=unname(col))
  
  
  if(label.points&!is.na(lab)){
    p <- p+geom_text_repel()
  }
  
  if(plot.it) print(p)
  if(return.data)return(list(p=p,pr=pr,perc=perc))
  return(p)
}






#' @rdname PCAPlots
#' @export
# Multiple PCA plots for top PCs ####
pca_multi <- function(df,subtype=NULL,n=6,col=colpal,title=NULL,plot.it=TRUE){
  pr=prcomp(df)
  if(n>ncol(pr$x)) {n=ncol(pr$x);message("n was set too large. Changed to number of maximum PCs.")}
  comp=data.frame(pr$x[,1:n])
  if(!is.null(subtype)){
    coldf <- col[as.numeric(as.factor(subtype))]
    p <- plot(comp,pch=20,col=coldf,main=title)
  } else plot(comp,main=title)
  if(plot.it) print(p)
  return(p)
}
#' @rdname PCAPlots
#' @export
# PCA plot with density as background ####
pca_bg <- function(df,subtype,col=colpal,legend.name=NULL,
                   point.size=6,text.size=26,mar=15,plot.it=TRUE){
  #1. prepare data
  subtype=as.character(subtype)
  mem.subtype=setNames(unique(subtype),unique(subtype))
  if(is.null(names(col)))col=setNames(col,mem.subtype)
  pr <- prcomp(t(df))
  df <- data.frame(x=pr$x[,1],y=pr$x[,2],
                   Subtype=subtype)
  nx <- round((max(df$x)-min(df$x)+2*mar)/5) #get the range of plot, add some margin
  ny <- round((max(df$y)-min(df$y)+2*mar)/5)
  pseup <- data.frame(x=rep(seq(min(df$x)-mar,max(df$x)+mar,length.out=nx),ny), #generate pseudo points within the range
                      y=rep(seq(min(df$y)-mar,max(df$y)+mar,length.out=ny),each=nx))
  #2. assign points to subtype based on euclidean distance
  cent <- lapply(mem.subtype,function(x){
    c(mean(df[df$Subtype==x,1]),mean(df[df$Subtype==x,2]))
  })
  centdf <- do.call(rbind,cent)
  centdist <- apply(pseup,1,function(i){
    mem.subtype[which.min(apply(centdf,1,function(x) dist_euc(x,i)))]
  })
  pseup$subtype <- centdist
  pseup.list <- lapply(mem.subtype,function(x){
    pseup[pseup$subtype==x,]
  })

  #3. plot as below
  perc <- round((pr$sdev)^2 / sum(pr$sdev^2)*100,digits=2) #Percentage of each PC
  p <- ggplot()+
    geom_point(data=df,aes(x,y,color=Subtype),
               shape=16,size=point.size,alpha=0.6)+
    scale_alpha(range = c(0.001, 0.25),guide=F)+# adjust shade intensity
    scale_fill_identity()+
    scale_color_manual(values=col,name=legend.name)+
    xlab(paste0("PC1 (",perc[1],"%)"))+
    ylab(paste0("PC2 (",perc[2],"%)"))+
    theme_minimal()+coord_cartesian(expand = FALSE)+
    theme(text=element_text(size=text.size),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size=text.size,color="black"),
          axis.title = element_text(size=text.size,face="plain"),
          axis.title.x = element_text(margin=margin(20,10,10,10)),
          axis.text.x = element_text(vjust=2),
          legend.position = "none",
          aspect.ratio=1,
          panel.border = element_rect(colour = "black", fill=NA, size=1)
          #legend.title = element_text(size=12,face="bold"),
          #legend.text = element_text(size=10,face="bold"),
    )

  for(i in mem.subtype){
    #Add shade according to pseudo points
    p <- p+
      stat_density_2d(geom = "raster",data=pseup.list[[i]],
                      aes(x,y,alpha=..density..),fill=col[i],contour = FALSE)
  }
  if(plot.it) print(p)
  return(p)
}


