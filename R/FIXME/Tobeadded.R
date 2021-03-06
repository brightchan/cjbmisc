# New features:
Update plots with ggstats
Add pca functions with multiple subtype colors in separate panels.


# To update:

p_fish.chi.t <- function(df,v1,v2,alt="two.sided",p.test="both",ws=2e6){
  # if any subtype has only one member, return 1 as p value
  if (length(unique(df[,v1]))==1|length(unique(df[,v2]))==1) {
    message("Less than 1 category. Returning p.value 1.")
	return(1)
  }
  if(p.test=="both"){
    tryCatch(fish.t(df,v1,v2,alt,ws)$p.value,
             error=function(e){
               warning(e)
               warning("Using ChiSQ test instead.")
               chisq.t(df,v1,v2)$p.value})
  }else if (p.test=="fisher"){
    fish.t(df,v1,v2,alt,ws)$p.value
  }else if (p.test=="chi"){
    chisq.t(df,v1,v2)$p.value
  }else stop('p.test must be "fisher" or "chi"')
}

subtypeChaHM_col <- function(df,subtype,subtype.col=NULL,cha,topbar,order_by=topbar,
                             colist=comp_hm_colist_full,
                             pic,outp,w,h){
  if(is.null(subtype.col)){
    lel <- levels(as.factor(df[,subtype]))
    colsubtype <- setNames(comp_hm_colist_full$disc[1:length(lel)],lel)
    } else colsubtype <- subtype.col
  df <- df[order(df[,order_by]),]
  colR <- circlize::colorRamp2(breaks = c(min(df[,topbar],na.rm=T),
                                median(df[,topbar],na.rm=T),
                                max(df[,topbar],na.rm=T)),
                     colors = subtype.col)
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
  colsubtype <- setNames(subtype.col,lel)
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
  tryCatch(
    { pic(outp,w,h)
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
    },
    error=function(e){
      warning(e)
      warning("Plotting Aborted.")
      dev.off()}
  )
  invisible(NULL)
}

# To add


api_parse <- function(request){
  text <- content(reques,"text",encoding="UTF-8")
  if(identical(text,""))warning("No output to parse")
  fromJSON(text)
}




