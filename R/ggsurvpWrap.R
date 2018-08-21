#' Survival plot with ggsurvplot
#' A warpper to plot survival plots
#' to further customize the plot,
#' @name ggsurvpWrap
#' @param survdf dataframe with "OS_Month" (character or numeric) and "OS_Status" (numeric of 0,1; or factor with Alive as base level)
#' @param fac character of the column name in the survdf to stratify the patients
#' @param discretize.method if \emph{fac} is numeric, how to discretizeit. Can be funcitons (Default median) or cutoff points (n-1 points)
#' @param col vector of colors
#' @param ... arguments to ggsurvplot
#' @return  A object of class ggsurvplot
#' @seealso  \code{\link[survminer]{ggsurvplot}}
#' @import survminer
#' @import survival
#' @export
ggsurvpWrap <- function(survdf,fac,title=NULL,col=comp_hm_colist_full$disc,risk.table = F,discretize.method=median,
                        xlim=c(0,120),break.time.by=60,textsize=10,pval.method=T,...){
  # if(is.null(col)) col <- c("#111111CF","#fd0d0dCF","#0fe00fEF","#090ee0CF","#4d1b7bCF",
  #                           "#f0d817CF","#483698FF","#832424FF","#005800FF","#3A3B98FF",
  #                           "#602995FF","#574700FF")
  if(is.numeric(survdf[,fac])) survdf[,fac] <- discretize(survdf[,fac],cont2dis) #transform
  columns <- c("OS_Month","OS_Status",fac)
  survdf <- survdf[complete.cases(survdf[,columns]),columns]
  survdf[,fac] <- as.factor(survdf[,fac])
  if(is.character(survdf[,"OS_Status"])){
    warning("OS_Status are characters. Converting to factors. Check the levels!")
    survdf[,"OS_Status"] <- as.numeric(as.factor(survdf[,"OS_Status"]))-1
  }else if(!is.numeric(survdf[,"OS_Status"])) survdf[,"OS_Status"] <- as.numeric(survdf[,"OS_Status"])-1
  fa <- as.formula(paste0("Surv(OS_Month,OS_Status)~as.factor(",fac,")"))
  fit <- surv_fit(fa,survdf)
  #print(survdiff(fa,survdf))
  #fit <- surv_fit(Surv(OS_Month,OS_Status)~get(fac),survdf)
  ggsurvplot(fit,data=survdf,
             title=title,palette = col[1:length(levels(survdf[,fac]))],
             legend.labs =levels(survdf[,fac]),legend.title="",
             #legend="none",
             risk.table = risk.table,conf.int = F ,pval = TRUE, xlab="Time (months)",
             ggtheme=
               theme(axis.title.y=element_text(size=textsize,face="plain",margin=margin(0,10,0,0)),
                     #axis.title.x=element_text(size=20,face="plain"),
                     #axis.title = element_text(size=18,face="plain"),
                     text = element_text(size=textsize,face="plain",color="black"),
                     #panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
                     panel.border = element_rect(linetype = "solid",fill=NA),
                     #panel.border = element_blank(),
                     axis.line = element_line(colour = "black")),
             xlim=xlim,
             break.time.by=break.time.by,
             ...)
}
