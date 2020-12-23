#' Survival plot with ggsurvplot
#' A warpper to plot survival plots
#' to further customize the plot,
#' @name ggsurvpWrap
#' @param survdf dataframe with \emph{surv.time} (character or numeric) and \emph{surv.status}  (numeric of 0,1; or factor with Alive as base level)
#' @param surv.time,surv.status column name of the survival time and status in the dataframe
#' @param fac character of the column name in the survdf to stratify the patients
#' @param cont.to.dis if \emph{fac} is numeric, how to discretize it. Can be a funciton (Default \emph{median}) or cutoff points (n-1 points)
#' @param col vector of colors
#' @param textsize pass to ggplot global text size
#' @param xlim,break.time.by,pval.method,risk.table,... arguments to \code{\link[survminer]{ggsurvplot}}
#' @return A object of class ggsurvplot
#' @import survminer
#' @import survival
#' @export
ggsurvpWrap <- function(survdf,fac,surv.time="OS_Month",surv.status="OS_Status",
                        cont.to.dis=NA,
                        title=fac,col=comp_hm_colist_full$disc,risk.table = F,
                        xlim=c(0,120),break.time.by=60,textsize=10,pval.method=T,...){
  # if(is.null(col)) col <- c("#111111CF","#fd0d0dCF","#0fe00fEF","#090ee0CF","#4d1b7bCF",
  #                           "#f0d817CF","#483698FF","#832424FF","#005800FF","#3A3B98FF",
  #                           "#602995FF","#574700FF")
  # avoid problems when using tibbles
  survdf <- as.data.frame(survdf)
  if(is.numeric(survdf[,fac])) {
    if(is.na(cont.to.dis)) cont.to.dis <- median
    survdf[,fac] <- discretize(survdf[,fac],cont.to.dis) #transform
    }
  columns <- c(surv.time,surv.status,fac)
  survdf <- survdf[complete.cases(survdf[,columns]),columns]
  survdf[,fac] <- as.factor(survdf[,fac])
  if(is.character(survdf[,surv.status])){
    warning(paste0(surv.status, " are characters. Converting to factors. Check the levels!"))
    survdf[,surv.status] <- as.numeric(as.factor(survdf[,surv.status]))-1
  }else if(!is.numeric(survdf[,surv.status])) survdf[,surv.status] <- as.numeric(survdf[,surv.status])-1
  fa <- as.formula(paste0("Surv(",surv.time,",",surv.status,")~as.factor(",fac,")"))
  fit <- surv_fit(fa,survdf)
  ggsurvplot(fit,data=survdf,
             title=title,palette = col[1:length(levels(survdf[,fac]))],
             legend.labs =levels(survdf[,fac]),legend.title="",
             #legend="none",
             risk.table = risk.table,conf.int = F ,pval = TRUE, xlab="Time",
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
