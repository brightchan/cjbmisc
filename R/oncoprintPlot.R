#' Oncoprint plots
#'
#' A group of functions for oncoprint plots using simplied single functions.
#' \itemize{
#'   \item maf2oncoprintdf: Input a maf, generate the dataframe of gene mutation types for each patients as the input for oncoprint plot.
#' }
#' @param inpmaf dataframe of the input maf, should have the columns of Variant_Classification, Hugo_Symbol, and Tumor_Sample_Barcode
#' @param subtypeOrder named character of the subtype info of the patients. If provided, the oncodf would be ordered by this subtype info.
#' @param gene vector of characters for the symbols of genes to be plotted. Default all genes.
#' @param ordergene boolean for whether the genes should be orderd by the order of the gene list
#' @param coding boolean for whether to subset only the coding mutations
#' @param noSilent boolean for whether to remove the silent mutations
#' @param mergeVariant boolean for whether to merge the "Del and Ins" to "Indel" and various mutations to "Translation_Start_Site"
#' @param subtypeOrder a cha vector specifying the subtype of the samples, the oncodf will be ordered by subtype first.
#' @return a plot saved in desinated path
#' @seealso \code{\link[ComplexHeatmap]{oncoPrint}}
#' @name oncoprint
NULL
#' @rdname oncoprint
#' @export
#' @import ComplexHeatmap
#' @importFrom circlize colorRamp2
#' @importFrom RColorBrewer brewer.pal
#
maf2oncoprintdf <- function(inpmaf,gene=NULL,ordergene=FALSE,coding=TRUE,
                            noSilent=FALSE,mergeVariant=TRUE,subtypeOrder=NULL){
  # subset by gene name
  if(!is.null(gene)) inpmaf <- inpmaf[inpmaf$Hugo_Symbol%in%gene,]
  # subset only coding mutations
  if(coding){
    coding.variation <- c("Silent","Missense_Mutation","Nonsense_Mutation","Splice_Site",
                          "Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Ins","In_Frame_Del",
                          "De_novo_Start_InFrame","De_novo_Start_OutOfFrame",
                          "Start_Codon_Del","Start_Codon_SNP","Start_Codon_Ins",
                          "Nonstop_Mutation","Stop_Codon_Del")
    inpmaf <- inpmaf[inpmaf$Variant_Classification%in%coding.variation,]
  }
  if(noSilent){inpmaf <- inpmaf[inpmaf$Variant_Classification!="Silent",]}
  # merge similar variants
  if(mergeVariant){
    inpmaf$Variant_Classification <- gsub("Frame_Shift_Del|Frame_Shift_Ins","Frame_Shift_Indel",inpmaf$Variant_Classification)
    inpmaf$Variant_Classification <- gsub("In_Frame_Ins|In_Frame_Del","In_Frame_Indel",inpmaf$Variant_Classification)
    inpmaf$Variant_Classification <- gsub("De_novo_Start_InFrame|De_novo_Start_OutOfFrame|Start_Codon_Ins|Start_Codon_Del|Start_Codon_SNP",
                                          "Translation_Start_Site",inpmaf$Variant_Classification)
  }
  # generate the dataframe for oncoprint
  genes <- unique(inpmaf$Hugo_Symbol)
  oncodf <- matrix("",ncol=length(genes),
                   nrow=length(levels(inpmaf$Tumor_Sample_Barcode)))
  row.names(oncodf) <- levels(inpmaf$Tumor_Sample_Barcode)
  colnames(oncodf) <- genes
  for (i in 1:nrow(inpmaf)){
    g <- as.character(inpmaf[i,"Hugo_Symbol"])
    s <- as.character(inpmaf[i,"Tumor_Sample_Barcode"])
    v <- as.character(inpmaf[i,"Variant_Classification"])
    oncodf[s,g] <- paste0(oncodf[s,g],v,";")
  }
  # order genes
  if (ordergene) {
    oncodf <- oncodf[,gene[gene%in%colnames(oncodf)]]
  } else{
    # order genes by patient frequency of a gene
    oncodf <- oncodf[,names(sort(apply(oncodf,2,function(x) length(grep(";",x,value = F))),decreasing = T))]
  }
  # order patients
  if (is.null(subtypeOrder)){
    #Order by mutations
    seq <- rev(do.call(order,as.data.frame(oncodf)))
    oncodf <- oncodf[seq,]
  } else{
    #order by subtype first then mutations. subset only patients in the subtype
    oncodf <- oncodf[row.names(oncodf)%in%names(subtypeOrder),]
    s <- rev(do.call(order,data.frame(
      Subtype=subtypeOrder[row.names(oncodf)],
      gsub(";.*$","",as.matrix(oncodf)))))
    oncodf <- oncodf[s,]
  }
  return(oncodf)
}
