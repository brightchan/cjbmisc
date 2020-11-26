#' A group of functions for RNAseq related analysis.
#'
#' \itemize{
#'   \item normSizeF: perform DESeq2 normalization
#'   \item nmfPre: pre-clustering for NMF \code{\link[NMF]{nmf}} to determine the optimal rank
#'   \item nmfCluster: NMF \code{\link[NMF]{nmf}} clustering output the NMF object and the assigned subtype
#'   \item write_gct: transform a expression matrix (gene in rows with rowname) into gct format and save file with name "gctfn".
#'   \item write_cls: generate a cls format from a vector of subtype and save file with name "clsfn".
#' }


#' @name RNAseq
NULL
#' @rdname RNAseq
#' @export
normSizeF <- function(countdata) {# ... for the integration into pipeline
  library("DESeq2")
  # Construct single column coldata (for DESeq2)
  coldata <- data.frame(Sample=colnames(countdata))
  dds <- DESeqDataSetFromMatrix(countdata, coldata, ~1) #~1:no design
  dds <- estimateSizeFactors(dds)
  norm_counts <- log2(counts(dds, normalized=TRUE)+1) # add 1 pseudocount
  colnames(norm_counts) <- colnames(countdata)
  return(norm_counts)
}

#' @rdname RNAseq
#' @export
nmfPre <- function(data,rank=3:6,r=50){
  set.seed(12345)
  doParallel::registerDoParallel(cores=12)
  result <- NMF::nmf(data,rank,nrun=r,.opt="v")
}

#' @rdname RNAseq
#' @export
nmfCluster <- function(data,rank=3,r=200){
  set.seed(12345)
  doParallel::registerDoParallel(cores=6)
  result <- NMF::nmf(data,rank,nrun=r,.opt="v")
  clusters <-matrix(apply(coef(result), 2, which.max))
  return(list(res=result,subtype=clusters))
}


#' @rdname RNAseq
#' @export
write_gct <- function(exprdf,gctfn,genelist=NULL){
  if(is.null(genelist)) genelist <- row.names(exprdf)
  write(paste0(c("#1.2",rep("",ncol(exprdf)+1)),collapse ="\t"),file=gctfn)
  write(paste(c(length(genelist),ncol(exprdf),rep("",ncol(exprdf))), collapse ="\t"),
        file=gctfn,append = T,ncolumns =2)
  write(paste(c("NAME","Description",colnames(exprdf)),collapse = "\t"),
        file=gctfn,append = T,ncolumns =ncol(exprdf)+2)
  write.table(cbind(genelist,rep(NA,length(genelist)),
                    exprdf),
              file = gctfn,append = T,
              quote = F,sep = "\t",row.names=F,col.names=F)
}

#' @rdname RNAseq
#' @export
write_cls <- function(subtype,clsfn){
  fileConn<-file(clsfn,"wb")
  writeLines(c(paste0(c(length(subtype),length(unique(subtype)),"1"),collapse=" "),
               paste0("# ",paste0(1:length(unique(subtype)),collapse=" ")),
               paste(as.numeric(as.factor(subtype)),collapse=" ")), fileConn)
  close(fileConn)
}





gageAna <- function(res,geneset=kegg.sets.hs,q=0.01){
  foldchanges = res$log2FoldChange*(1-padj)
  names(foldchanges) = res$entrezgene
  foldchanges <- foldchanges[!is.na(names(foldchanges))]
  foldchanges <- foldchanges[order(abs(foldchanges),decreasing = T)]
  foldchanges <- foldchanges[!duplicated(foldchanges)]
  keggres = gage::gage(foldchanges, gsets=geneset, same.dir=T,
                       compare = "unpaired",use.stouffer=TRUE)
  up <- na.omit(keggres$greater[,"q.val"])
  up <- up[up<q]
  #substring pathway names into two parts
  up <- list(id=names(up),path.name=names(up),q.val=up)
  down <- na.omit(keggres$less[,"q.val"])
  down <- down[down<q]
  down <- list(id=names(down),path.name=names(down),q.val=down)
  message("up: ",paste0(up$id,"; \n"))
  message("down: ",paste0(down$id,"; \n"))
  list(up=up,down=down)
}

plot_pathview = function(res,pid,title) {
  genes <- res$log2FoldChange
  names(genes) <- row.names(res)
  pathview(gene.data=genes, pathway.id=pid,
           species="hsa", out.suffix=title, kegg.dir=paste0(outp,"/Pathview/"))
}

