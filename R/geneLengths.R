#' Retrieve gene lengths
#'
#' A group of functions to gain length of genes, exons and cds.
#' \itemize{
#'   \item hg19ExonLengths: Gain total exon lengths (bp)
#'   \item hg19cdsLengths: Gain total cds lengths (bp)
#' }
#' @param symbols vector of gene symbols. Will be converted into character internally by the function.
#' @return a vector of lengths
#' @seealso Adapted from \url{https://www.biostars.org/p/62583/}
#' @import GenomicFeatures
#' @import GenomicRanges
#' @import AnnotationDbi
#@import TxDb.Hsapiens.UCSC.hg19.knownGene
#@import org.Hs.eg.db
#' @name geneLengths
NULL
#' @rdname geneLengths
#' @export
# Get exon length based on gene symbols https://www.biostars.org/p/62583/
hg19ExonLengths <- function(symbols){
  require("org.Hs.eg.db")
  require("TxDb.Hsapiens.UCSC.hg19.knownGene")
  symbols <- as.character(symbols)
  exons.db = exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by='gene')
  egs    = unlist(mget(symbols[ symbols %in% keys(org.Hs.egSYMBOL2EG) ],org.Hs.egSYMBOL2EG) )
  sapply(egs,function(eg)
  {
    exons = exons.db[[eg]]
    exons = reduce(exons)
    sum( width(exons) )
  })
}
#' @rdname geneLengths
#' @export
# Get cds length based on gene symbols
hg19cdsLengths <- function(symbols){
  require("org.Hs.eg.db")
  require("TxDb.Hsapiens.UCSC.hg19.knownGene")
  symbols <- as.character(symbols)
  cds.db = cdsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by='gene')
  egs    = unlist(mget(symbols[ symbols %in% keys(org.Hs.egSYMBOL2EG) ],org.Hs.egSYMBOL2EG) )
  sapply(egs,function(eg)
  {
    cds = cds.db[[eg]]
    cds = reduce(cds)
    sum( width(cds) )
  })
}
