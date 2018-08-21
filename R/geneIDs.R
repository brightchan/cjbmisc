#' Annotate gene symbol, geneID or EnsemblID etc
#'
#' FIXME
#'
#' A group of functions to annotate geneIDs
#' \itemize{
#'   \item del_dupcol: remove duplicated columns by name.
#'   \item del_alldupcol: remove columns that are not unique by name
#'   \item subcol_tcgaTnN: Subset tumor/normal from tcgacount
#'   \item subset_tcgasam: Subset one dataframe according to another, by colnames of patient name.
#' }
#' @return dataframe
#' @param df,x,y dataframe
#' @param tum TRUE for tumor and FALSE for normal
#' @name geneIDs
NULL
#' @rdname geneIDs
#' @export
#'
# Annotate gene info of this STAR-RSEM pipeline, since genes used are the same btw TCGA and GSK
EN2other <- function(x){
  library("mygene")
  EnID <- gsub("\\..*$","",row.names(x))
  EAnno <- queryMany(EnID, scopes="ensembl.gene",
                     fields=c("symbol","name","ensembl.gene","entrezgene","HGNC","type_of_gene"),
                     species="human") # DF with query,_id,notfound
  EAnno <- EAnno[!(duplicated(EAnno$query)),]
  return(as.data.frame(EAnno))
}
#EnAnno <- EN2other(srtum)#Store the anno of this gene set #use this EnAnno as reference

# Return annotation DF of query ENSEMBL as vector
ENquery <- function(x,EnAn=EnAnno){
  EnID <- gsub("\\..*$","",x)
  QAnno <- EnAn[match(EnID,EnAn$query),c("query","symbol","name","entrezgene","HGNC","type_of_gene")]
  return(as.data.frame(QAnno))
}

# Return annotation DF of query ENSEMBL as vector
Symb2other <- function(x){
  library("mygene")
  Anno <- queryMany(x, scopes="symbol",fields=c("symbol","name","ensembl.gene","entrezgene","HGNC","type_of_gene"),species="human") # DF with query,_id,notfound
  Anno <- Anno[!(duplicated(Anno$query)),]
  return(as.data.frame(Anno))
}

# Change rowname from ENSEMBL to geneID. Omit unmatched row. Sum and Delete rows with same geneID.
# Input dataframe output the same.
rname_EN2gID <- function(x,EngID=EnAnno) {
  library("mygene")
  #EnID <- gsub("\\..*$","",row.names(x))
  #temp <- queryMany(EnID, fields="entrezgene",species="human") # DF with query,_id,notfound
  #temp <- temp[!(duplicated(temp$query)),] #delete duplicated results
  x <- x[!is.na(EngID$"entrezgene"),]
  gID <- EngID[!is.na(EngID$"entrezgene"),"entrezgene"]
  x <- cbind(gID,x)
  dup_gID <- which(duplicated(gID)) # Row No. of duplicated w/o first element
  while (length(dup_gID)>0) {
    a <- which(gID==gID[dup_gID[1]]) # All row No. of duplicated
    x <- sumNdelete(x,a)
    dup_gID <- dup_gID[!dup_gID %in% a]
  }
  x <- na.omit(x)
  row.names(x) <- x[,1]
  x <- x[,-1]
  tempn <- sum(na.omit(EngID$notfound=="TRUE"))
  message(as.character(tempn)," rows ommited without hit. Data dim: ",nrow(x),"x",ncol(x))
  #print(temp[complete.cases(temp$notfound),"query"]) # Echo the ENSEMBLID of non-hit
  return(x)
}


# Change rowname from ENSEMBL to symbol.
rname_EN2sym <- function(x,Ensym=EnAnno) {
  row.names(x) <- gsub("\\..*$","",row.names(x))
  x <- x[!is.na(Ensym$"symbol"),]
  sym <- Ensym[!is.na(Ensym$"symbol"),"symbol"]
  x <- cbind(sym,x)
  dup_sym <- which(duplicated(sym)) # Row No. of duplicated w/o first element
  while (length(dup_sym)>0) {
    a <- which(sym==sym[dup_sym[1]]) # All row No. of duplicated
    x <- sumNdelete(x,a)
    dup_sym <- dup_sym[!dup_sym %in% a]
  }
  x <- na.omit(x)
  row.names(x) <- x[,1]
  x <- x[,-1]
  tempn <- sum(na.omit(Ensym$notfound=="TRUE"))
  message(as.character(tempn)," rows ommited without hit. Data dim: ",nrow(x),"x",ncol(x))
  #print(temp[complete.cases(temp$notfound),"query"]) # Echo the ENSEMBLID of non-hit
  return(x)
}

# Change rowname from symbol to geneID. Input dataframe output the same.
rname_symb2gID <- function(x){
  library("mygene")
  temp <- queryMany(row.names(x),scopes="symbol", fields="entrezgene",species="human") # DF with query,_id,notfound
  temp <- temp[! duplicated(temp$query),] # remove duplicated gene (some symbol match to multiple id)
  x <- x[complete.cases(temp$"entrezgene"),] # omit non-hit genes
  row.names(x) <- temp[complete.cases(temp$"entrezgene"),"entrezgene"]
  tempn <- sum(na.omit(temp$notfound=="TRUE"))
  message(as.character(tempn)," rows ommited without hit.")
  #print(temp[complete.cases(temp$notfound),"query"]) # Echo the ENSEMBLID of non-hit
  return(x)
}
