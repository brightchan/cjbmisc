#' A group of functions for Differential gene expression and pathway analysis after that.
#'
#' \itemize{
#'   \item get_sam4grp: generate tibble of library IDs according to the groups to be compared for DEG analysis
#'   \item get_compdf: wrapper using get_sam4grp to generate all pairwise or one-vs-all comparisons
#'   \item get_deg_EdgeR: generate DEG output with EdgeR
#'   \item get_deg_Deseq2: generate DEG output with DESeq2
#'   \item fgseaWrap: wrapper to run fgsea
#'   \item DEG_fgsea_Wrap: wrapper to run DESeq2 and fgsea
#'  }

#' @rdname DEG_Pathway
#' @export
get_sam4grp <- function(metadf,feat.coln,group1, group2, sample.coln,sep.grp="+"){

  # feat.coln:  column name of a categorical value of the feature to be compared
  # group1,group2: string of the groups to be included in group1 and group2, separated by "sep.grp"
  # sample.coln: column name of sample ID
  # metadf: must contain feat.coln and sample.coln

  # output: a tibble with five columns: Feature, group1, group2, grp1.sample, grp2.sample
  # the last two columns are lists of chr
  grp1 <- str_split_fixed(group1,paste0("[",sep.grp,"]"),Inf) %>% unlist
  grp2 <- str_split_fixed(group2,paste0("[",sep.grp,"]"),Inf) %>% unlist
  samID.grp1 <- metadf %>% filter(.data[[feat.coln]]%in%grp1) %>% pull({{sample.coln}})
  samID.grp2 <- metadf %>% filter(.data[[feat.coln]]%in%grp2) %>% pull({{sample.coln}})


  return(tibble(name.feature=feat.coln,
                name.group1=group1,
                name.group2=group2,
                samID.grp1=list(samID.grp1),
                samID.grp2=list(samID.grp2)))
}


#' @rdname DEG_Pathway
#' @export
get_compdf <- function(metadf,feat.coln,sample.coln,comp="1vA",sep.grp="+"){

  # comp: comparisons to be generated, "1vA" is one vs All the rest, "pair" is all possible pairwise comparison
  # feat.coln:  column name of a categorical value of the feature to be compared
  # sample.coln: column name of sample ID
  # metadf: must contain feat.coln and sample.coln

  # return a df with five columns: Feature, group1, group2, grp1.sample, grp2.sample
  # the last two columns are lists of chr

  # get all unique features in the feat.coln
  feat <- metadf[[feat.coln]] %>% unique %>% na.omit()


  if (!comp%in%c("1vA","pair")){
    stop('"comp" must be one of "1vA","pair"')
  }else if (comp=="1vA"){
    # generate tibble of all possible combinations of one element vs all the rest
    df.para <- tibble(
      feat.coln=feat.coln,
      group1=feat
    ) %>% rowwise %>%
      mutate(group2=paste0(setdiff(feat,group1),
                           collapse = sep.grp)) %>%
      ungroup()

  }else if (comp=="pair"){
    # generate all combinations of the features
    feat.comb <- combn(feat,2,simplify = F)

    df.para <- tibble(
      feat.coln=feat.coln,
      group1=sapply(feat.comb,"[",1),
      group2=sapply(feat.comb,"[",2))
  }

  # get the sample list from the group names
  df.out <- pmap_dfr(df.para,get_sam4grp,metadf=metadf,sample.coln=sample.coln,sep.grp=sep.grp)

  return(df.out)
}


#' @rdname DEG_Pathway
#' @export
fgseaWrap <- function(fn.gmt,ID.gene,stat.gene,nperm=1000,
                      minSize=20,maxSize=500,n.top=10,
                      fn.outp=NULL,dir.outp=NULL,outp.return=F){

  # "fn.gmt" is the path to gmt file
  # ID.gene = character, vector of genes
  # stat.gene = the log2 fold change or pvalues etc that is used to do the fgsea

  stat.gene=setNames(stat.gene,ID.gene)
  pathways<-fgsea::gmtPathways(fn.gmt)

  fgseaRes<-fgsea(pathways=pathways, stats=stat.gene, nperm = nperm,
                  minSize = minSize,maxSize = maxSize)
  fgseaRes<-fgseaRes[order(pval), ] %>%
    filter(!is.na(pval))

  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=n.top), pathway] # top n.top up pathways
  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=n.top), pathway] # top n.top down
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

  # Saving result table and plot
  if(!is.null(fn.outp)&!is.null(dir.outp)&nrow(fgseaRes)>0){

    data.table::fwrite(fgseaRes, file=paste0(dir.outp,"/fgsea_",fn.outp,".tsv"),
                       sep="\t", sep2=c("", " ", ""))

    pdf(paste0(dir.outp,"/",fn.outp,".pdf"), width = 12, height = 8)
    plotGseaTable(pathways[topPathways], stat.gene, fgseaRes)
    dev.off()
  }

  if (outp.return) return(fgseaRes)
  else return(NULL)
}


#' @rdname DEG_Pathway
#' @export
get_deg_EdgeR <- function(countMat,name.feature,name.group1,name.group2,
                          samID.grp1,samID.grp2,f.design=NULL,
                          metadf=NULL,samID.coln=NULL,deseq.coln.trans=T,
                          df.gene.anno=NULL,
                          dir.outp=NULL,fn.outp=NULL){

  # countMat: rownames as geneID and colnames as sampleID
  # df.gene.anno: dataframe with one column name "gene" and other annotation columns

  # group1 is numerator, group2 is denominator
  name.feature <- make.names(name.feature)
  name.group1 <- make.names(name.group1)
  name.group2 <- make.names(name.group2)

  ### make original DGE list
  dge = DGEList(countMat[,c(samID.grp1,samID.grp2)],
                group=c(rep(name.group1,length(samID.grp1)),
                        rep(name.group2,length(samID.grp2))))

  ### calculate normalized factor
  dge = calcNormFactors(dge, method="TMM")


  ### define design matrix
  # Not that the feature of interest is the last one
  # the group2 would be the baseline group
  if(is.null(f.design)){
    f.design <- as.formula(paste0("~0+",name.feature))
  }else{
    f.design <- as.formula(paste0(f.design,"+",name.feature))
  }

  if(is.null(metadf)) {
    metadf <- data.frame(
      factor(c(rep(name.group1,length(samID.grp1)),
               rep(name.group2,length(samID.grp2))),
             # so that group2 would be the base
             levels = c(name.group2,name.group1))) %>%
      `colnames<-`(name.feature)
  } else {
    metadf <- metadf[metadf[[samID.coln]]%in%c(samID.grp1,samID.grp2),]
    metadf[[name.feature]]=factor(c(rep(name.group1,length(samID.grp1)),
                                    rep(name.group2,length(samID.grp2))),
                                  # so that group2 would be the base
                                  levels = c(name.group2,name.group1))
  }

  mat.design = model.matrix(f.design,metadf)

  # set the contrast to be group1-group2
  glm.contrast <- setNames(rep(0,ncol(mat.design)),colnames(mat.design))
  glm.contrast[paste0(name.feature,name.group1)] <- 1
  # if no-intercept design used, need to put -1 for the baseline values,
  # otherwise the baseline factor won't be in the matrix and no need to add that -1 .
  if(grepl("0.*\\+|\\+.*0",as.character(f.design)[2])) glm.contrast[paste0(name.feature,name.group2)] <- -1

  ### estimate the dispersion using GLM Robust model
  disp = estimateGLMRobustDisp(dge, mat.design)

  ### fit the model using GLM quasi-likelihood and set robust=T
  fit = glmQLFit(disp, mat.design, robust=TRUE)

  ### use glmQLFTest
  qlf = glmQLFTest(fit, contrast = glm.contrast)
  res = as_tibble(qlf$table)
  ## add gene symbol and FDR
  res$gene = rownames(countMat)
  res$FDR = p.adjust(res$PValue, method = 'BH')
  res <- res %>% select(gene,everything())


  ### Transform the column names to be the same as DESeq2
  if(deseq.coln.trans) res <- res %>% rename(log2FoldChange=logFC,baseMean=logCPM,padj=FDR)


  ### Add gene name annotations
  if(!is.null(df.gene.anno)){
    # remove duplicated gene ID
    df.gene.anno <- df.gene.anno %>% filter(!duplicated(gene))
    res <- left_join(res,df.gene.anno)
  }

  ### Save the results
  if(!is.null(dir.outp)){
    # create subfolder based on feature
    dir.create(dir.outp,recursive = T,showWarnings = F)
    # use default file name if not provided
    if(is.null(fn.outp))  fn.outp <- paste0(dir.outp,"/DEG.EdgeR--",name.feature,"--",name.group1,"--Vs--",name.group2,".tsv") # file name
    else fn.outp=paste0(dir.outp,"/",fn.outp)
    write_tsv(res,fn.outp)

  }

  return(res)

}



#' @rdname DEG_Pathway
#' @export
get_deg_Deseq2 <- function(countMat,name.feature,name.group1,name.group2,
                   samID.grp1,samID.grp2,f.design=NULL,
                   metadf=NULL,samID.coln=NULL,
                   outpdir=NULL){

  # Input:
  # countMat: data frame of raw counts, genes in rows, samples in columns, has 2 extra cols: 1st col = ensembl_id; 2nd col = gene symbol
  # name.feature,name.group1,name.group2: string, just the names
  # samID.grp1,samID.grp2: vectors of sample IDs (matching the countMat)
  # f.design: a string of the design such as "~Age+Gender". A design matrix and formula will be built based on this formula
  # metadf: dataframe that contains sample ID (samID.coln),

  # Output:
  # DESeq2 object?

  ##########################################
  # Generate coldata from df.para and design
  ##########################################

  ## df for group1
  grp1=data.frame(samID = samID.grp1,
                  group = rep(name.group1, length(samID.grp1)))
  ## df for group2
  grp2=data.frame(samID = samID.grp2,
                  group = rep(name.group2, length(samID.grp2)))
  ## Join df1 & df2
  grp12=rbind(grp1,grp2) %>%
    rename_at(vars(group),~name.feature) %>%
    rename_at(vars(samID),~samID.coln)

  if(is.null(f.design)){
    ## without a design formula

    coldata <- grp12
    design <- as.formula(paste0("~",name.feature))

  }else if (!is.null(metadf)&!is.null(samID.coln)){
    ## with a design formula

    coldata <- grp12 %>%
      left_join(metadf %>%
                  select(-all_of(name.feature)),
                by=samID.coln)

    design <- as.formula(f.design)


  }else stop("Please provide metadata to 'metadf' with a column specifying sample ID in 'samID.coln' for the experiment design")


  ##########################################
  # Load, filter, reorder countdata
  ##########################################
  # -- Reorder the col (sample ID) so it will be the same order with coldata

  countdata=countMat %>%
    select(-2,1,coldata[[samID.coln]]) %>% # remove the 2nd col containing gene symbol
    `rownames<-`(.[[1]]) %>%  # make the 1st col containing ensembl id as rownames
    select(coldata[[samID.coln]]) # order the countdata column

  all(colnames(countdata)==coldata$RNA_lib) # check if the order of sample ID between coldata and countdata is the same

  #####################################
  # Perform DGE
  #####################################
  # -- Check first if the is custom design

  message("Performing DGE using DESeq2 ...")
  message(paste0("Analysing ", nrow(countdata), " genes in ", ncol(countdata), " samples"))

  dds <- DESeqDataSetFromMatrix(countData= round(countdata),
                                colData = coldata,
                                design = design)
  dds <- DESeq(dds) # DGE
  contrast <- c(name.feature,name.group1,name.group2) # specify group that want to compare
  res.dds <- results(dds, contrast) # result
  res.dds <- res.dds[order(res.dds$padj),] # order the result table based on p-adjust
  head(res.dds)
  message("DGE analysis finished, saving results ...")

  ###########################################
  # Save the result
  ###########################################

  if(!is.null(outpdir)){

    dir.create(paste0(outpdir,"/",name.feature)) # create subfolder based on feature
    resdir <- paste0(outpdir,"/",name.feature)

    # R object
    fn.rds <- paste0(resdir,"/S01_dge_",name.feature,"_",name.group1,"_vs_",name.group2,".RDS") # file name
    saveRDS(res.dds,fn.rds)

    # data.frame

    sym=countMat %>%
      select(1:2) %>% `colnames<-`(c("Geneid","gene_name")) # extract the gene symbol and ensembl id from count matrix so I can put back the gene symbol in DGE result

    res.dds.df=as.data.frame(res.dds) %>%
      tibble::rownames_to_column(var="Geneid") %>%
      left_join(sym, by="Geneid") %>% # join with gene symbol
      select(Geneid,gene_name,everything())

    fn.tsv <- paste0(resdir,"/S01_dge_",name.feature,"_",name.group1,"_vs_",name.group2,".tsv")
    write_tsv(res.dds.df,fn.tsv) # save dge result

    fn.col <- paste0(resdir,"/S01_dge_coldata_",name.feature,"_",name.group1,"_vs_",name.group2,".tsv")
    write_tsv(coldata,fn.col) # save coldata

    system(paste0("touch ",resdir, "/Design_",as.character(design) %>% paste0(collapse = ""))) # save design

    if(!is.null(metadf)){

      fn.meta <- paste0(resdir,"/S01_dge_metadata_",name.feature,"_",name.group1,"_vs_",name.group2,".tsv")
      write_tsv(metadf,fn.meta) # save metadata
    }

  }

  return(list(des=res.dds,df=res.dds.df))

}



#' @rdname DEG_Pathway
#' @export
DEG_fgsea_Wrap <- function(countMat,name.feature,name.group1,name.group2,
                        samID.grp1,samID.grp2,fn.gmt,f.design=NULL,
                        metadf=NULL,samID.coln=NULL,outpdir=NULL){

  deseqRes=getDGE(countMat=countMat,
                  name.feature=name.feature,
                  name.group1=name.group1,
                  name.group2=name.group2,
                  samID.grp1=samID.grp1,
                  samID.grp2=samID.grp2,
                  f.design=f.design,
                  metadf=master,
                  samID.coln=samID.coln,
                  outpdir=outpdir)


  # Make a new col for rank in deseq result
  dds.res=deseqRes$df %>%
    mutate(rank=-log10(.$padj)*sign(.$log2FoldChange)) %>%
    filter(!is.na(rank),!is.na(gene_name))

  fgseaRes=fgseaWrap(fn.gmt = fn.gmt,
                     ID.gene = dds.res$gene_name,
                     stat.gene = dds.res$rank,
                     fn.outp = paste0(name.feature,"_",name.group1,"_vs_",name.group2),
                     outpdir = paste0(outpdir,"/",name.feature))

  return(list("des"=deseqRes,"fgsea"=fgseaRes))

}
