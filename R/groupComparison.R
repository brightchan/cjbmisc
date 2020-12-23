
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
      group1=map_chr(feat.comb,1),
      group2=map_chr(feat.comb,2))
  }

  # get the sample list from the group names
  df.out <- pmap_dfr(df.para,get_sam4grp,metadf=metadf,sample.coln=sample.coln,sep.grp=sep.grp)

  return(df.out)
}


get_DegDESeq2 <- function(countMat,name.feature,name.group1,name.group2,
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

    fn.tsv <- paste0(resdir,"/S01_dge_",name.feature,"_",name.group1,"_vs_",name.group2,".tsv") # file name
    write_tsv(res.dds.df,fn.tsv)

    fn.col <- paste0(resdir,"/S01_dge_coldata_",name.feature,"_",name.group1,"_vs_",name.group2,".tsv")
    write_tsv(coldata,fn.col)

    system(paste0("touch ",resdir, "/Design_",as.character(design) %>% paste0(collapse = "")))

    if(!is.null(metadf)){

      fn.meta <- paste0(resdir,"/S01_dge_metadata_",name.feature,"_",name.group1,"_vs_",name.group2,".tsv")
      write_tsv(metadf,fn.meta)
    }

  }

  return(list(des=res.dds,df=res.dds.df))

}

fgseaWrap <- function(fn.gmt,ID.gene,stat.gene,minSize=20,maxSize=200,
                      fn.outp=NULL,...){
  pathways <-  fgsea::gmtPathways(fn.gmt)
  fgseaRes <-  fgsea(pathways=pathways, stats=scores %>% sort(.,decreasing = T),
                     minSize = minSize,maxSize = maxSize,...)

  if(!is.null(fn.outp)) write_tsv(fgseaRes,fn.outp)

  return(fgseaRes)
}



