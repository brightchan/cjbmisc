
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


get_deg_DESeq2 <- function(countMat,name.feature,name.group1,name.group2,
                           samID.grp1,samID.grp2,f.design=NULL,
                           metadf=NULL,samID.coln=NULL,
                           outpdir=NULL){

  # Input:
  # countMat: matrix of raw counts, genes in rows, samples in columns
  # name.feature,name.group1,name.group2: string, just the names
  # samID.grp1,samID.grp2: vectors of sample IDs (matching the countMat)
  # f.design: a string of the design such as "~Age+Gender". A design matrix and formula will be built based on this formula
  # metadf: dataframe that contains sample ID (samID.coln),

  # Output:
  # DESeq2 object?


  if(!is.null(f.design)){
  ## without a design formula


  }else if (!is.null(metadf)&!is.null(samID.coln)){
  ## with a design formula


  }else stop("Please provide metadata to 'metadf' with a column specifying sample ID in 'samID.coln' for the experiment design")


  # Save the result to outpdir
  if(!is.null(outpdir)){
    fn <- paste0(outpdir,"/",name.feature,"_",name.group1,"_vs_",name.group2,".tsv")
    write_tsv(xxx,fn)
  }

  return(xxx)

}

