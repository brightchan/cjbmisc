% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getPvalue.R
\name{getPvalue}
\alias{getPvalue}
\alias{p_fish.chi.t}
\alias{p_aov.t}
\alias{p_ContCont}
\alias{p_ContDisc}
\alias{p_lm}
\alias{p_xVsAll}
\alias{p_dfAll}
\alias{p_uniCox}
\alias{p_adjust_mat}
\title{Returning only the p values}
\usage{
p_fish.chi.t(df, v1, v2, alt = "two.sided", p.test = "both", ws = 2e+06)

p_aov.t(df, v1, v2)

p_ContCont(df, v1, v2, method = "spearman")

p_ContDisc(df, v1, v2, method = "kruskal.test")

p_lm(df, v1, v2)

p_xVsAll(
  df,
  x.coln,
  y.coln = NULL,
  num.num.test = "spearman",
  cat.num.test = "kruskal.test",
  cat.cat.test = "both"
)

p_dfAll(df, x.coln, y.coln = NULL, ...)

p_uniCox(
  survdf,
  feat,
  surv.time = "PFS.days",
  surv.status = "PFS.status",
  signif.cutoff = 0.05,
  keep.all.in.barplot = F,
  plot.surv = T,
  survp.xlim = c(0, 1400),
  survp.timebreak = 365
)

p_adjust_mat(pvaldf, p.adjust.method = "BH")
}
\arguments{
\item{df}{dataframe with rows of samples and columns of features, one column should contain the subtype info.}

\item{v1, v2}{character of the column name in the df}

\item{alt}{passing to alternative in \code{\link[stats]{fisher.test}}, should be one of "two.sided", "greater" or "less". You can specify just the initial letter. Only used in the 2 by 2 case.}

\item{p.test}{Should be "fisher" or "chi", to use only that method to get p values.}

\item{ws}{passing to workspace size in \code{\link[stats]{fisher.test}},  integer specifying the size of the workspace used in the network algorithm. In units of 4 bytes.}

\item{method}{p_ContDisc: the test method to-be-used, can be t.test, kruskal.test or wilcox.test}

\item{...}{for p_feat_subtype, pass to \code{\link[cjbmisc]{p_fish.chi.t}}}

\item{feat}{character vector of the features to be calculated}

\item{subtype}{column name of the subtype}

\item{cont.test}{the significance test to be used for continuous variables. Use the name of the r basic tests.}

\item{disc.test}{the significance test to be used for discrete variables. Should be "fisher","chi" or"both"}
}
\description{
A group of functions to return only the p values from various of tests.
\itemize{
  \item p_fish.chi.t: Try first use \code{\link[cjbmisc]{fish.t}} to get p value,\cr
  if failed (usually due to too many catergories), then use \code{\link[cjbmisc]{chisq.t}} with popup warning. \cr
  If only Fisher's exact test or ChiSQ is wanted, specify it in "p.test" parameter
  \item p_aov.t: p-value from \code{\link[stats]{aov}}
  \item p_ContDisc: p-value from \code{\link[stats]{t.test}},\code{\link[stats]{wilcox.test}},\code{\link[stats]{kruskal.test}}
  \item p_ContCont: p-value from \code{\link[stats]{cor.test}} or from\code{\link[cjbmisc]{p_lm}}
  \item p_lm: p-value of the coefficient from the univariate \code{\link[stats]{lm}}
  \item p_adjust_mat: adjust a matrix of p-values
  \item p_xVsAll: generate a vector of pvalues by contrasting x vs all the rest of the variables in a dataframe.
  \item p_dfAll: generate a df of pvalues by contrasting all x.coln vs all y.coln. parameters are the same as in p_xVsAll
  \item p_uniCox: generate a list of univariate Cox models and a vector of the minimal pvalue from each model
  \item p_feat_subtype: p-value of a set of features according to a subtype in a manner similar to \code{\link[cjbmisc]{plot_feat_subtype}}
}
}
\seealso{
\code{\link{chisq.t}}, \code{\link{fish.t}}, \code{\link[stats]{fisher.test}},\code{\link[stats]{chisq.test}},
\code{\link[stats]{aov}},\code{\link[stats]{kruskal.test}},\code{\link[stats]{lm}}
}
