% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subsetRows.R
\name{subsetRows}
\alias{subsetRows}
\alias{topvar}
\alias{topvar2}
\alias{topvar3}
\alias{subset_nonzero}
\alias{subset_nonzero2}
\alias{subset_rowname}
\title{Subset rows from a dataframe by various criteria}
\usage{
topvar(df, p, v = cv)

topvar2(df, thres, v = cv)

topvar3(df, n, v = cv)

subset_nonzero(df, p)

subset_nonzero2(df, n)

subset_rowname(x, y, order = T)
}
\arguments{
\item{df}{numeric dataframe}

\item{p}{percentage, 0-100}

\item{v}{method to calculate the variance (\code{\link{cv}}, \code{\link[stats]{mad}}, \code{\link[stats]{mean}},etc)}

\item{thres}{a threshold of variance for the cutoff}

\item{n}{integer}
}
\value{
dataframe
}
\description{
A group of functions to subset rows from a dataframe 
\itemize{
  \item topvar: subset top p% rows by variance (or other methods like mean/median) 
  \item topvar2: subset rows with at least thres variance.
  \item topvar3: subset top n rows by variance.
  \item \emph{NOTE: for topvars, rows with zero sum are automatically removed}
  \item subset_nonzero: subset rows with at least p% non-zero numbers. (For filtering non-expressed genes)
  \item subset_nonzero2: subset rows with at least n non-zero numbers
  \item subset_rowname: Subset one dataframe according to another, by rownames. Return df has same order as the reference df
}
}
