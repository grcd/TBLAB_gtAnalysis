\name{rungt}
\alias{rungt}
\title{Run the global test on the list of mirs}
\usage{
  rungt(mirs, X, Y, path, dbName, tables, numOverlapping,
    top = -1, ...)
}
\arguments{
  \item{mirs}{list of mirs}

  \item{X}{mRNA expression}

  \item{Y}{miRNA expression}

  \item{path}{path to database}

  \item{dbName}{database name}

  \item{tables}{prediction databases}

  \item{numOverlapping}{number of at least overlapping
  targets between databases}

  \item{top}{number of significant targets returned -1 is
  all}
}
\value{
  list of microRNAs and targets
}
\description{
  rungt is a wrapper for the gt() function
}
\details{
  Details follow.
}
\author{
  Maarten van Iterson, Sander Bervoets
}

