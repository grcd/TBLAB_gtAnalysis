\name{overlap}
\alias{overlap}
\title{overlapping targets}
\usage{
  overlap(path, dbName, tables = NULL, mir,
    numOverlapping = 2, full = TRUE)
}
\arguments{
  \item{path}{path to database}

  \item{dbName}{database name}

  \item{tables}{tables in the database}

  \item{mir}{microRNA identifier}

  \item{numOverlapping}{number of at least overlapping
  targets}

  \item{full}{return full overlapping table or reduced only
  overlap}
}
\value{
  matrix with overlapping targets
}
\description{
  Determine overlapping targets between prediction tools
}
\details{
  Details follows
}
\author{
  Maarten van Iterson
}

