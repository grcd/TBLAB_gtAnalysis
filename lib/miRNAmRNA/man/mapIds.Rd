\name{mapIds}
\alias{mapIds}
\title{Map identifiers}
\usage{
  mapIds(x, from, to = "ENTREZ", organism = c("Hs", "Mm"))
}
\arguments{
  \item{x}{vector of ids type 'from'}

  \item{from}{name of the current identifier}

  \item{to}{name of identifier mapped to}

  \item{organism}{organism}
}
\value{
  vector of ids type 'to'
}
\description{
  Mapping between various gene identifiers using
  .db-package.
}
\details{
  Details follow.
}
\author{
  Maarten van Iterson
}

