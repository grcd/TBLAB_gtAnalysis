\name{toTable}
\alias{toTable}
\title{toTable}
\usage{
  toTable(results, method = "BH", alpha = 0.05,
    level = c("global", "local"),
    org = c("Hs", "Mm", "None"))
}
\arguments{
  \item{results}{gtrun results}

  \item{method}{see p.adjust methods}

  \item{alpha}{singificance level}

  \item{level}{local or global}

  \item{org}{either Hs or Mm are implemented yet}
}
\value{
  data.frame
}
\description{
  toTable generates table from gtrun result list
}
\details{
  Details follow
}
\author{
  Maarten van Iterson
}

