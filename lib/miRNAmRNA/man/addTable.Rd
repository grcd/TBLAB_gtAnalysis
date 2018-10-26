\name{addTable}
\alias{addTable}
\title{add table to the database}
\usage{
  addTable(file,
    tableName = c("pita", "microcosm", "targetscan"),
    parser = function(file, Org) stop("Parser not defined!"),
    dbName, path, Org = c("Mm", "Hs"))
}
\arguments{
  \item{file}{file and path}

  \item{tableName}{either "pita", "microcosm",
  "targetscan"}

  \item{dbName}{database name}

  \item{path}{directory where database should be stored}

  \item{parser}{user-defined parser function}

  \item{Org}{organism use abbreviation e.g. Mm and Hs}
}
\value{
  NULL
}
\description{
  Add table to the database
}
\details{
  PITA, microcosm and targetscan data for Mm were
  downloaded from
  \url{http://genie.weizmann.ac.il/pubs/mir07/catalogs/PITA_targets_mm9_0_0_TOP.tab.gz}
  \url{ftp://ftp.ebi.ac.uk/pub/databases/microcosm/v5/arch.v5.gff.mus_musculus.zip}
  \url{http://www.targetscan.org//mmu_60//mmu_60_data_download/Conserved_Site_Context_Scores.txt.zip}
}
\author{
  Maarten
}

