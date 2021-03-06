% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{removeNonExpressed}
\alias{removeNonExpressed}
\title{Remove counts below threshold}
\usage{
removeNonExpressed(counts, countCutoff = 5)
}
\arguments{
\item{counts}{Counts from a nano object}

\item{countCutoff}{The expression threshold}
}
\description{
Removes the annotation columns from a nano count matrix
}
\examples{
removeNonExpressed()
}
\keyword{utility}

