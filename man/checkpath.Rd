% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/YFinder.R
\name{checkpath}
\alias{checkpath}
\title{Check if subclade is +ve and return path}
\usage{
checkpath(
  hg,
  jsonnest,
  calltable,
  callquality,
  count = 1,
  depth = 2,
  path = NULL,
  yfullsnp,
  numsnptotal
)
}
\arguments{
\item{hg}{The subclade eg R-Z2124}

\item{jsonnest}{the loaded JSON yfull file (use jsonlite::fromJSON() in package 
'jsonlite')}

\item{calltable}{object provided by calling function}

\item{callquality}{value provided by calling function, default 0.9}

\item{count}{value provided by calling function, count of ancestors who 
are +ve as well}

\item{depth}{value provided by calling function, default 2}

\item{path}{list provided by calling function}

\item{yfullsnp}{table provided by calling function}

\item{numsnptotal}{value provided by calling function. default 2. 
Min number of SNPs needed to give +ve call. 1 gives false calls.}
}
\value{
List object with Path and derived status TRUE or FALSE
}
\description{
Check if subclade is +ve and return path
}
