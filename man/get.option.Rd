% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.option.R
\name{get.option}
\alias{get.option}
\title{Helper function to recursively get an VariTAS option}
\usage{
get.option(name, varitas.options = NULL, nesting.character = "\\\\.")
}
\arguments{
\item{name}{Option name}

\item{varitas.options}{Optional list of options to search in}

\item{nesting.character}{String giving Regex pattern of nesting indication string. Defaults to '\\.'}
}
\value{
value Requested option
}
\description{
Helper function to recursively get an VariTAS option
}
