% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{system.ls}
\alias{system.ls}
\title{Run ls command}
\usage{
system.ls(pattern = "", directory = "", error = FALSE)
}
\arguments{
\item{pattern}{pattern to match files}

\item{directory}{base directory command should be run from}

\item{error}{logical indicating whether to throw an error if no matching founds found. Defaults to False.}
}
\value{
paths returned by ls command
}
\description{
Runs ls command on system. This is a workaround since list.files can not match patterns based on subdirectory structure.
}
