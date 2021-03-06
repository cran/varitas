% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/overwrite.varitas.options.R
\name{overwrite.varitas.options}
\alias{overwrite.varitas.options}
\title{overwrite.varitas.options}
\usage{
overwrite.varitas.options(config.file)
}
\arguments{
\item{config.file}{Path to config file that should be used to overwrite options}
}
\value{
None
}
\description{
Overwrite VariTAS options with options provided in config file.
}
\examples{
\dontrun{
config <- file.path(path.package('varitas'), 'config.yaml')
overwrite.varitas.options(config)
}

}
