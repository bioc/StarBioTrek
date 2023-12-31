% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/path_star.R
\name{IPPI}
\alias{IPPI}
\title{Multilayer analysis Cava et al. BMC Genomics 2017}
\usage{
IPPI(pathax, netwa)
}
\arguments{
\item{pathax}{pathway matrix Please see example path for format}

\item{netwa}{a dataframe Please see example path for format netw}
}
\value{
a list with driver genes for each pathway
}
\description{
IPPI function takes as input pathway and network data in order to select genes with central role in that pathway. Please see Cava et al. 2017 BMC Genomics
}
\examples{
\dontrun{
DRIVER_SP<-IPPI(pathax=pathway_matrix[,1:3],netwa=netw_IPPI[1:50000,])}
}
