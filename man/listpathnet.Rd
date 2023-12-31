% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/path_star.R
\name{listpathnet}
\alias{listpathnet}
\title{Get human KEGG pathway data and the output of list_path_net  define the common genes.}
\usage{
listpathnet(lista_net, pathway_exp)
}
\arguments{
\item{lista_net}{output of path_net}

\item{pathway_exp}{pathway data as provided by getKEGGdata}
}
\value{
a list of genes for each pathway (interacting genes belong to that pathway)
}
\description{
listpathnet creates a list of interacting genes for each human pathway.
}
\examples{
lista_network<-pathnet(genes.by.pathway=pathway[1:5],data=netw)
list_path<-listpathnet(lista_net=lista_network,pathway=pathway[1:5])
}
