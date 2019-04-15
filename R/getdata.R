#' @title Get general information inside pathways.
#' @description GetData creates a list with genes inside the pathways. 
#' @param species  variable. The user can select the species of interest from SELECT_path_species(path_spec) 
#' @param pathwaydb variable. The user can select the pathway database of interest from SELECT_path_graphite(path_spec) 
#' @export
#' @importFrom graphite pathways pathwayTitle
#' @return a list of pathways
#' @examples
#' \dontrun{
#' species="hsapiens"
#' pathwaydb="pharmgkb"
#' path<-GetData(species,pathwaydb)}
GetData<-function(species,pathwaydb){
  humanpath <- pathways(species, pathwaydb)
  humanReactome<-humanpath
  le<-list()
  for (j in  1:length(humanReactome)){
    e<-humanReactome[[j]]
    print(paste0("Querying.............  ",pathwayTitle(e),"   ", j, " of ",length(humanReactome)," pathways"))
    le[[j]]<-e
  }
  names(le)<- names(humanReactome)
  return(le)
}


#' @title Get  genes inside pathways.
#' @description GetPathData creates a list of genes inside the pathways. 
#' @param path_ALL  variable. The user can select the variable as obtained by  GetData function
#' @export
#' @importFrom graphite nodes pathwayTitle
#' @return a list of pathways
#' @examples
#' pathway_ALL_GENE<-GetPathData(path_ALL=path[1:3])
GetPathData<-function(path_ALL){
  le<-list()
  for (j in  1:length(path_ALL)){
    e<-path_ALL[[j]]
    genes<-nodes(e,which = "proteins")
    print(paste0("Downloading.............  ",pathwayTitle(e),"   ", j, " of ",length(path_ALL)," pathways"))
    le[[j]]<-genes
  }
  names(le)<- names(path_ALL)
  return(le)
}




#' @title Get  interacting genes inside pathways.
#' @description GetPathNet creates a list of genes inside the pathways. 
#' @param path_ALL  variable. The user can select the variable as obtained by  GetData function
#' @export
#' @importFrom graphite edges pathwayTitle
#' @return a list of pathways
#' @examples
#' pathway_net<-GetPathNet(path_ALL=path[1:3])
GetPathNet<-function(path_ALL){
  le<-list()
  for (j in  1:length(path_ALL)){
    e<-path_ALL[[j]]
    genes<-edges(e,which = "proteins")
    print(paste0("Downloading.............  ",pathwayTitle(e),"   ", j, " of ",length(path_ALL)," pathways"))
    le[[j]]<-genes
  }
  names(le)<- names(path_ALL)
  return(le)
}



#' @title Get  interacting genes inside pathways.
#' @description GetPathNet creates a list of genes inside the pathways. 
#' @param path_ALL  variable. The user can select the variable as obtained by  GetData function
#' @export
#' @importFrom graphite nodes pathwayTitle convertIdentifiers
#' @return a list of pathways
#' @examples
#' pathway<-ConvertedIDgenes(path_ALL=path[1:3])
ConvertedIDgenes<-function(path_ALL){
  le<-list()
  for (j in  1:length(path_ALL)){
    e<-path_ALL[[j]]
    s1<-convertIdentifiers(e, "symbol")
    genes<-nodes(s1,which = "proteins")
    er <- sapply(strsplit(genes, split=':', fixed=TRUE), function(x) (x[2]))
    print(paste0("Mapping Uniprot ID to Gene Symbol, using convertIdentifiers of graphite package..........  ",pathwayTitle(e),"   ", j, " of ",length(path_ALL)," pathways"))
    #attr(mm, "names")<-NULL
    le[[j]]<-er
  }
  names(le)<- names(path_ALL)
  return(le)
}
  
  
 

#' @title Get network data from GeneMania.
#' @description getNETdata creates a data frame with network data. 
#' Network category can be filtered among: physical interactions, co-localization, genetic interactions and shared protein domain.
#' @param network  variable. The user can use the following parameters 
#' based on the network types to be used. PHint for Physical_interactions,
#' COloc for Co-localization, GENint for Genetic_interactions and
#' SHpd for Shared_protein_domains
#' @param organismID organism==NULL default value is homo sapiens.
#' @export
#' @importFrom SpidermiR  SpidermiRquery_spec_networks SpidermiRdownload_net SpidermiRprepare_NET
#' @return list with gene-gene (or protein-protein interactions)
#' @examples
#' \dontrun{
#' organismID="Saccharomyces_cerevisiae"
#' netw<-getNETdata(network="SHpd",organismID)}
getNETdata<-function(network,organismID=NULL){
  if( is.null(organismID) ){
    prr<-SpidermiRprepare_NET(organismID = 'Homo_sapiens',
                         data = SpidermiRdownload_net(data = SpidermiRquery_spec_networks(organismID = 'Homo_sapiens',network 
                                                                                          )))
  }
  if( !is.null(organismID) ){
    prr<-SpidermiRprepare_NET(organismID,
                              data = SpidermiRdownload_net(data = SpidermiRquery_spec_networks(organismID , 
                                                                                               network)))
  }
  return(prr)
}








