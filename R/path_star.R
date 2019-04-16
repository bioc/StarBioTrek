#' @title Get human KEGG pathway data and creates a network data.
#' @description pathnet creates a list of network data for each human pathway. The network data will be generated when interacting genes belong to that pathway.  
#' @param data  a list of network data as provided by getNETdata
#' @param genes.by.pathway  a list of pathway data as provided by ConvertedIDgenes
#' @importFrom igraph graph.data.frame induced.subgraph get.data.frame
#' @export
#' @return a list of network data for each pathway (interacting genes belong to that pathway)
#' @examples
#' lista_net<-pathnet(genes.by.pathway=pathway[1:5],data=netw)
pathnet<-function(genes.by.pathway,data){
  geneSymb_net_shar_pro<-data
  ds_shar_pro<-do.call("rbind", geneSymb_net_shar_pro)
  data_shar_pro<-as.data.frame(ds_shar_pro[!duplicated(ds_shar_pro), ]) 
  sdc_shar_pro<-unlist(data_shar_pro$gene_symbolA,data_shar_pro$gene_symbolB)
  m_shar_pro<-c(data_shar_pro$gene_symbolA)
  m2_shar_pro<-c(data_shar_pro$gene_symbolB)
  ss_shar_pro<-cbind(m_shar_pro,m2_shar_pro)
  data_pr_shar_pro<-as.data.frame(ss_shar_pro[!duplicated(ss_shar_pro), ]) 
  pathwayx<-  Uniform(genes.by.pathway)
  data<-data_pr_shar_pro

  lista_int<-list()
  for (k in 1:ncol(pathwayx)){
    print(colnames(pathwayx)[k])
    currentPathway_genes<-pathwayx[,k]
    colnames(data) <- c("gene_symbolA", "gene_symbolB")
    i <- sapply(data, is.factor)
    data[i] <- lapply(data[i], as.character)
    ver<-unlist(data)
    n<-unique(ver)
    s<-intersect(n,currentPathway_genes)
    g <- graph.data.frame(data,directed=FALSE)
    g2 <- induced.subgraph(graph=g,vids=s)
    aaa<-get.data.frame(g2)
    colnames(aaa)[1] <- 'V1'
    colnames(aaa)[2] <- 'V2'
    lista_int[[k]]<-aaa
    names(lista_int)[k]<-colnames(pathwayx)[k] 
  }
  return(lista_int)
}




#' @title Get human KEGG pathway data and the output of list_path_net  define the common genes.
#' @description listpathnet creates a list of interacting genes for each human pathway.   
#' @param lista_net  output of path_net
#' @param pathway_exp  pathway data as provided by getKEGGdata
#' @export
#' @return a list of genes for each pathway (interacting genes belong to that pathway)
#' @examples
#' lista_network<-pathnet(genes.by.pathway=pathway[1:5],data=netw)
#' list_path<-listpathnet(lista_net=lista_network,pathway=pathway[1:5])
listpathnet<-function(lista_net,pathway_exp){
  top3<-  Uniform(pathway_exp)
  pathway<-top3
v=list()
bn=list()
for (j in 1:length(lista_net)){
  cf<-lista_net[[j]]
  i <- sapply(cf, is.factor) 
  cf[i] <- lapply(cf[i], as.character)
  colnames(cf) <- c("m_shar_pro", "m2_shar_pro")
  m<-c(cf$m_shar_pro)
  m2<-c(cf$m2_shar_pro)
  s<-c(m,m2)
  fr<- unique(s)
  n<-as.data.frame(fr)
  if(length(n)==0){
    v[[j]]<-NULL
  }
  i <- sapply(n, is.factor) 
  n[i] <- lapply(n[i], as.character)
  #for (k in  1:ncol(pathway)){
  if (length(intersect(n$fr,pathway[,j]))==nrow(n)){
    print(paste("List of genes interacting in the same pathway:",colnames(pathway)[j]))
    aa<-intersect(n$fr,pathway[,j])
    v[[j]]<-aa
    names(v)[j]<-colnames(pathway)[j]
  }
}
return(v)}




#' @title Get human KEGG pathway data and a gene expression matrix in order to obtain a list with the gene expression for only pathways given in input .
#' @description GE_matrix creates a list of gene expression for pathways given by the user.   
#' @param DataMatrix  gene expression matrix (eg.TCGA data)
#' @param genes.by.pathway a list of pathway data as provided by GetData and ConvertedID_genes
#' @export
#' @return a list for each pathway ( gene expression level belong to that pathway)
#' @examples
#' list_path_gene<-GE_matrix(DataMatrix=tumo[,1:2],genes.by.pathway=pathway[1:5])
GE_matrix<-function(DataMatrix,genes.by.pathway) {
  pathwayfrom<-genes.by.pathway
  top3<-  Uniform(pathwayfrom)
  pathway<-top3
zz<-as.data.frame(DataMatrix)
v<-list()
for ( k in 1: ncol(pathway)){
  #k=2
  if (length(intersect(rownames(zz),pathway[,k])!=0)){
    print(colnames(pathway)[k])
  currentPathway_genes_list_common <- intersect(rownames(zz), currentPathway_genes<-pathway[,k])
  currentPathway_genes_list_commonMatrix <- as.data.frame(zz[currentPathway_genes_list_common,])
  rownames(currentPathway_genes_list_commonMatrix)<-currentPathway_genes_list_common
  v[[k]]<- currentPathway_genes_list_commonMatrix
  names(v)[k]<-colnames(pathway)[k]
  }
}  
return(v)
}



#' @title Get human KEGG pathway data and a gene expression matrix in order to obtain a matrix with the mean gene expression for only pathways given in input .
#' @description GE_matrix creates a matrix of mean gene expression levels for pathways given by the user.   
#' @param DataMatrix  gene expression matrix (eg.TCGA data)
#' @param genes.by.pathway  list of pathway data as provided by getKEGGdata
#' @export
#' @return a matrix for each pathway (mean gene expression level belong to that pathway)
#' @examples
#' list_path_plot<-GE_matrix_mean(DataMatrix=tumo[,1:2],genes.by.pathway=pathway[1:5])
GE_matrix_mean<-function(DataMatrix,genes.by.pathway) {
  pathwayfrom<-genes.by.pathway
  top3<-  Uniform(pathwayfrom)
  pathway<-top3
  path_name<-sub(' ', '_',colnames(pathway))
  d_pr<- gsub(" - Homo sapiens (human)", "", path_name, fixed="TRUE")
  colnames(pathway)<-d_pr
  zz<-as.data.frame(rowMeans(DataMatrix))
  v<-list()
  for ( k in 1: ncol(pathway)){
    #k=2
    if (length(intersect(rownames(zz),pathway[,k])!=0)){
      print(colnames(pathway)[k])
      currentPathway_genes_list_common <- intersect(rownames(zz), currentPathway_genes<-pathway[,k])
      currentPathway_genes_list_commonMatrix <- as.data.frame(zz[currentPathway_genes_list_common,])
      rownames(currentPathway_genes_list_commonMatrix)<-currentPathway_genes_list_common
      v[[k]]<- currentPathway_genes_list_common
      names(v)[k]<-colnames(pathway)[k]
    }
  }  
  PEAmatrix <- matrix( 0,nrow(DataMatrix),ncol(pathway))
  rownames(PEAmatrix) <- as.factor(rownames(DataMatrix))
  colnames(PEAmatrix) <-  as.factor(colnames(pathway))
  for (i in 1:length(v)){
  PEAmatrix[v[[i]],i]<-zz[v[[i]],]
  }
  PEAmatrix<-PEAmatrix[which(rowSums(PEAmatrix) > 0),]
  return(PEAmatrix)
}







#' @title For TCGA data get human pathway data and creates a matrix with the average of genes for each pathway.
#' @description average creates a matrix with a summarized value for each pathway  
#' @param pathwayexpsubset list of pathway data
#' @export
#' @return a matrix value for each pathway 
#' @examples
#' list_path_gene<-GE_matrix(DataMatrix=Data_CANCER_normUQ_fil,genes.by.pathway=pathway[1:50])
#' score_mean<-average(pathwayexpsubset=list_path_gene)
average<-function(pathwayexpsubset){
#pathwayexpsubset_fil<-pathwayexpsubset[-which(lapply(pathwayexpsubset,is.null) == T)]
  pathwayexpsubset_fil<-pathwayexpsubset[!sapply(pathwayexpsubset, is.null)]
  
  gslist<-lapply(pathwayexpsubset_fil,function(x) colMeans(x))
gsmatrix<-do.call("rbind",gslist)
return(gsmatrix)
}

#my_hclust_gene <- hclust(dist(score_mean), method = "complete")
#install.packages("dendextend")

# load package
#library(dendextend)

#as.dendrogram(my_hclust_gene) %>%
#  plot(horiz = TRUE)
#my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = 2) 
#my_gene_col <- data.frame(cluster = ifelse(test = my_gene_col == 1, yes = "cluster 1", no = "cluster 2"))


#my_sample_col <- data.frame(sample = rep(c("normal", "tumour"), c(113,113)))
#rownames(my_sample_col) <- colnames(score_mean)
#library(pheatmap)
#pheatmap(score_mean, annotation_row = my_gene_col, annotation_col = my_sample_col)

#' @title For TCGA data get human pathway data and creates a measure of cross-talk among pathways 
#' @description eucdistcrtlk creates a matrix with euclidean distance for pairwise pathways  
#' @param dataFilt TCGA matrix
#' @param pathway_exp list of pathway data
#' @export
#' @return a matrix value for each pathway 
#' @examples
#' score_euc_dista_t<-eucdistcrtlk(dataFilt=tumo[,1:2],pathway_exp=pathway[1:5])
eucdistcrtlk <- function(dataFilt,pathway_exp){
  listpathgenees<-GE_matrix(dataFilt,pathway_exp)
  PEAmatrix<-average(listpathgenees)
  #step 5 distance
  # EUCLIDEA DISTANCE
  df=combn(rownames(PEAmatrix),2) # possibili relazioni tra i pathway
  df=t(df)
  ma_d<-matrix(0,nrow(df),ncol(PEAmatrix)) # creo matrix che conterr? le distanze
  colnames(ma_d)<-colnames(PEAmatrix) # colnames conterr? il nome dei pazienti
  for ( p in 1: ncol(PEAmatrix)){ # per ogni paziente
    print(paste0("Calculating cross-talk........... for patient n ", p, " of ", ncol(PEAmatrix) ))
    patients <- (PEAmatrix)[,p] 
    distance<-dist(patients) # calcolo distanza EUCLIDEA tra le possibile combinazioni
    ma_d[,p]<-distance
  }
  euc_dist<-cbind(df,ma_d) # inserisco label con le relazioni tra i pathway
  
  PathwaysPair <- paste( as.matrix(euc_dist[,1] ), as.matrix(euc_dist[,2] ),sep="_" )
  pp<-as.data.frame(euc_dist[,3:ncol(euc_dist)])
  rownames(pp)<-PathwaysPair
  for( i in 1: ncol(pp)){
    pp[,i] <- as.numeric(as.character(pp[,i]))
  }
  
  return(pp)
}




#heatmap(as.matrix(score_euc_dista_t),scale="column",margins =c(12,9))


#listpathgenees<-GE_matrix(dataFilt,pathway_exp)
#listpathgenees<-listpathgenees[!sapply(listpathgenees, is.null)]
#listpathgenees<-lapply(listpathgenees,function(x) apply(x,2,sd))


#listpathgenees<-listpathgenees[!sapply(listpathgenees, function(x) all(is.na(x)))]
#gsmatrix_sd<-do.call("rbind",listpathgenees)


#' @title For TCGA data get human pathway data and creates a measure of standard deviations among pathways 
#' @description stdv creates a matrix with standard deviation for pathways  
#' @param gslist pathway data
#' @export
#' @return a matrix value for each pathway 
#' @examples
#' list_path_gene<-GE_matrix(DataMatrix=tumo[,1:2],genes.by.pathway=pathway[1:5])
#' score_stdev<-stdv(gslist=list_path_gene)
stdv<-function(gslist){
  gslist<-gslist[!sapply(gslist, function(x) all(is.na(x)))]
  gslist<-lapply(gslist,function(x) apply(x,2,sd))
  gslist<-gslist[!sapply(gslist, function(x) all(is.na(x)))]
  gsmatrix_sd<-do.call("rbind",gslist)
  return(gsmatrix_sd)
}





#' @title For TCGA data get human pathway data and creates a measure of discriminating score among pathways 
#' @description dsscorecrtlk creates a matrix with  discriminating score for pathways  
#' @param dataFilt TCGA matrix
#' @param pathway_exp a list of pathway data
#' @export
#' @return a matrix value for each pathway 
#' @examples
#' cross_talk_st_dv<-dsscorecrtlk(dataFilt=tumo[,1:2],pathway_exp=pathway[1:5])
dsscorecrtlk<-function(dataFilt,pathway_exp){
  listpathgenees<-GE_matrix(dataFilt,pathway_exp)
  PEAmatrix<-average(listpathgenees)
  
  # standard deviation
  PEAmatrix_sd<-stdv(listpathgenees)
  r<-intersect(rownames(PEAmatrix),rownames(PEAmatrix_sd))
  PEAmatrix<-PEAmatrix[r,]
  PEAmatrix_sd<-PEAmatrix_sd[r,]
  
  #step 5 distance
  # EUCLIDEA DISTANCE
  df=combn(rownames(PEAmatrix),2) # possibili relazioni tra i pathway
  df=t(df)
  ma_d<-matrix(0,nrow(df),ncol(PEAmatrix)) # creo matrix che conterr? le distanze
  colnames(ma_d)<-colnames(PEAmatrix) # colnames conterr? il nome dei pazienti
  for ( p in 1: ncol(PEAmatrix)){ # per ogni paziente
    print(paste0("Calculating cross-talk........... for patient n ", p, " of ", ncol(PEAmatrix) ))
    patients <- (PEAmatrix)[,p] 
    distance<-dist(patients) # calcolo distanza EUCLIDEA tra le possibile combinazioni
    ma_d[,p]<-distance
  }
  
  df=combn(rownames(PEAmatrix_sd),2) 
  df=t(df)
  ma<-matrix(0,nrow(df),ncol(PEAmatrix_sd)) # creo matrix che conterr le somme delle dev st
  colnames(ma)<-colnames(PEAmatrix_sd) # colnames conterra il nome dei pazienti
  for ( p in 1: ncol(PEAmatrix_sd)){ # per ogni paziente
    patients <- (PEAmatrix_sd)[,p] 
    out <- apply(df, 1, function(x) sum(patients[x])) # calcolo somma delle dev standard tra le possibili combinazioni
    ma[,p]<-out
  }
  score<-ma_d/ma # discriminating score M1-M2/S1+S2
  score<- cbind(df,score)  
  
  PathwaysPair <- paste( as.matrix(score[,1] ), as.matrix(score[,2] ),sep="_" )
  pp2<-as.data.frame(score[,3:ncol(score)])
  rownames(pp2)<-PathwaysPair
  for( i in 1: ncol(pp2)){
    pp2[,i] <- as.numeric(as.character(pp2[,i]))
  

  }
  return(pp2)
}


#' @title SVM classification for each feature
#' @description svm class creates a list with AUC, Accuracy, Sensitivity, Specificity values 
#' @param TCGA_matrix gene expression matrix where the first two columns represent the interacting pathways.
#' @param nfs nfs split data into a training  and test set
#' @param tumour barcode samples for a class
#' @param normal barcode samples for another class
#' @param Target label for the classes
#' @export
#' @importFrom e1071 tune svm 
#' @importFrom ROCR prediction performance 
#' @importFrom  MLmetrics  Accuracy Sensitivity Specificity
#' @return a list with AUC value for pairwise pathway 
#' @examples
#' \dontrun{
#' nf <- 60
#' res_class<-svm_classification(TCGA_matrix=score_euc_dista[1:30,],nfs=nf,
#' normal=colnames(norm[,1:10]),tumour=colnames(tumo[,1:10]))}
svm_classification<-function(TCGA_matrix,tumour,normal,nfs){
  #library("e1071")
  #library(ROCR)
  #TCGA_matrix<-TCGA_matrix[complete.cases(TCGA_matrix), ]
  #scoreMatrix <- as.data.frame(TCGA_matrix[,3:ncol(TCGA_matrix)])
  #scoreMatrix <-as.data.frame(scoreMatrix)
  #for( i in 1: ncol(scoreMatrix)){
   # scoreMatrix[,i] <- as.numeric(as.character(scoreMatrix[,i]))
  #}

  
  #PathwaysPair <- paste( as.matrix(TCGA_matrix[,1] ), as.matrix(TCGA_matrix[,2] ),sep="_" )
  
  #rownames(scoreMatrix) <-PathwaysPair

  scoreMatrix<-TCGA_matrix
  tDataMatrix<-as.data.frame(t(scoreMatrix))
  #tDataMatrix$Target[,1]<-0
  
  tDataMatrix<-cbind(Target=0,tDataMatrix )

  tum<-intersect(rownames(tDataMatrix),tumour)
  nor<-intersect(rownames(tDataMatrix),normal)
  #tDataMatrix$
    
  Dataset_g1<-tDataMatrix[nor,]
  Dataset_g3<- tDataMatrix[tum,]
    
  
#training=read.table('C:/Users/UserInLab05/Desktop/trai.txt',header = TRUE)
#testset=read.table('C:/Users/UserInLab05/Desktop/test.txt',header = TRUE)

  Dataset_g1$Target <- 0
  Dataset_g3$Target<-1
#Dataset_g3 <- Dataset_g3[Dataset_g3$Target <- 1, ]
 set.seed(42) 
tab_g1_training <- sample(rownames(Dataset_g1),round(nrow(Dataset_g1) / 100 * nfs ))
tab_g3_training <- sample(rownames(Dataset_g3),round(nrow(Dataset_g3) / 100 * nfs ))
tab_g1_testing <- setdiff(rownames(Dataset_g1),tab_g1_training)
tab_g3_testing <- setdiff(rownames(Dataset_g3),tab_g3_training)

FR<-intersect(rownames(Dataset_g1),tab_g1_training)

#rownames(Dataset_g1)<-Dataset_g1[,1]
G1<-Dataset_g1[FR,]

FR1<-intersect(rownames(Dataset_g3),tab_g3_training)
#rownames(Dataset_g3)<-Dataset_g3$ID

G3<-Dataset_g3[FR1,]
training<-rbind(G1,G3)
training$Target<-as.factor(training$Target) #deve essere factor per avere confusion matrix



inter1<-intersect(rownames(Dataset_g1),tab_g1_testing)
#rownames(Dataset_g1)<-Dataset_g1$ID

G1_testing<-Dataset_g1[inter1,]

inter2<-intersect(rownames(Dataset_g3),tab_g3_testing)
#rownames(Dataset_g3)<-Dataset_g3$ID
G3_testing<-Dataset_g3[inter2,]

testing<-rbind(G1_testing,G3_testing)

x <- subset(training, select=-Target)
y <- as.factor(training$Target)
#testing[,2]<-NULL
z<-subset(testing, select=-Target)

zi<-testing$Target#factor for classification
#svm_tune <- tune(svm, train.x=x, train.y=y, 
                 #kernel="radial", ranges=list(cost=10^(-1:2), gamma=c(.5,1,2)))
PEAmatrix <- matrix( 0, ncol(training)-1,4)
rownames(PEAmatrix) <- colnames(x)
colnames(PEAmatrix) <-  c("auc","acc","sens","spec")
#prediction_labels=list()
#training2<-rbind(G1,G3)
for( k in 2: ncol(training)){
  print(paste0(colnames(training)[k]," pathway n. ", k-1, " of ",ncol(training)-1))
  

  svm_tune <- tune(svm, train.x=as.numeric(as.character(x[,k-1])), train.y=y, 
                   kernel="radial", ranges=list(cost=10^(-1:2), gamma=c(.5,1,2)))

  #print(svm_tune)
  
  svm_model_after_tune <- svm(Target ~ ., data=training[,c(1,k)], kernel="radial", cost=svm_tune$best.parameters$cost, gamma=svm_tune$best.parameters$gamma,probability = TRUE)

  j=k-1
  z2=z[,j]
  z3<-as.data.frame(z2)
  rownames(z3)<-rownames(z)
  #colnames(z3)<-as.character(paste("X",j,sep = ""))
  colnames(z3)<-colnames(z)[j]
  #pred3 <- predict(svm_model_after_tune, newdata = z3,decision.values = TRUE,probability = FALSE)
  
  
  
  #bb2<-attr(pred3,"decision.values")
  #bb<-attr(pred4,"probabilities")
  
  ### per avere matrice di confusione ma non per AUC ####################
  pred3 <- predict(svm_model_after_tune, newdata = z3,decision.values = TRUE,probability = FALSE)
  #bb2<-attr(pred3,"decision.values")
  #svm.roc <- prediction(bb2, zi)
  #svm.auc <- performance(svm.roc, 'tpr', 'fpr')
  #perf <- performance(svm.roc, "auc")
  #auc2<-perf@y.values[[1]]
  a<-table(zi,pred3)
  acc<-Accuracy(pred3, zi)
  sens<-Sensitivity(zi, pred3)
  spec<-Specificity(zi, pred3)
  ########################################################
  
  
  
  ########################################################
  pred4 <- predict(svm_model_after_tune, newdata = z3,probability = TRUE,decision.values = FALSE)
  bb<-attr(pred4,"probabilities")
  svm.roc <- prediction(bb[,2], zi)
  
  svm.auc <- performance(svm.roc, 'tpr', 'fpr')
  perf <- performance(svm.roc, "auc")
  auc3<-perf@y.values[[1]]
  
  
  
  
  #a<-table(zi,pred4)
  #library(MLmetrics)
  #auc<-AUC(pred3, zi)
  #acc<-Accuracy(pred4, zi)
  #sens<-Sensitivity(zi, pred4)
  #spec<-Specificity(zi, pred3)

  #bgr<-as.data.frame(pred3)

  # https://stats.stackexchange.com/questions/142255/svm-in-r-package-e1071-predicting-class-using-predict
  # r predict probability true and decision values
  
  PEAmatrix[k-1,1] <- auc3
  PEAmatrix[k-1,2] <- acc
  PEAmatrix[k-1,3] <- sens
  PEAmatrix[k-1,4] <- spec
}
  #prediction_labels[[k]]<-bgr
  #prediction_labels_performance<-list(PEAmatrix,prediction_labels)
return(PEAmatrix)
}




#' @title Multilayer analysis Cava et al. BMC Genomics 2017
#' @description IPPI function takes as input pathway and network data in order to select genes with central role in that pathway. Please see Cava et al. 2017 BMC Genomics 
#' @param pathax pathway matrix Please see example path for format
#' @param netwa a dataframe Please see example path for format netw
#' @export
#' @importFrom SpidermiR SpidermiRanalyze_degree_centrality 
#' @importFrom igraph graph.data.frame induced.subgraph get.data.frame
#' @return a list with driver genes for each pathway 
#' @examples 
#' \dontrun{
#' DRIVER_SP<-IPPI(pathax=pathway_matrix[,1:3],netwa=netw_IPPI[1:50000,])}
IPPI<-function(pathax,netwa){
  topwhol_net<-SpidermiRanalyze_degree_centrality(netwa)
  colnames(netwa)<-c("m_shar_pro","m2_shar_pro")
  m<-c(netwa$m_shar_pro)
  m2<-c(netwa$m2_shar_pro)
  s<-c(m,m2)
  fr<- unique(s)
  lista_nett<-IPPIpath_net(pathway=pathax,data=netwa)
  lista_netw<- delete.NULLs(lista_nett)
  a<-intersect(names(lista_netw),colnames(pathax))
  pat<-pathax[,a]
  list_path<-IPPIlist_path_net(lista_net=lista_netw,pathway=pat)
  topwhol=list()
  df=list()
  for (i in 1:length(lista_netw)){
    #  for (i in 1:4){
    #i=1
    print(paste(i,"INTERACTION",names(lista_netw)[i]))
    #if(is.data.frame(lista_net[[i]])){
    #if(nrow(lista_net[[i]])){
    topwhol[[i]]<-SpidermiRanalyze_degree_centrality(lista_netw[[i]])
    topwhol[[i]]<-topwhol[[i]][which(topwhol[[i]]$Freq!=0),]
    topwhol[[i]]$topwhol_net<-""
    topwhol[[i]]$wi<-""
    topwhol[[i]]$d_expected_path<-""
    topwhol[[i]]$driver<-""
    p=list()
    for (j in 1:nrow(topwhol_net)){
      #j=644
      if (length(intersect(topwhol_net[j,1],topwhol[[i]]$dfer))!=0)  {
        d_expected_path<-topwhol_net[j,2]*(length(list_path[[i]]))/(length(fr))
        
        index<- which(as.character(topwhol_net[j,1])==as.character(topwhol[[i]]$dfer))
        sa<-topwhol[[i]][which(as.character(topwhol_net[j,1])==as.character(topwhol[[i]]$dfer)),]
        wi<-sa$Freq/length(list_path[[i]])
        #wi=sa$Freq-d_expected_path
        topwhol[[i]][index,3]<-topwhol_net[j,2]
        topwhol[[i]][index,4]<-wi
        topwhol[[i]][index,5]<-d_expected_path
        if(d_expected_path<wi){
          topwhol[[i]][index,6]<-"driver"
        }
      }
      if (length(intersect(topwhol_net[j,1],topwhol[[i]]$dfer))==0)  {
        df[[i]]<-topwhol_net[j,1]
      } 
    }
  }   
  for (k in 1:length(lista_netw)){
    names(topwhol)[k]<-names(lista_netw)[k]
  }
  return(topwhol)
}







#' @title Preparation for plotting cross-talk
#' @description plot_cross_talk function takes as input pathway data and prepares the data to visualize (e.g. ggplot2, qqgraph, igraph)
#' @param pathway_plot pathway 
#' @param gs_expre a gene expression matrix
#' @export
#' @return a list with correlation matrix and gene set for each gene
#' @examples 
#' formatplot<-plotcrosstalk(pathway_plot=pathway[1:6],gs_expre=tumo)
plotcrosstalk<-function(pathway_plot,gs_expre){
informationplot<-list()
  vv<-list()
for (i in 1:length(pathway_plot)){
iim<-intersect(pathway_plot[[i]],rownames(gs_expre))
vv[[i]]<-as.factor(iim)
}
names(vv)<-names(pathway_plot)
#select genes presents only in one pathway
vvs<-lapply(1:length(vv), function(n) as.factor(setdiff(vv[[n]], unlist(vv[-n]))))
names(vvs)<-names(pathway_plot)
ff<-unique(unlist(vvs))
tt<-gs_expre[rownames(gs_expre) %in% ff,]
# Find the number of elements
num.el <- sapply(vvs, length)
# Generate the matrix
res <- cbind(as.character(unlist(vvs)), rep(names(vvs), num.el))
rownames(res)<-res[,1]
order_matrix<-merge(res, tt, by = "row.names", all = TRUE)
#correlation
rownames(order_matrix)<-order_matrix[,1]
corr_matrix<-order_matrix[,4:ncol(order_matrix)]
bg<-cor(t(corr_matrix))

informationplot[[1]]<-bg
informationplot[[2]]<-order_matrix[,3]
names(informationplot)<-c("correlation","groups")
#qgraph(formatplot[[1]], minimum = 0.25, cut = 0.6, vsize = 5, groups = formatplot[[2]], 
 #      legend = TRUE, borders = FALSE,layoutScale=c(0.8,0.8))

#qgraph(formatplot[[1]],groups=formatplot[[2]], layout="spring", diag = FALSE, cut = 0.6,legend.cex = 0.5,vsize = 6,layoutScale=c(0.8,0.8))

return(informationplot)
}




#' @title Preparation for circle plot
#' @description circleplot function takes as input data derived by the function plotcrosstalk and plOt a circle plot. 
#' @param preplot a list as obtained from the function plotcrosstalk
#' @param scoregene a score for each gene with values included between -10 e +10
#' @export
#' @importFrom reshape2 dcast 
#' @return a list with correlation matrix and gene set for each gene
#' @examples 
#' formatplot<-plotcrosstalk(pathway_plot=pathway[1:6],gs_expre=tumo)
#' score<-runif(length(formatplot[[2]]), min=-10, max=+10)
#' circleplot(preplot=formatplot,scoregene=score)
circleplot<-function(preplot,scoregene){
mme<-as.data.frame(preplot$groups)
  mme[,2]<-rownames(mme)
prova<-dcast(mme[,2]~mme[,1], data=mme,fun.aggregate = length)
rownames(prova)<-prova[,1]
prova[,1]<-NULL
asd<-as.matrix(prova)
#score<-runif(99, min=-4, max=+4)
asda<-as.matrix(cbind(asd,scoregene))
colnames(asda)[ncol(asda)]<-c("logFC")
GOChord(asda, limit = c(0, 5))
}