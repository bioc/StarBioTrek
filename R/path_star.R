#' @title Get human KEGG pathway data and network data in order to define the common gene.
#' @description path_net creates a list of network data for each human pathway. The network data will be generated when interacting genes belong to that pathway.  
#' @param net_type  network data as provided by getNETdata
#' @param pathway  pathway data as provided by getKEGGdata
#' @export
#' @return a list of network data for each pathway (interacting genes belong to that pathway)
#' @examples
#' lista_net<-path_net(pathway=path,net_type=netw)
path_net<-function(pathway,net_type){
  lista_int<-list()
  for (k in 1:ncol(pathway)){
    #k=1 
    print(paste(k,"PATHWAY",colnames(pathway)[k]))
    currentPathway_genes<-pathway[,k]
    common1 <- intersect( net_type$m_shar_pro, currentPathway_genes)
    common2 <- intersect( net_type$m2_shar_pro, currentPathway_genes)
    if (length(common1)==0 & length(common2)==0 ){
      mago2<-character(length = 0)
    }
    if (length(common1)!=0 | length(common2)!=0 ){
      b=list()
      for (i in 1:length(common1)){
        x<-common1[i]
        n<-overlap(net_type,x,currentPathway_genes)
        b[[i]]<-n
      }
      v<-do.call("rbind", b)
      c=list()
      for (i in 1:length(common2)){
        x<-common1[i]
        n<-overlap(net_type,x,currentPathway_genes)
        c[[i]]<-n
      }
      v2<-do.call("rbind", b)
      mago<-rbind(v,v2)
      mago2<-mago[!duplicated(mago), ]
    }
    
    if (length(mago2)!=0){
      lista_int[[k]]<-mago2
    }
    if (length(mago2)==0){
      lista_int[[k]]<-"0"} 
    
    names(lista_int)[k]<-colnames(pathway)[k] 
  }   
  return(lista_int)
}


#' @title Get human KEGG pathway data and output of path_net in order to define the common genes.
#' @description list_path_net creates a list of interacting genes for each human pathway.   
#' @param lista_net  output of path_net
#' @param pathway  pathway data as provided by getKEGGdata
#' @export
#' @return a list of genes for each pathway (interacting genes belong to that pathway)
#' @examples
#' lista_netw<-path_net(pathway=path,net_type=netw)
#' list_path<-list_path_net(lista_net=lista_netw,pathway=path)
list_path_net<-function(lista_net,pathway){
v=list()
bn=list()
for (j in 1:length(lista_net)){
  cf<-lista_net[[j]]
  i <- sapply(cf, is.factor) 
  cf[i] <- lapply(cf[i], as.character)
  m<-c(cf$m_shar_pro)
  m2<-c(cf$m2_shar_pro)
  s<-c(m,m2)
  fr<- unique(s)
  n<-as.data.frame(fr)
  if(length(n)==0){
    v[[j]]<-NULL
    
  }
  if(length(n)!=0){
  i <- sapply(n, is.factor) 
  n[i] <- lapply(n[i], as.character)
  #for (k in  1:ncol(pathway)){
  if (length(intersect(n$fr,pathway[,j]))==nrow(n)){
    print(paste("List of genes interacting in the same pathway:",colnames(pathway)[j]))
    aa<-intersect(n$fr,pathway[,j])
    v[[j]]<-aa
    names(v)[j]<-colnames(pathway)[j]
  }
}}
return(v)}




#' @title Get human KEGG pathway data and a gene expression matrix in order to obtain a matrix with the gene expression for only pathways given in input .
#' @description GE_matrix creates a matrix of gene expression for pathways given by the user.   
#' @param DataMatrix  gene expression matrix (eg.TCGA data)
#' @param pathway  pathway data as provided by getKEGGdata
#' @export
#' @return a matrix for each pathway ( gene expression level belong to that pathway)
#' @examples
#' list_path_gene<-GE_matrix(DataMatrix=tumo[,1:2],pathway=path)
GE_matrix<-function(DataMatrix,pathway) {
  path_name<-sub(' ', '_',colnames(pathway))
d_pr<- gsub(" - Homo sapiens (human)", "", path_name, fixed="TRUE")
colnames(pathway)<-d_pr
#zz<-as.data.frame(rowMeans(DataMatrix))
zz<-as.data.frame(DataMatrix)
v<-list()
for ( k in 1: ncol(pathway)){
  #k=2
  if (length(intersect(rownames(zz),pathway[,k])!=0)){
    print(colnames(path)[k])
  currentPathway_genes_list_common <- intersect(rownames(zz), currentPathway_genes<-pathway[,k])
  currentPathway_genes_list_commonMatrix <- as.data.frame(zz[currentPathway_genes_list_common,])
  rownames(currentPathway_genes_list_commonMatrix)<-currentPathway_genes_list_common
  v[[k]]<- currentPathway_genes_list_commonMatrix
  names(v)[k]<-colnames(pathway)[k]
  }
}  
#PEAmatrix <- matrix( 0,nrow(DataMatrix),ncol(pathway))
#rownames(PEAmatrix) <- as.factor(rownames(DataMatrix))
#colnames(PEAmatrix) <-  as.factor(colnames(pathway))
#for (i in 1:length(v)){
#PEAmatrix[v[[i]],i]<-zz[v[[i]],]
#}
#PEAmatrix<-PEAmatrix[which(rowSums(PEAmatrix) > 0),]
return(v)
}



#' @title Get human KEGG pathway data and a gene expression matrix in order to obtain a matrix with the mean gene expression for only pathways given in input .
#' @description GE_matrix creates a matrix of mean gene expression for pathways given by the user.   
#' @param DataMatrix  gene expression matrix (eg.TCGA data)
#' @param pathway  pathway data as provided by getKEGGdata
#' @export
#' @return a matrix for each pathway (mean gene expression level belong to that pathway)
#' @examples
#' list_path_plot<-matrix_plot(DataMatrix=tumo[,1:2],pathway=path)
matrix_plot<-function(DataMatrix,pathway) {
  path_name<-sub(' ', '_',colnames(pathway))
  d_pr<- gsub(" - Homo sapiens (human)", "", path_name, fixed="TRUE")
  colnames(pathway)<-d_pr
  zz<-as.data.frame(rowMeans(DataMatrix))
  v<-list()
  for ( k in 1: ncol(pathway)){
    #k=2
    if (length(intersect(rownames(zz),pathway[,k])!=0)){
      print(colnames(path)[k])
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













#' @title Get human KEGG pathway data and a gene expression matrix we obtain a matrix with the gene expression for only pathways given in input .
#' @description plotting_matrix creates a matrix of gene expression for pathways given by the user.   
#' @param DataMatrix  gene expression matrix (eg.TCGA data)
#' @param pathway  pathway data as provided by getKEGGdata
#' @param path_matrix  output of the function matrix_plot
#' @export
#' @return a plot for pathway cross talk
#' @examples
#' mt<-plotting_cross_talk(DataMatrix=tumo[,1:2],pathway=path,path_matrix=list_path_plot)
plotting_cross_talk<-function(DataMatrix,pathway,path_matrix){
  zz<-as.data.frame(rowMeans(DataMatrix))
  v<-list()
  for ( k in 1: ncol(pathway)){
    path_name<-sub(' ', '_',colnames(pathway))
    d_pr<- gsub(" - Homo sapiens (human)", "", path_name, fixed="TRUE")
    colnames(pathway)<-d_pr
    if (length(intersect(rownames(zz),pathway[,k])!=0)){
      print(colnames(path)[k])
      currentPathway_genes_list_common <- intersect(rownames(zz), currentPathway_genes<-pathway[,k])
      currentPathway_genes_list_commonMatrix <- as.data.frame(zz[currentPathway_genes_list_common,])
      rownames(currentPathway_genes_list_commonMatrix)<-currentPathway_genes_list_common
      v[[k]]<- as.factor(currentPathway_genes_list_common)
      names(v)[k]<-colnames(pathway)[k]
    }
  }
  vv<-list()
  mi<-t(path_matrix)
  
  dc<-cor(mi)
  for ( k in 1: length(v)){
    currentPathway_genes_list_common <- intersect(rownames(dc), v[[k]])
    a<-match(currentPathway_genes_list_common,rownames(dc))
    vv[[k]]<- a
    names(vv)[k]<-colnames(pathway)[k]
  }
  list_plt=list(corr=dc,gruppi=vv)
 #r<-qgraph(list_plt$corr, groups=list_plt$gruppi, mar=c(1,1,1,1),minimum=0.6)
  return(list_plt)
}




#' @title For TCGA data get human pathway data and creates a matrix with the average of genes for each pathway.
#' @description average creates a matrix with a summarized value for each pathway  
#' @param dataFilt TCGA matrix
#' @param pathway pathway data
#' @export
#' @return a matrix value for each pathway 
#' @examples
#' score_mean<-average(dataFilt=tumo[,1:2],path)
average<-function(dataFilt,pathway){
  DataMatrix<-dataFilt
  #dataFilt[ , "new.col"] <- gsub("\\|.*", "", rownames(dataFilt))
  #DataMatrix<-dataFilt[which(dataFilt$new.col!="?"),]
  #DataMatrix <- subset(DataMatrix, !duplicated(DataMatrix$new.col)) 
  #rownames(DataMatrix)<-DataMatrix$new.col
  #DataMatrix$new.col<-NULL

PEAmatrix <- matrix( 0, ncol(pathway),ncol(DataMatrix))
rownames(PEAmatrix) <- colnames(pathway)
colnames(PEAmatrix) <-  colnames(DataMatrix)
listIPA_pathways<-colnames(pathway)
for ( k in 1: nrow(PEAmatrix)){
  #k=1
  currentPathway <- colnames(pathway)[k]
  currentPathway_genes_list_common <- intersect(rownames(DataMatrix), currentPathway_genes<-pathway[,k])
  currentPathway_genes_list_commonMatrix <- DataMatrix[currentPathway_genes_list_common,]
  SumGenes <- colSums(currentPathway_genes_list_commonMatrix)
  AverageGenes <- SumGenes / length(currentPathway_genes_list_common)
  PEAmatrix[k,] <- AverageGenes
}
return(PEAmatrix)
}



  








#' @title For TCGA data get human pathway data and creates a measure of cross-talk among pathways 
#' @description euc_dist_crtlk creates a matrix with euclidean distance for pairwise pathways  
#' @param dataFilt TCGA matrix
#' @param pathway pathway data
#' @export
#' @return a matrix value for each pathway 
#' @examples
#' score_euc_dista<-euc_dist_crtlk(dataFilt=tumo[,1:2],path)
euc_dist_crtlk <- function(dataFilt,pathway){
  PEAmatrix<-average(dataFilt,pathway)
  #step 5 distance
  # EUCLIDEA DISTANCE
  df=combn(rownames(PEAmatrix),2) # possibili relazioni tra i pathway
  df=t(df)
  ma_d<-matrix(0,nrow(df),ncol(PEAmatrix)) # creo matrix che conterr? le distanze
  colnames(ma_d)<-colnames(PEAmatrix) # colnames conterr? il nome dei pazienti
  for ( p in 1: ncol(PEAmatrix)){ # per ogni paziente
    patients <- (PEAmatrix)[,p] 
    distance<-dist(patients) # calcolo distanza EUCLIDEA tra le possibile combinazioni
    ma_d[,p]<-distance
  }
  euc_dist<-cbind(df,ma_d) # inserisco label con le relazioni tra i pathway
  return(euc_dist)
}




#' @title For TCGA data get human pathway data and creates a measure of standard deviations among pathways 
#' @description st_dv creates a matrix with standard deviation for pathways  
#' @param DataMatrix TCGA matrix
#' @param pathway pathway data
#' @export
#' @return a matrix value for each pathway 
#' @examples
#' stand_dev<-st_dv(DataMatrix=tumo[,1:2],pathway=path)
st_dv<-function(DataMatrix,pathway){
#DataMatrix<-dataFilt

#dataFilt[ , "new.col"] <- gsub("\\|.*", "", rownames(dataFilt))
#DataMatrix<-dataFilt[which(dataFilt$new.col!="?"),]
#DataMatrix <- subset(DataMatrix, !duplicated(DataMatrix$new.col)) 
#rownames(DataMatrix)<-DataMatrix$new.col
#DataMatrix$new.col<-NULL

PEAmatrix_sd <- matrix( 0, ncol(pathway),ncol(DataMatrix))
rownames(PEAmatrix_sd) <- colnames(pathway)
colnames(PEAmatrix_sd) <-  colnames(DataMatrix)
for ( k in 1: nrow(PEAmatrix_sd)){
  print(colnames(pathway)[k])
  currentPathway <- colnames(pathway)[k]
  currentPathway_genes_list_common <- intersect( rownames(DataMatrix), currentPathway_genes<-pathway[,k])
  currentPathway_genes_list_commonMatrix <- DataMatrix[currentPathway_genes_list_common,]
  stdev<-apply(currentPathway_genes_list_commonMatrix,2,sd) #deviazione standard dei pathway
  PEAmatrix_sd[k,] <- stdev
  }
return(PEAmatrix_sd)
}






#' @title For TCGA data get human pathway data and creates a measure of discriminating score among pathways 
#' @description ds_score_crtlk creates a matrix with  discriminating score for pathways  
#' @param dataFilt TCGA matrix
#' @param pathway pathway data
#' @export
#' @return a matrix value for each pathway 
#' @examples
#' cross_talk_st_dv<-ds_score_crtlk(dataFilt=tumo[,1:2],pathway=path)
ds_score_crtlk<-function(dataFilt,pathway){
  PEAmatrix<-average(dataFilt,pathway)
  #step 5 distance
  # EUCLIDEA DISTANCE
  df=combn(rownames(PEAmatrix),2) # possibili relazioni tra i pathway
  df=t(df)
  ma_d<-matrix(0,nrow(df),ncol(PEAmatrix)) # creo matrix che conterr? le distanze
  colnames(ma_d)<-colnames(PEAmatrix) # colnames conterr? il nome dei pazienti
  for ( p in 1: ncol(PEAmatrix)){ # per ogni paziente
    patients <- (PEAmatrix)[,p] 
    distance<-dist(patients) # calcolo distanza EUCLIDEA tra le possibile combinazioni
    ma_d[,p]<-distance
  }
  PEAmatrix_sd<-st_dv(dataFilt,pathway)
  df=combn(rownames(PEAmatrix_sd),2) 
  df=t(df)
  ma<-matrix(0,nrow(df),ncol(PEAmatrix_sd)) # creo matrix che conterr? le somme delle dev st
  colnames(ma)<-colnames(PEAmatrix_sd) # colnames conterr? il nome dei pazienti
  for ( p in 1: ncol(PEAmatrix_sd)){ # per ogni paziente
    patients <- (PEAmatrix_sd)[,p] 
    out <- apply(df, 1, function(x) sum(patients[x])) # calcolo somma delle dev standard tra le possibili combinazioni
    ma[,p]<-out
  }
  score<-ma_d/ma # discriminating score M1-M2/S1+S2
  score<- cbind(df,score)  
return(score)
}



#' @title SVM classification for each feature
#' @description svm class creates a list with auc value  
#' @param TCGA_matrix gene expression matrix
#' @param nfs nfs split data into a training  and test set
#' @param tumour barcode samples for a class
#' @param normal barcode samples for another class
#' @export
#' @importFrom e1071 tune svm 
#' @importFrom ROCR prediction performance 
#' @importFrom  grDevices rainbow
#' @return a list with AUC value for pairwise pathway 
#' @examples
#' nf <- 60
#' res_class<-svm_classification(TCGA_matrix=score_euc_dist,nfs=nf,
#' normal=colnames(norm[,1:10]),tumour=colnames(tumo[,1:10]))
svm_classification<-function(TCGA_matrix,tumour,normal,nfs){
  #library("e1071")
  #library(ROCR)

  scoreMatrix <- as.data.frame(TCGA_matrix[,3:ncol(TCGA_matrix)])
  scoreMatrix <-as.data.frame(scoreMatrix)
  for( i in 1: ncol(scoreMatrix)){
    scoreMatrix[,i] <- as.numeric(as.character(scoreMatrix[,i]))
  }

  TCGA_matrix[,1] <- gsub(" ", "_", TCGA_matrix[,1])
  d<-sub('_-_Homo_sapiens_*', '', TCGA_matrix[,1])
  #d_pr<-sub(')*', '', DataMatrix[,1])
  
  d_pr<- gsub("(human)", "", d, fixed="TRUE")
  d_pr <- gsub("_", "", d_pr)
  d_pr <- gsub("-", "", d_pr)
  
  TCGA_matrix[,2] <- gsub(" ", "_", TCGA_matrix[,2])
  d2<-sub('_-_Homo_sapiens_(human)*', '', TCGA_matrix[,2])
  d_pr2<- gsub("(human)", "", d2, fixed="TRUE")
  d_pr2 <- gsub("_", "", d_pr2)
  d_pr2 <- gsub("-", "", d_pr2)
  
  PathwaysPair <- paste( as.matrix(d_pr), as.matrix(d_pr2),sep="_" )
  
  rownames(scoreMatrix) <-PathwaysPair

  
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

inter1<-intersect(rownames(Dataset_g1),tab_g1_testing)
#rownames(Dataset_g1)<-Dataset_g1$ID

G1_testing<-Dataset_g1[inter1,]

inter2<-intersect(rownames(Dataset_g3),tab_g3_testing)
#rownames(Dataset_g3)<-Dataset_g3$ID
G3_testing<-Dataset_g3[inter2,]

testing<-rbind(G1_testing,G3_testing)

x <- subset(training, select=-Target)
y <- training$Target
#testing[,2]<-NULL
z<-subset(testing, select=-Target)

zi<-testing$Target

auc.df<-list()
svm_model_after_tune_COMPL<-list()
for( k in 2: ncol(training)){
  print(colnames(training)[k])
  svm_tune <- tune(svm, train.x=x, train.y=y, 
                   kernel="radial", ranges=list(cost=10^(-1:2), gamma=c(.5,1,2)),cross=10)
  #print(svm_tune)
  
  svm_model_after_tune <- svm(Target ~ ., data=training[,c(1,k)], kernel="radial", cost=svm_tune$best.parameters$cost, gamma=svm_tune$best.parameters$gamma,cross=10,probability = TRUE)
  
  
  #svm_model_after_tune <- svm(Target ~ ., data=training[,c(1,k)], kernel="radial", cost=svm_tune$best.parameters[1], gamma=svm_tune$best.parameters[2],cross=10,probability = TRUE)
  #summary(svm_model_after_tune)

  j=k-1
  z2=z[,j]
  z3<-as.data.frame(z2)
  #rownames(z3)<-rownames(z)
  #colnames(z3)<-as.character(paste("X",j,sep = ""))
  colnames(z3)<-colnames(z)[j]
  #classifiersMatrix <- c(classifiersMatrix,svm_model_after_tune)
  pred <- predict(svm_model_after_tune,z3,decision.values=TRUE,cross=10)

  #a<-table(pred,zi)
  svm.roc <- prediction(attributes(pred)$decision.values, zi)
  svm.auc <- performance(svm.roc, 'tpr', 'fpr')

  perf <- performance(svm.roc, "auc")
  auc<-perf@y.values[[1]]
  
  auc.df[[j]]<- auc
  svm_model_after_tune_COMPL[[j]]<-svm_model_after_tune
  
  palette <- as.matrix(rainbow(ncol(z)))
  #print(j)
  if (j >1 & j < 6) {
    plot(svm.auc,col=palette[j], add=TRUE)
    legend('bottomright', colnames(z), 
           lty=1, col=palette, bty='n', cex=.90,pch = 20,ncol=1)
    

  }
  else {
    plot(svm.auc, col=palette[j])

    
  }
  
}
names(auc.df) <- colnames(z)
return(auc.df)
}

