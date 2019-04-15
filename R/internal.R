

Uniform<-function(pathwayfrom){
  mapped_genes<-pathwayfrom
  xx <- unique(unlist(mapped_genes))
  top3 <- matrix(0, length(xx), length(pathwayfrom))
  rownames(top3) <- xx
  colnames(top3)<- names(pathwayfrom)
  for (j in  1:length(xx)){
    for (k in  1:length(pathwayfrom)){
      if (length(intersect(xx[[j]],pathwayfrom[[k]])!=0)){
        
        top3[j,k]<-xx[j]
      }
    }
  }
  top3[top3 == 0] <- " "
return(top3)
  }














delete.NULLs  <-  function(xlist){   # delele null/empty entries in a list
  xlist[unlist(lapply(xlist, nrow) != 0)]
}



#' @title Select the class of TCGA data
#' @description select two labels from ID barcode
#' @param Dataset gene expression matrix
#' @param typesample the labels of the samples (e.g. tumor,normal)
#' @export
#' @return a gene expression matrix of the samples with specified label
#' @examples
#' tumo<-SelectedSample(Dataset=Data_CANCER_normUQ_fil,typesample="tumour")[,2]
SelectedSample <- function(Dataset,typesample){
  if( typesample =="tumour"){
    Dataset <- Dataset[,which( as.numeric(substr(colnames(Dataset), 14, 15)) == 01) ]
  }
  
  if( typesample =="normal"){
    Dataset <- Dataset[,which( as.numeric(substr(colnames(Dataset), 14, 15)) >= 10) ]
  }
  
  return(Dataset)
  
}


#' @title Select the class of TCGA data
#' @description select best performance
#' @param performance_matrix  list of AUC value
#' @param cutoff cut-off for AUC value
#' @return a gene expression matrix with only pairwise pathway with a particular cut-off
select_class<-function(performance_matrix,cutoff){
  tmp_ordered <- as.data.frame(performance_matrix[order(performance_matrix[,1],decreasing=TRUE),])

  er<-tmp_ordered[tmp_ordered[,1]>cutoff,]
  #ase<-tmp_ordered[tmp_ordered$pathway>cutoff,]
  #rownames(er)<-rownames(tmp_ordered)
  #er[,2]<-tmp_ordered$pathway
  #lipid_metabolism<-er[1:length(ase),]
  return(er)
}




IPPIpath_net<-function(pathway,data){
  lista_int<-list()
  for (k in 1:ncol(pathway)){
    print(colnames(pathway)[k])
    currentPathway_genes<-pathway[,k]
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
    names(lista_int)[k]<-colnames(pathway)[k] 
  }
  return(lista_int)
}





IPPIlist_path_net<-function(lista_net,pathway){
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





check_chord <- function(mat, limit){
  
  if(all(colSums(mat) >= limit[2]) & all(rowSums(mat) >= limit[1])) return(mat)
  
  tmp <- mat[(rowSums(mat) >= limit[1]),]
  mat <- tmp[,(colSums(tmp) >= limit[2])]
  
  mat <- check_chord(mat, limit)
  return(mat)
}

bezier <- function(data, process.col){
  x <- c()
  y <- c()
  Id <- c()
  sequ <- seq(0, 1, by = 0.01)
  N <- dim(data)[1]
  sN <- seq(1, N, by = 2)
  if (process.col[1] == '') col_rain <- grDevices::rainbow(N) else col_rain <- process.col
  for (n in sN){
    xval <- c(); xval2 <- c(); yval <- c(); yval2 <- c()
    for (t in sequ){
      xva <- (1 - t) * (1 - t) * data$x.start[n] + t * t * data$x.end[n]
      xval <- c(xval, xva)
      xva2 <- (1 - t) * (1 - t) * data$x.start[n + 1] + t * t * data$x.end[n + 1]
      xval2 <- c(xval2, xva2)
      yva <- (1 - t) * (1 - t) * data$y.start[n] + t * t * data$y.end[n]  
      yval <- c(yval, yva)
      yva2 <- (1 - t) * (1 - t) * data$y.start[n + 1] + t * t * data$y.end[n + 1]
      yval2 <- c(yval2, yva2)			
    }
    x <- c(x, xval, rev(xval2))
    y <- c(y, yval, rev(yval2))
    Id <- c(Id, rep(n, 2 * length(sequ)))
  }
  df <- data.frame(lx = x, ly = y, ID = Id)
  return(df)
}


#' @name GOChord
#' @title Displays the relationship between genes and terms.
#' @description The GOChord function generates a circularly composited overview 
#'   of selected/specific genes and their assigned processes or terms. More 
#'   generally, it joins genes and processes via ribbons in an intersection-like
#'   graph. 
#' @param data The matrix represents the binary relation (1= is related to, 0= 
#'   is not related to) between a set of genes (rows) and processes (columns); a
#'   column for the logFC of the genes is optional
#' @param title The title (on top) of the plot
#' @param space The space between the chord segments of the plot
#' @param gene.order A character vector defining the order of the displayed gene
#'   labels
#' @param gene.size The size of the gene labels
#' @param gene.space The space between the gene labels and the segement of the 
#'   logFC
#' @param nlfc Defines the number of logFC columns (default=1)
#' @param lfc.col The fill color for the logFC specified in the following form: 
#'   c(color for low values, color for the mid point, color for the high values)
#' @param lfc.min Specifies the minimium value of the logFC scale (default = -3)
#' @param lfc.max Specifies the maximum value of the logFC scale (default = 3)
#' @param ribbon.col The background color of the ribbons
#' @param border.size Defines the size of the ribbon borders
#' @param process.label The size of the legend entries
#' @param limit A vector with two cutoff values (default= c(0,0)). 
#' @import ggplot2
#' @import grDevices
GOChord <- function(data, title, space, gene.order, gene.size, gene.space, nlfc = 1, lfc.col, lfc.min, lfc.max, ribbon.col, border.size, process.label, limit){
  y <- id <- xpro <- ypro <- xgen <- ygen <- lx <- ly <- ID <- logFC <- NULL
  Ncol <- dim(data)[2]
  
  if (missing(title)) title <- ''
  if (missing(space)) space = 0
  if (missing(gene.order)) gene.order <- 'none'
  if (missing(gene.size)) gene.size <- 3
  if (missing(gene.space)) gene.space <- 0.2
  if (missing(lfc.col)) lfc.col <- c('brown1', 'azure', 'cornflowerblue')
  if (missing(lfc.min)) lfc.min <- -10
  if (missing(lfc.max)) lfc.max <- 10
  if (missing(border.size)) border.size <- 0.5
  if (missing (process.label)) process.label <- 11
  if (missing(limit)) limit <- c(0, 0)
  
  if (gene.order == 'logFC') data <- data[order(data[, Ncol], decreasing = T), ]
  if (gene.order == 'alphabetical') data <- data[order(rownames(data)), ]
  if (sum(!is.na(match(colnames(data), 'logFC'))) > 0){
    if (nlfc == 1){
      cdata <- check_chord(data[, 1:(Ncol - 1)], limit)
      lfc <- sapply(rownames(cdata), function(x) data[match(x,rownames(data)), Ncol])
    }else{
      cdata <- check_chord(data[, 1:(Ncol - nlfc)], limit)
      lfc <- sapply(rownames(cdata), function(x) data[, (Ncol - nlfc + 1)])
    }
  }else{
    cdata <- check_chord(data, limit)
    lfc <- 0
  }
  if (missing(ribbon.col)) colRib <- grDevices::rainbow(dim(cdata)[2]) else colRib <- ribbon.col
  nrib <- colSums(cdata)
  ngen <- rowSums(cdata)
  Ncol <- dim(cdata)[2]
  Nrow <- dim(cdata)[1]
  colRibb <- c()
  for (b in 1:length(nrib)) colRibb <- c(colRibb, rep(colRib[b], 202 * nrib[b]))
  r1 <- 1; r2 <- r1 + 0.1
  xmax <- c(); x <- 0
  for (r in 1:length(nrib)){
    perc <- nrib[r] / sum(nrib)
    xmax <- c(xmax, (pi * perc) - space)
    if (length(x) <= Ncol - 1) x <- c(x, x[r] + pi * perc)
  }
  xp <- c(); yp <- c()
  l <- 50
  for (s in 1:Ncol){
    xh <- seq(x[s], x[s] + xmax[s], length = l)
    xp <- c(xp, r1 * sin(x[s]), r1 * sin(xh), r1 * sin(x[s] + xmax[s]), r2 * sin(x[s] + xmax[s]), r2 * sin(rev(xh)), r2 * sin(x[s]))
    yp <- c(yp, r1 * cos(x[s]), r1 * cos(xh), r1 * cos(x[s] + xmax[s]), r2 * cos(x[s] + xmax[s]), r2 * cos(rev(xh)), r2 * cos(x[s]))
  }
  df_process <- data.frame(x = xp, y = yp, id = rep(c(1:Ncol), each = 4 + 2 * l))
  xp <- c(); yp <- c(); logs <- NULL
  x2 <- seq(0 - space, -pi - (-pi / Nrow) - space, length = Nrow)
  xmax2 <- rep(-pi / Nrow + space, length = Nrow)
  for (s in 1:Nrow){
    xh <- seq(x2[s], x2[s] + xmax2[s], length = l)
    if (nlfc <= 1){
      xp <- c(xp, (r1 + 0.05) * sin(x2[s]), (r1 + 0.05) * sin(xh), (r1 + 0.05) * sin(x2[s] + xmax2[s]), r2 * sin(x2[s] + xmax2[s]), r2 * sin(rev(xh)), r2 * sin(x2[s]))
      yp <- c(yp, (r1 + 0.05) * cos(x2[s]), (r1 + 0.05) * cos(xh), (r1 + 0.05) * cos(x2[s] + xmax2[s]), r2 * cos(x2[s] + xmax2[s]), r2 * cos(rev(xh)), r2 * cos(x2[s]))
    }else{
      tmp <- seq(r1, r2, length = nlfc + 1)
      for (t in 1:nlfc){
        logs <- c(logs, data[s, (dim(data)[2] + 1 - t)])
        xp <- c(xp, (tmp[t]) * sin(x2[s]), (tmp[t]) * sin(xh), (tmp[t]) * sin(x2[s] + xmax2[s]), tmp[t + 1] * sin(x2[s] + xmax2[s]), tmp[t + 1] * sin(rev(xh)), tmp[t + 1] * sin(x2[s]))
        yp <- c(yp, (tmp[t]) * cos(x2[s]), (tmp[t]) * cos(xh), (tmp[t]) * cos(x2[s] + xmax2[s]), tmp[t + 1] * cos(x2[s] + xmax2[s]), tmp[t + 1] * cos(rev(xh)), tmp[t + 1] * cos(x2[s]))
      }}}
  if(lfc[1] != 0){
    if (nlfc == 1){
      df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:Nrow), each = 4 + 2 * l), logFC = rep(lfc, each = 4 + 2 * l))
    }else{
      df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:(nlfc*Nrow)), each = 4 + 2 * l), logFC = rep(logs, each = 4 + 2 * l))  
    }
  }else{
    df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:Nrow), each = 4 + 2 * l))
  }
  aseq <- seq(0, 180, length = length(x2)); angle <- c()
  for (o in aseq) if((o + 270) <= 360) angle <- c(angle, o + 270) else angle <- c(angle, o - 90)
  df_texg <- data.frame(xgen = (r1 + gene.space) * sin(x2 + xmax2/2),ygen = (r1 + gene.space) * cos(x2 + xmax2 / 2),labels = rownames(cdata), angle = angle)
  df_texp <- data.frame(xpro = (r1 + 0.15) * sin(x + xmax / 2),ypro = (r1 + 0.15) * cos(x + xmax / 2), labels = colnames(cdata), stringsAsFactors = FALSE)
  cols <- rep(colRib, each = 4 + 2 * l)
  x.end <- c(); y.end <- c(); processID <- c()
  for (gs in 1:length(x2)){
    val <- seq(x2[gs], x2[gs] + xmax2[gs], length = ngen[gs] + 1)
    pros <- which((cdata[gs, ] != 0) == T)
    for (v in 1:(length(val) - 1)){
      x.end <- c(x.end, sin(val[v]), sin(val[v + 1]))
      y.end <- c(y.end, cos(val[v]), cos(val[v + 1]))
      processID <- c(processID, rep(pros[v], 2))
    }
  }
  df_bezier <- data.frame(x.end = x.end, y.end = y.end, processID = processID)
  df_bezier <- df_bezier[order(df_bezier$processID,-df_bezier$y.end),]
  x.start <- c(); y.start <- c()
  for (rs in 1:length(x)){
    val<-seq(x[rs], x[rs] + xmax[rs], length = nrib[rs] + 1)
    for (v in 1:(length(val) - 1)){
      x.start <- c(x.start, sin(val[v]), sin(val[v + 1]))
      y.start <- c(y.start, cos(val[v]), cos(val[v + 1]))
    }
  }	
  df_bezier$x.start <- x.start
  df_bezier$y.start <- y.start
  df_path <- bezier(df_bezier, colRib)
  if(length(df_genes$logFC) != 0){
    tmp <- sapply(df_genes$logFC, function(x) ifelse(x > lfc.max, lfc.max, x))
    logFC <- sapply(tmp, function(x) ifelse(x < lfc.min, lfc.min, x))
    df_genes$logFC <- logFC
  }
  
  g<- ggplot() +
    geom_polygon(data = df_process, aes(x, y, group=id), fill='gray70', inherit.aes = F,color='black') +
    geom_polygon(data = df_process, aes(x, y, group=id), fill=cols, inherit.aes = F,alpha=0.6,color='black') +	
    geom_point(aes(x = xpro, y = ypro, size = factor(labels, levels = labels), shape = NA), data = df_texp) +
    guides(size = guide_legend("Pathway", ncol = 4, byrow = T, override.aes = list(shape = 22, fill = unique(cols), size = 8))) +
    theme(legend.text = element_text(size = process.label)) +
    geom_text(aes(xgen, ygen, label = labels, angle = angle), data = df_texg, size = gene.size) +
    geom_polygon(aes(x = lx, y = ly, group = ID), data = df_path, fill = colRibb, color = 'black', size = border.size, inherit.aes = F) +		
    labs(title = title) + theme(axis.line = element_blank(), axis.text.x = element_blank(),
                                axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(),
                                axis.title.y = element_blank(), panel.background = element_blank(), panel.border = element_blank(),
                                panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
  
  
  if (nlfc >= 1){
    g + geom_polygon(data = df_genes, aes(x, y, group = id, fill = logFC), inherit.aes = F, color = 'black') +
      scale_fill_gradient2('score', space = 'Lab', low = lfc.col[3], mid = lfc.col[2], high = lfc.col[1], guide = guide_colorbar(title.position = "top", title.hjust = 0.5), 
                           breaks = c(min(df_genes$logFC), max(df_genes$logFC)), labels = c(round(min(df_genes$logFC)), round(max(df_genes$logFC)))) +
      theme(legend.position = 'bottom', legend.background = element_rect(fill = 'transparent'), legend.box = 'horizontal', legend.direction = 'horizontal')
  }else{
    g + geom_polygon(data = df_genes, aes(x, y, group = id), fill = 'gray50', inherit.aes = F, color = 'black')+
      theme(legend.position = 'bottom', legend.background = element_rect(fill = 'transparent'), legend.box = 'horizontal', legend.direction = 'horizontal')
  }
}


