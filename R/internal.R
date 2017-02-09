


select_path_carb<-function(Carbohydrate){
species<-c("- Homo sapiens (human)")  
a<-paste("Glycolysis / Gluconeogenesis", species)
b<-paste("Citrate cycle (TCA cycle)", species)
c<-paste("Pentose phosphate pathway", species)
d<-paste("Pentose and glucuronate interconversions", species)
e<-paste("Fructose and mannose metabolism", species)
f<-paste("Galactose metabolism", species)
g<-paste("Ascorbate and aldarate metabolism", species)
h<-paste("Starch and sucrose metabolism", species)
i<-paste("Amino sugar and nucleotide sugar metabolism", species)
l<-paste("Pyruvate metabolism", species)
m<-paste("Glyoxylate and dicarboxylate metabolism", species)
n<-paste("Propanoate metabolism", species)
o<-paste("Butanoate metabolism", species)
p<-paste("C5-Branched dibasic acid metabolism", species)
q<-paste("Inositol phosphate metabolism", species)
r<-paste("Enzymes", species)
s<-paste("Compounds with biological roles",species)
mer<-c(a,b,c,d,e,f,g,h,i,l,m,n,o,p,q,r,s)
return(mer)
}

select_path_en<-function(Energy){
  species<-c("- Homo sapiens (human)")  
  r<-paste("Oxidative phosphorylation", species)
  s<-paste("Photosynthesis", species)
  t<-paste("Photosynthesis - antenna proteins", species)
  v<-paste("Carbon fixation in photosynthetic organisms", species)
  u<-paste("Carbon fixation pathways in prokaryotes", species)
  z<-paste("Methane metabolism", species)
  aa<-paste("Nitrogen metabolism", species)
  ab<-paste("Sulfur metabolism", species)
  mer<-c(r,s,t,v,u,z,aa,ab)
  return(mer)
}  
  

select_path_lip<-function(Lipid){ 
  species<-c("- Homo sapiens (human)")  
ac<-paste("Fatty acid biosynthesis", species)
ad<-paste("Fatty acid elongation", species)
ae<-paste("Fatty acid degradation", species)
af<-paste("Synthesis and degradation of ketone bodies", species)
ag<-paste("Cutin, suberine and wax biosynthesis", species)
ah<-paste("Steroid biosynthesis", species)
ai<-paste("Primary bile acid biosynthesis", species)
al<-paste("Secondary bile acid biosynthesis", species)
am<-paste("Steroid hormone biosynthesis", species)
an<-paste("Glycerolipid metabolism", species)
ao<-paste("Glycerophospholipid metabolism", species)
ap<-paste("Ether lipid metabolism", species)
aq<-paste("Sphingolipid metabolism", species)
ar<-paste("Arachidonic acid metabolism", species)
as<-paste("Linoleic acid metabolism", species)
at<-paste("alpha-Linolenic acid metabolism", species)
av<-paste("Biosynthesis of unsaturated fatty acids", species)

mer<-c(ac,ad,ae,af,ag,ah,ai,al,am,an,ao,ap,aq,ar,as,at,av)
return(mer)
}




select_path_amn<-function(Aminoacid){ 
  species<-c("- Homo sapiens (human)")  
ac<-paste("Alanine, aspartate and glutamate metabolism", species)
ad<-paste("Glycine, serine and threonine metabolism", species)
ae<-paste("Cysteine and methionine metabolism", species)
af<-paste("Valine, leucine and isoleucine degradation", species)
ag<-paste("Valine, leucine and isoleucine biosynthesis", species)
ah<-paste("Lysine biosynthesis", species)
ai<-paste("Lysine degradation", species)
al<-paste("Arginine biosynthesis", species)
am<-paste("Arginine and proline metabolism", species)
an<-paste("Histidine metabolism", species)
ao<-paste("Tyrosine metabolism", species)
ap<-paste("Phenylalanine metabolism", species)
aq<-paste("Tryptophan metabolism", species)
ar<-paste("Phenylalanine, tyrosine and tryptophan biosynthesis", species)
as<-paste("beta-Alanine metabolism", species)
at<-paste("Taurine and hypotaurine metabolism", species)
av<-paste("Phosphonate and phosphinate metabolism", species)
au<-paste("Selenocompound metabolism", species)
az<-paste("Cyanoamino acid metabolism", species)
a<-paste("D-Glutamine and D-glutamate metabolism", species)
b<-paste("D-Arginine and D-ornithine metabolism", species)
c<-paste("D-Alanine metabolism", species)
d<-paste("Glutathione metabolism", species)

mer<-c(ac,ad,ae,af,ag,ah,ai,al,am,an,ao,ap,aq,ar,as,at,av,au,az,a,b,c,d)
return(mer)
}

select_path_gly<-function(Glybio_met){ 
ac<-paste("N-Glycan biosynthesis", species)
ad<-paste("Various types of N-glycan biosynthesis", species)
ae<-paste("Mucin type O-Glycan biosynthesis", species)
af<-paste("Other types of O-glycan biosynthesis", species)
ag<-paste("Glycosaminoglycan biosynthesis - CS/DS", species)
ah<-paste("Glycosaminoglycan biosynthesis - HS/Hep", species)
ai<-paste("Glycosaminoglycan biosynthesis - KS", species)
al<-paste("Glycosaminoglycan degradation", species)
am<-paste("Glycosylphosphatidylinositol(GPI)-anchor biosynthesis", species)
an<-paste("Glycosphingolipid biosynthesis - lacto and neolacto series", species)
ao<-paste("Glycosphingolipid biosynthesis - globo series", species)
ap<-paste("Glycosphingolipid biosynthesis - ganglio series", species)
aq<-paste("Lipopolysaccharide biosynthesis", species)
ar<-paste("Peptidoglycan biosynthesis", species)
as<-paste("Other glycan degradation", species)
mer<-c(ac,ad,ae,af,ag,ah,ai,al,am,an,ao,ap,aq,ar,as)
return(mer)
}



select_path_cofa<-function(Cofa_vita_met){ 
  species<-c("- Homo sapiens (human)")  
ac<-paste("Thiamine metabolism", species)
ad<-paste("Riboflavin metabolism", species)
ae<-paste("Vitamin B6 metabolism", species)
af<-paste("Nicotinate and nicotinamide metabolism", species)
ag<-paste("Pantothenate and CoA biosynthesis", species)
ah<-paste("Biotin metabolism", species)
ai<-paste("Lipoic acid metabolism", species)
al<-paste("Folate biosynthesis", species)
am<-paste("One carbon pool by folate", species)
an<-paste("Retinol metabolism", species)
ao<-paste("Porphyrin and chlorophyll metabolism", species)
ap<-paste("Ubiquinone and other terpenoid-quinone biosynthesis", species) 	
mer<-c(ac,ad,ae,af,ag,ah,ai,al,am,an,ao,ap)
return(mer)
}

select_path_transc<-function(Transcription){ 
  species<-c("- Homo sapiens (human)")  
ac<-paste("RNA polymerase", species)
ad<-paste("Basal transcription factors", species)
ae<-paste("Spliceosome", species)
af<-paste("Transcription factors", species)
ag<-paste("Transcription machinery", species)
mer<-c(ac,ad,ae,af,ag)
return(mer)
}



select_path_transl<-function(Translation){ 
  species<-c("- Homo sapiens (human)")  
ac<-paste("Ribosome", species)
ad<-paste("Aminoacyl-tRNA biosynthesis", species)
ae<-paste("RNA transport", species)
af<-paste("mRNA surveillance pathway", species)
ag<-paste("Ribosome biogenesis in eukaryotes", species)
ah<-paste("Ribosomal proteins", species)
ai<-paste("Ribosome biogenesis", species)
al<-paste("Transfer RNA biogenesis", species)
am<-paste("Translation factors", species)
mer<-c(ac,ad,ae,af,ag,ah,ai,al,am)
return(mer)
}

select_path_fold<-function(Folding_sorting_and_degradation){ 
  species<-c("- Homo sapiens (human)")  
ac<-paste("Protein export", species)
ad<-paste("Protein processing in endoplasmic reticulum", species)
ae<-paste("SNARE interactions in vesicular transport", species)
af<-paste("Ubiquitin mediated proteolysis", species)
ag<-paste("Sulfur relay system", species)
ah<-paste("RNA degradation", species)
ai<-paste("Chaperones and folding catalysts", species)
al<-paste("SNAREs", species)
am<-paste("Ubiquitin system", species)
an<-paste("Proteasome", species)
mer<-c(ac,ad,ae,af,ag,ah,ai,al,am,an)
return(mer)
}




select_path_repl<-function(Replication_and_repair){ 
  species<-c("- Homo sapiens (human)")  
ac<-paste("DNA replication", species)
ad<-paste("Base excision repair", species)
ae<-paste("Nucleotide excision repair", species)
af<-paste("Mismatch repair", species)
ag<-paste("Homologous recombination", species)
ah<-paste("Non-homologous end-joining", species)
ai<-paste("Fanconi anemia pathway", species)
al<-paste("DNA replication proteins", species)
am<-paste("Chromosome", species)
an<-paste("DNA repair and recombination", species)
ao<-paste("proteins", species)
mer<-c(ac,ad,ae,af,ag,ah,ai,al,am,an,ao)
return(mer)
}



select_path_sign<-function(Signal_transduction){ 
  species<-c("- Homo sapiens (human)")  
a<-paste("Ras signaling pathway", species)
b<-paste("Rap1 signaling pathway", species)
c<-paste("MAPK signaling pathway", species)
d<-paste("ErbB signaling pathway", species)
e<-paste("Wnt signaling pathway", species)
f<-paste("Notch signaling pathway", species)
g<-paste("Hedgehog signaling pathway", species)
h<-paste("TGF-beta signaling pathway", species)
i<-paste("Hippo signaling pathway", species)
l<-paste("VEGF signaling pathway", species)
m<-paste("Jak-STAT signaling pathway", species)
n<-paste("NF-kappa B signaling pathway", species)
o<-paste("TNF signaling pathway", species)
p<-paste("HIF-1 signaling pathway", species)
q<-paste("FoxO signaling pathway", species)
r<-paste("Calcium signaling pathway", species)
s<-paste("Phosphatidylinositol signaling system", species)
t<-paste("Phospholipase D signaling pathway", species)
v<-paste("Sphingolipid signaling pathway", species)
u<-paste("cAMP signaling pathway", species)
z<-paste("cGMP-PKG signaling pathway", species)
ab<-paste("PI3K-Akt signaling pathway", species)
ac<-paste("AMPK signaling pathway", species)
ad<-paste("mTOR signaling pathway", species)
mer<-c(a,b,c,d,e,f,g,h,i,l,m,n,o,p,q,r,s,t,v,u,z,ab,ac,ad)
return(mer)
}


select_path_sign_mol<-function(Signaling_molecules_and_interaction){ 
  species<-c("- Homo sapiens (human)")  
a<-paste("Neuroactive ligand-receptor interaction", species)
b<-paste("Cytokine-cytokine receptor interaction", species)
c<-paste("ECM-receptor interaction", species)
d<-paste("Cell adhesion molecules (CAMs)", species)
mer<-c(a,b,c,d)
return(mer)
}


select_path_transp_ca<-function(Transport_and_catabolism){ 
  species<-c("- Homo sapiens (human)")  
a<-paste("Endocytosis", species)
b<-paste("Phagosome", species)
c<-paste("Lysosome", species)
d<-paste("Peroxisome", species)
e<-paste("Regulation of autophagy", species)
mer<-c(a,b,c,d,e)
return(mer)
}

select_path_cell_grow<-function(Cell_growth_and_death){ 
  species<-c("- Homo sapiens (human)")  
  a<-paste("Cell cycle", species)
b<-paste("Apoptosis", species)
c<-paste("p53 signaling pathway", species)
mer<-c(a,b,c)
return(mer)
}


select_path_cell_comm<-function(Cellular_community){ 
  species<-c("- Homo sapiens (human)")  
  a<-paste("Focal adhesion", species)
b<-paste("Adherens junction", species)
c<-paste("Tight junction", species)
d<-paste("Gap junction", species)
e<-paste("Signaling pathways regulating pluripotency of stem cells ", species)
mer<-c(a,b,c,d,e)
return(mer)
}


select_path_imm_syst<-function(Immune_system){
  species<-c("- Homo sapiens (human)")  
a<-paste("Hematopoietic cell lineage", species)
b<-paste("Complement and coagulation cascades", species)
c<-paste("Platelet activation", species)
d<-paste("Toll-like receptor signaling pathway", species)
e<-paste("Toll and Imd signaling pathway", species)
f<-paste("NOD-like receptor signaling pathway", species)
g<-paste("RIG-I-like receptor signaling pathway", species)
h<-paste("Cytosolic DNA-sensing pathway", species)
i<-paste("Natural killer cell mediated cytotoxicity", species)
l<-paste("Antigen processing and presentation", species)
m<-paste("T cell receptor signaling pathway", species)
n<-paste("B cell receptor signaling pathway", species)
o<-paste("Fc epsilon RI signaling pathway", species)
p<-paste("Fc gamma R-mediated phagocytosis", species)
q<-paste("Leukocyte transendothelial migration", species)
r<-paste("Intestinal immune network for IgA production", species)
s<-paste("Chemokine signaling pathway", species)

mer<-c(a,b,c,d,e,f,g,h,i,l,m,n,o,p,q,r,s)
return(mer)
}




select_path_end_syst<-function(Endocrine_system){ 
  species<-c("- Homo sapiens (human)")  
a<-paste("Insulin secretion", species)
b<-paste("Insulin signaling pathway", species)
c<-paste("Glucagon signaling pathway", species)
d<-paste("Regulation of lipolysis in adipocytes", species)
e<-paste("Adipocytokine signaling pathway", species)
f<-paste("PPAR signaling pathway", species)
g<-paste("GnRH signaling pathway", species)
h<-paste("Ovarian steroidogenesis", species)
i<-paste("Estrogen signaling pathway", species)
l<-paste("Progesterone-mediated oocyte maturation", species)
m<-paste("Prolactin signaling pathway", species)
n<-paste("Oxytocin signaling pathway", species)
o<-paste("Thyroid hormone synthesis", species)
p<-paste("Thyroid hormone signaling pathway", species)
q<-paste("Melanogenesis", species)
r<-paste("Renin secretion", species)
s<-paste("Renin-angiotensin system", species)
t<-paste("Aldosterone synthesis and secretion", species)


mer<-c(a,b,c,d,e,f,g,h,i,l,m,n,o,p,q,r,s,t)
return(mer)
}


select_path_circ_syst<-function(Circulatory_system){ 
  species<-c("- Homo sapiens (human)")  
  a<-paste("Cardiac muscle contraction", species)
b<-paste("Adrenergic signaling in cardiomyocytes", species)
c<-paste("Vascular smooth muscle contraction", species)
mer<-c(a,b,c)
return(mer)
}


select_path_dig_syst<-function(Digestive_system){ 
  species<-c("- Homo sapiens (human)")  
  a<-paste("Salivary secretion", species)
b<-paste("Gastric acid secretion", species)
c<-paste("Pancreatic secretion", species)
d<-paste("Bile secretion", species)
e<-paste("Carbohydrate digestion and absorption", species)
f<-paste("Protein digestion and absorption", species)
g<-paste("Fat digestion and absorption", species)
h<-paste("Vitamin digestion and absorption", species)
i<-paste("Mineral absorption", species)

mer<-c(a,b,c,d,e,f,g,h,i)
return(mer)
}



select_path_exc_syst<-function(Excretory_system){ 
  species<-c("- Homo sapiens (human)")  
  a<-paste("Vasopressin-regulated water reabsorption", species)
b<-paste("Aldosterone-regulated sodium reabsorption", species)
c<-paste("Endocrine and other factor-regulated calcium reabsorption", species)
d<-paste("Proximal tubule bicarbonate reclamation", species)
e<-paste("Collecting duct acid secretion", species)


mer<-c(a,b,c,d,e)
return(mer)
}


select_path_ner_syst<-function(Nervous_system){
  species<-c("- Homo sapiens (human)")  
a<-paste("Glutamatergic synapse", species)
b<-paste("GABAergic synapse", species)
c<-paste("Cholinergic synapse", species)
d<-paste("Dopaminergic synapse", species)
e<-paste("Serotonergic synapse", species)
f<-paste("Long-term potentiation", species)
g<-paste("Long-term depression", species)
h<-paste("Retrograde endocannabinoid signaling", species)
i<-paste("Synaptic vesicle cycle", species)
l<-paste("Neurotrophin signaling pathway", species)

mer<-c(a,b,c,d,e,f,g,h,i,l)
return(mer)
}


select_path_sens_syst<-function(Sensory_system){ 
  species<-c("- Homo sapiens (human)")  
  a<-paste("Phototransduction", species)
b<-paste("Olfactory transduction", species)
c<-paste("Taste transduction", species)
d<-paste("Inflammatory mediator regulation of TRP channels", species)
mer<-c(a,b,c,d)
return(mer)
}



#' @title Select the class of TCGA data
#' @description select two labels from ID barcode
#' @param Dataset gene expression matrix
#' @param typesample the labels of the samples (e.g. tumor,normal)
#' @export
#' @return a gene expression matrix of the samples with specified label
#' @examples
#' tumo<-SelectedSample(Dataset=Data_CANCER_normUQ_filt,typesample="tumor")[,2]
SelectedSample <- function(Dataset,typesample){
  if( typesample =="tumor"){
    Dataset <- Dataset[,which( as.numeric(substr(colnames(Dataset), 14, 15)) == 01) ]
  }
  
  if( typesample =="normal"){
    Dataset <- Dataset[,which( as.numeric(substr(colnames(Dataset), 14, 15)) >= 10) ]
  }
  
  return(Dataset)
  
}


#' @title Select the class of TCGA data
#' @description select two labels from ID barcode
#' @param cutoff cut-off for AUC value
#' @param auc.df list of AUC value
#' @return a gene expression matrix with only pairwise pathway with a particular cut-off
select_class<-function(auc.df,cutoff){
ds<-do.call("rbind", auc.df)
tmp_ordered <- as.data.frame(ds[order(ds,decreasing=TRUE),])
colnames(tmp_ordered)<-'pathway'
er<-as.data.frame(tmp_ordered$pathway>cutoff)
ase<-tmp_ordered[tmp_ordered$pathway>cutoff,]
rownames(er)<-rownames(tmp_ordered)
er[,2]<-tmp_ordered$pathway
lipid_metabolism<-er[1:length(ase),]
return(lipid_metabolism)
}




#' @title Process matrix TCGA data after the selection of pairwise pathway
#' @description processing gene expression matrix
#' @param measure matrix with measure of cross-talk among pathways
#' @param list_perf output of the function select_class 
#' @return a gene expression matrix for case study 1
process_matrix<-function(measure,list_perf){
scoreMatrix <- as.data.frame(measure[,3:ncol(measure)])
for( i in 1: ncol(scoreMatrix)){
  scoreMatrix[,i] <- as.numeric(as.character(scoreMatrix[,i]))
}
measure[,1] <- gsub(" ", "_", measure[,1])
d<-sub('_-_Homo_sapiens_*', '', measure[,1])
d_pr<- gsub("(human)", "", d, fixed="TRUE")
d_pr <- gsub("_", "", d_pr)
d_pr <- gsub("-", "", d_pr)
measure[,2] <- gsub(" ", "_", measure[,2])
d2<-sub('_-_Homo_sapiens_(human)*', '', measure[,2])
d_pr2<- gsub("(human)", "", d2, fixed="TRUE")
d_pr2 <- gsub("_", "", d_pr2)
d_pr2 <- gsub("-", "", d_pr2)
PathwaysPair <- paste( as.matrix(d_pr), as.matrix(d_pr2),sep="_" )
rownames(scoreMatrix) <-PathwaysPair
intera<-intersect(rownames(scoreMatrix),rownames(list_perf))
path_bestlipd<-scoreMatrix[intera,]
return(path_bestlipd)
}



process_matrix_cell_process<-function(measure_cell_process){
score__cell_grow_d <- as.data.frame(measure_cell_process[,3:ncol(measure_cell_process)])
for( i in 1: ncol(score__cell_grow_d)){
  score__cell_grow_d[,i] <- as.numeric(as.character(score__cell_grow_d[,i]))
}

measure_cell_process[,1] <- gsub(" ", "_", measure_cell_process[,1])
d<-sub('_-_Homo_sapiens_*', '', measure_cell_process[,1])

d_pr<- gsub("(human)", "", d, fixed="TRUE")
d_pr <- gsub("_", "", d_pr)
d_pr <- gsub("-", "", d_pr)

measure_cell_process[,2] <- gsub(" ", "_", measure_cell_process[,2])
d2<-sub('_-_Homo_sapiens_(human)*', '', measure_cell_process[,2])
d_pr2<- gsub("(human)", "", d2, fixed="TRUE")
d_pr2 <- gsub("_", "", d_pr2)
d_pr2 <- gsub("-", "", d_pr2)

PathwaysPair <- paste( as.matrix(d_pr), as.matrix(d_pr2),sep="_" )
rownames(score__cell_grow_d) <-PathwaysPair
return(score__cell_grow_d)
}


#' @title Get human KEGG pathway data.
#' @description getKEGGdata creates a data frame with human KEGG pathway. Columns are the pathways and rows the genes inside those pathway 
#' @param mer  output for example of select_path_carb
#' @export
#' @importFrom KEGGREST keggList
#' @return dataframe with human pathway data
proc_path<-function(mer){
pathways.list <- keggList("pathway", "hsa")## returns the list of human pathways
common<-intersect(pathways.list,mer)
lo<-list()
for (i in 1:length(pathways.list)){
  if (length(intersect(pathways.list[[i]],common)!=0)){
    lo[[i]]<-pathways.list[[i]]
    names(lo)[[i]]<-names(pathways.list)[[i]]
  }
}
pathways.list<-lo[lapply(lo,length)!=0] 
pathway.codes <- sub("path:", "", names(pathways.list))
b<-do.call("rbind", pathways.list)
list_pathkegg<-list(pathway.codes,b)
return(list_pathkegg)
}

