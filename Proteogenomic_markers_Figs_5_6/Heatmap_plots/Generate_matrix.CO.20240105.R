#2024-01-05: update using marker filtering with Fold Change in the cancer cohort.
#2023-12-19: update analysis using latest markers.

library(plyr)
library(dplyr)
library(reshape)
library(data.table)

library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(tidyverse)
library(ConsensusClusterPlus)

############
#Main Body:#
############

can="COAD"
can_s="CO"
group="MSI"


#NOTE: one sample data is missing in RNA-seq! 03BR011
#Read id mapping table:
ids_map=read.table('Expression_data/mapping_metadata.tsv.gz',sep='\t',header=T)
ids_map=ids_map[ids_map$CASE_ID!='03BR011' & ids_map$cohort==can,]

#Remove POLE sample; other POLE samples are among UCEC and also 1 OV.
ids_map=ids_map[ids_map$CASE_ID!='11CO070',]

#Read the list of protein-coding genes:
gene_list=read.table('gencode_gene_info_v22/gencode.gene.info.v22.tsv',
sep='\t',header=T)
gene_list=gene_list[gene_list$gene_type=='protein_coding',]


#Read Fold changes:
expr_fch=read.table('Expression_data/Fold_changes/Gene_expresion_Fch_byCancer_byGroup.20231221.tsv',sep='\t',header=T)
expr_fch=expr_fch[expr_fch$can==can & expr_fch$group==group,]
expr_fch=expr_fch[,c('gene_id','gene_name','f_ch')]
colnames(expr_fch)=c('ID','Gene','Log2_Fch_cohort')

prot_fch=read.table('Expression_data/Fold_changes/Protein_expresion_Fch_byCancer_byGroup.20231221.tsv',sep='\t',header=T)
prot_fch=prot_fch[prot_fch$can==can & prot_fch$group==group,]
prot_fch=prot_fch[,c('gene_name','f_ch')]
colnames(prot_fch)=c('Gene','Log2_Fch_cohort')

ph_fch=read.table('Expression_data/Fold_changes/Phosphosite_expresion_Fch_byCancer_byGroup.20231221.tsv',sep='\t',header=T)
ph_fch=ph_fch[ph_fch$can==can & ph_fch$group==group,]
ph_fch=ph_fch[,c('site_id','site','gene_name','f_ch')]
colnames(ph_fch)=c('ID','Site','Gene','Log2_Fch_cohort')



#Read markers:
expr_m=read.table('Markers/Gene_expresion_markers.20231204.tsv',header=TRUE,sep="\t")
expr_m=expr_m[expr_m$Group==group & expr_m$FDR<0.05,]
expr_m=merge(expr_m,expr_fch)
expr_m$Abs_Log2_Fch_cohort=abs(expr_m$Log2_Fch_cohort)
expr_m=expr_m[order(-expr_m$Abs_Log2_Fch_cohort),]


prot_m=read.table('Markers/Protein_markers.20231204.tsv',header=T,sep='\t')
prot_m=prot_m[prot_m$Group==group & prot_m$FDR<0.05,]
prot_m=merge(prot_m,prot_fch)
prot_m$Abs_Log2_Fch_cohort=abs(prot_m$Log2_Fch_cohort)
prot_m=prot_m[order(-prot_m$Abs_Log2_Fch_cohort),]


phospho_m=read.table('Markers/Phosphosites_markers.20231204.tsv',header=T,sep='\t')
phospho_m=phospho_m[phospho_m$Group==group & phospho_m$FDR<0.05,]
phospho_m=merge(phospho_m,ph_fch)
phospho_m$Abs_Log2_Fch_cohort=abs(phospho_m$Log2_Fch_cohort)
phospho_m=phospho_m[order(-phospho_m$Abs_Log2_Fch_cohort),]


###Read annotations:
#######Generate matrice for column-annotations, and check cases ids:
sign=read.table(paste(group,'.txt',sep=''),sep='\t',header=T)
sign=sign[sign$Disease==can,]

sign=sign[,c(2,5,6)]
colnames(sign)=c('Case','Mut_count',paste('Count_',group,'_MutSignature',sep=''))
rownames(sign)=sign$Case


msi=read.table('DATA/MSI_score/MSI_status_All_1063_samples.20230818.tsv',sep='\t',header=T)
msi=msi[msi$Cancer==can_s,]
msi$MSI_status=ifelse(msi$Sites_fraction>=3.5,"MSI-High","MSS")



msi_genes=c('PMS1','PMS2','MLH1','MSH2','MSH3','MSH6')

#Also need to import germline and somatic mutations information:
#Germline:
g_path=read_delim('DATA/Germline_var/CPTAC.PanCan.charged.manuallyReviewed.cancerLabeled.gencode34.v1.1.tsv',delim='\t')
g_path=as.data.frame(g_path)


g_path=g_path[g_path$Cancer_Predisposition_Gene=='TRUE' & g_path$ManualReview=='PASS',]
g_path=g_path[g_path$Sample %in% ids_map$CASE_ID & g_path$HUGO_Symbol %in% msi_genes & g_path$Overall_Classification!='Prioritized VUS' & g_path$Disease==can_s,]


#Somatic:
som=fread('DATA/Somatic_MAF/PanCan_Union_Maf_Broad_WashU_v1.1.Lite.20230718.maf')
som=as.data.frame(som)
som_s=som[,c('Hugo_Symbol','Variant_Classification','Tumor_Sample_Barcode','COHORT')]
som_s=som_s[!(som_s$Variant_Classification %in% c("3'UTR","5'UTR","5'Flank","COULD_NOT_DETERMINE","IGR","Intron","Silent","RNA")),]
som_s$Sample=gsub('(.*)_T$','\\1',som_s$Tumor_Sample_Barcode)

som_s=som_s[som_s$Sample %in% ids_map$CASE_ID & som_s$Hugo_Symbol %in% msi_genes & som_s$COHORT==can,]


#Add annotations of germline and somatic mutations in the MSI-associated genes.
sign$Germline=ifelse(sign$Case %in% g_path$Sample, 1, 0)
sign$Somatic=ifelse(sign$Case %in% som_s$Sample, 1, 0)

msi_s=msi[,c('Sample','Sites_fraction','MSI_status')]
colnames(msi_s)=c('Case','MSI_score','MSI_status')
sign=merge(sign,msi_s)




#####################
#Read gene exression#
#####################


path_f='DATA/Expression_data/Gene_expression/ALL_RNA-Seq_Expr_WashU_FPKM_UQ.tsv.gz'

tab=fread(path_f)
tab=as.data.frame(tab)

tab=tab[tab$gene_id %in% c(expr_m$ID[1:150]),]
rownames(tab)=tab$gene_name
tab=tab[,3:ncol(tab)]


tab=tab[,colnames(tab) %in% ids_map$RNA_Sample_ID]

rownames(ids_map)=ids_map$RNA_Sample_ID
ids_map=ids_map[colnames(tab),]

colnames(tab)=ids_map$CASE_ID

tab=log2(tab+1)
tab=t(tab)
gexpr_sel=scale(tab)


########################
#Read protein exression#
########################

path_f=paste('DATA/Expression_data/Proteome/Proteome_Broad_updated_tumor_NAT_raw_gene_level.tsv.gz', sep='')

tab=fread(path_f)
tab=as.data.frame(tab)

tab=tab[tab$geneSymbol %in% c(prot_m$ID[1:150]),]
rownames(tab)=tab$geneSymbol
tab=tab[,6:ncol(tab)]


tab=tab[,colnames(tab) %in% ids_map$Proteome_Sample_ID]

rownames(ids_map)=ids_map$Proteome_Sample_ID
ids_map=ids_map[colnames(tab),]

colnames(tab)=ids_map$CASE_ID

count_NA=apply(tab, 1, function(x) length(x[is.na(x)]))
fraction_NA=count_NA/ncol(tab)
tab=tab[which(fraction_NA<0.25),]

k <- which(is.na(tab), arr.ind=TRUE)

if (nrow(k)>0){
	tab[k] <- rowMeans(tab, na.rm=TRUE)[k[,1]]
}

tab=t(tab)
prot_sel=scale(tab)


############################
#Read phosphosite exression#
############################

path_f=paste('DATA/Expression_data/Phosphoproteome/Phosphoproteome_Broad_v1_tumor_NAT_raw.tsv.gz',sep='')


tab=fread(path_f)
tab=as.data.frame(tab)

tab=tab[tab$id %in% c(phospho_m$ID[1:150]),]
rownames(tab)=tab$py_sites
tab=tab[,7:ncol(tab)]


tab=tab[,colnames(tab) %in% ids_map$Phosphoproteome_Sample_ID]

rownames(ids_map)=ids_map$Phosphoproteome_Sample_ID
ids_map=ids_map[colnames(tab),]

colnames(tab)=ids_map$CASE_ID

count_NA=apply(tab, 1, function(x) length(x[is.na(x)]))
fraction_NA=count_NA/ncol(tab)
tab=tab[which(fraction_NA<0.25),]

k <- which(is.na(tab), arr.ind=TRUE)

if (nrow(k)>0){
	tab[k] <- rowMeans(tab, na.rm=TRUE)[k[,1]]
}

tab=t(tab)
phospho_sel=scale(tab)




########################
####Bind the tables#####
########################


gexpr_sel=t(gexpr_sel)
prot_sel=t(prot_sel)
phospho_sel=t(phospho_sel)


gexpr_sel=gexpr_sel[,ids_map$CASE_ID]
prot_sel=prot_sel[,ids_map$CASE_ID]
phospho_sel=phospho_sel[,ids_map$CASE_ID]


prot_gexpr_sel=rbind(gexpr_sel,prot_sel,phospho_sel)
rownames(prot_gexpr_sel)=gsub('(.*) ','\\1',rownames(prot_gexpr_sel))

write.table(prot_gexpr_sel,paste('Matrices/',group,'_',can,'_geneExpr_prot_phospo_markers_matrix.txt',sep=""),sep="\t",quote=FALSE)



######Use ConsensuClusterPlus to infer clusters:
title=paste('ConsensusClustering/',group,'_',can,sep='')
results=ConsensusClusterPlus(prot_gexpr_sel,maxK=10,reps=1000,pItem=0.8,pFeature=1,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png",title=title)
mod_kclus=results[[2]][["consensusClass"]]
mod_kclus=as.character(unlist(mod_kclus))
mod_kclus=ifelse(mod_kclus==1, 'MSI_low','MSI_high')

split <- factor(paste0(mod_kclus), levels=c("MSI_low","MSI_high"))
x=as.data.frame(rbind(colnames(prot_gexpr_sel),mod_kclus))
x=as.data.frame(t(x))
colnames(x)=c('Sample','Subtype')
write.table(x,paste('Subtype_annot/',group,'_',can,'markers_based_Subtype.txt',sep=""),sep="\t",row.names=FALSE,quote=FALSE)
###############



col_annot=sign
colnames(x)=c('Case','Subtype_by_markers')
col_annot=merge(col_annot,x)


row_annot=as.data.frame(rownames(prot_gexpr_sel))
colnames(row_annot)='ID'
row_annot$Marker_type=c(replicate(nrow(gexpr_sel),'Gene_expression'), replicate(nrow(prot_sel),'Protein_expression'), replicate(nrow(phospho_sel),'Phosphosite_expression'))


write.table(row_annot,paste('Matrix_annotations/',group,'_',can,'_row_annot.tsv',sep=''),sep='\t',quote=F,row.names=F)
write.table(col_annot,paste('Matrix_annotations/',group,'_',can,'_col_annot.tsv',sep=''),sep='\t',quote=F,row.names=F)



###END