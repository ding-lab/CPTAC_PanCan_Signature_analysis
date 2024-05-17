#2024-01-05: update using marker filtering with Fold Change in the cancer cohort.
#2023-12-19: update analysis using latest markers.

library(plyr)
library(dplyr)
library(tibble)
library(reshape)


library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(tidyverse)
library(data.table)



can="BRCA"
group="HRD"

#Read the list of cancer-associated genes
c_genes=read.table('624_Cancer_genes.txt',sep='\t',header=F)

#Read row and column annotations:
col_annot=read.table(paste('Matrix_annotations/',group,'_',can,'_col_annot.tsv',sep=''),sep='\t',header=T)
row_annot=read.table(paste('Matrix_annotations/',group,'_',can,'_row_annot.tsv',sep=''),sep='\t',header=T)

row_annot$Gene=gsub('(.*)\\_.*','\\1',row_annot$ID)


#Signature annotation (not used right now)
s_annot=read.table(paste('../outs/',group,'.txt',sep=''),sep='\t',header=T)
s_annot=s_annot[s_annot$Disease==can,]
high_c=s_annot$CASE_ID[s_annot$High_candidate=='Yes']
col_annot$Signature_annot=ifelse(col_annot$Case %in% high_c, "High_candidate","Other")



#Read matrix with scaled expression of selected markers:
prot_gexpr_sel=fread(paste('Matrices/',group,'_',can,'_geneExpr_prot_phospo_markers_matrix.txt',sep=""))
prot_gexpr_sel=as.data.frame(prot_gexpr_sel)
prot_gexpr_sel=as.matrix(prot_gexpr_sel)

rownames(prot_gexpr_sel)=prot_gexpr_sel[,1]
prot_gexpr_sel=prot_gexpr_sel[,-1]

class(prot_gexpr_sel)="numeric"


rownames(col_annot)=col_annot$Case
col_annot=col_annot[as.character(colnames(prot_gexpr_sel)),]

####identify ouliers:
#outliers=col_annot$Case[col_annot$Count_MSI_MutSignature<200 & col_annot$Subtype_by_markers=='MSI_high']

col_annot$Fraction=col_annot[,3]/col_annot[,2]

col_annot$empty="NA"

colors = brewer.pal(5, "Set1")
subt_col=colors
names(subt_col)=c("Basal","Her2","LumA","LumB","Normal-like") 


column_ha = HeatmapAnnotation(PAM50Subtype=col_annot$PAM50, MutCount = anno_barplot(col_annot$Mut_count,bar_width = 0.5, height = unit(1, "cm")), MutBRCA1_BRCA2=anno_simple(col_annot$empty,col=c("NA"="white"),pch=18,pt_gp = gpar(col = ifelse(col_annot$Germline==1,"#984ea3","#4daf4a")),pt_size = unit(ifelse(col_annot$Germline==1 | col_annot$Somatic==1,4, NA), "mm")), HRD_score=col_annot$HRD_score,
HRD_MutSign = anno_points(col_annot$Count_HRD_MutSignature, 
#        gp = gpar(col = ifelse(col_annot$Count_HRD_MutSignature >=30,"red","black")), 
        gp = gpar(col = ifelse( col_annot$Signature_annot=='High_candidate',"red","black")), 
        height = unit(1.5, "cm")),
col = list(PAM50Subtype=subt_col,HRD_score=colorRamp2(c(0, 25,50,75, 100), c("#ffffcc","#fed976", "#fd8d3c","#e31a1c", "#800026"))))

split_row=factor(row_annot$Marker_type, levels=c('Gene_expression','Protein_expression','Phosphosite_expression'))
split_column=factor(col_annot$Subtype_by_markers,levels=c('HRD_low','HRD_high'))


genes_annot=row_annot$Gene %in% c(c_genes$V1,"PARP1","BRCA1") 

ht = Heatmap(prot_gexpr_sel, name = paste("GeneExpression",sep=""), top_annotation = column_ha,	col= colorRamp2(c(-1.7, 0, 1.7), c("#0571b0", "white", "#ca0020")),
show_column_names = FALSE,show_parent_dend_line = FALSE, column_split=split_column,row_split=split_row, cluster_column_slices = FALSE,cluster_row_slices = FALSE,show_column_dend = FALSE, 
show_row_dend = FALSE,use_raster = TRUE) + 
Heatmap(genes_annot + 0, name='Genes', col = c("0" = "white", "1" = "red"), show_heatmap_legend = FALSE, width = unit(5, "mm")) +
 rowAnnotation(link = anno_mark(at = which(genes_annot), labels = rownames(prot_gexpr_sel)[genes_annot],labels_gp = gpar(fontsize = 10)))

pdf(paste('Heatmap_plots/Heatmap_',group,'_',can,'.20240105.pdf',sep=""), width=10.4, height=8,useDingbats=FALSE)
print(ht)
for (slice in 1:2){
for(an in c('PAM50Subtype','MutBRCA1_BRCA2','HRD_score')) {
    decorate_annotation(an, {
        grid.rect(gp = gpar(fill = NA, col = "black"))
    },slice=slice)
}
}
dev.off()

############################################################END
