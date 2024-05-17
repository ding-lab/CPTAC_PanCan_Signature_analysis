#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp

library(plyr)
library(dplyr)

library(ggplot2)
library(ggrepel)
library(reshape)
library(reshape2)
library(data.table)



#################MAIN CODE:################

#Read cases list:

samples_tab=read.table('GO_case_list_CPTAC_germline.tsv',
sep='\t',header=T)
colnames(samples_tab)=c('Sample','Cancer')

#Read id mapping table:

ids_map=read.table('/diskmnt/Projects/TSAP/data/harmonized_data_freeze/mapping_metadata.tsv.gz',sep='\t',header=T)


#Read the list of protein-coding genes:

gene_list=read.table('DATA/gencode_gene_info_v22/gencode.gene.info.v22.tsv',
sep='\t',header=T)
gene_list=gene_list[gene_list$gene_type=='protein_coding',]


#Prepare assay data:

path_f=paste('/diskmnt/Projects/TSAP/data/harmonized_data_freeze/Expression_data/Gene_expression/ALL_RNA-Seq_Expr_WashU_FPKM_UQ.tsv.gz',
sep='')

ph=fread(path_f)
ph=as.data.frame(ph)

#Filter by gene_id, so that we get the 19,814 same gene_id genes.
ph=ph[ph$gene_id %in% gene_list$gene_id,]


gene_tab=ph[,c(1,2)]
ph=ph[,c(3:ncol(ph))]

ph=ph[,colnames(ph) %in% ids_map$RNA_Sample_ID]

#Check this:
#ids_map[!(ids_map$RNA_Sample_ID %in% colnames(ph)),]
#NOTE: one sample data is missing in RNA-seq! 03BR011

ids_map_s=ids_map[ids_map$RNA_Sample_ID %in% colnames(ph),]
rownames(ids_map_s)=ids_map_s$RNA_Sample_ID
ids_map_s=ids_map_s[colnames(ph),]

colnames(ph)=ids_map_s$CASE_ID

ids_map_s1=ids_map_s[,c('CASE_ID','cohort')]
colnames(ids_map_s1)=c('Sample','Cancer')

#We will use row N as an ID:
tab=reshape2::melt(as.matrix(ph))
colnames(tab)=c('ID','Sample','Expr')

#No NAs in the expression table, so no need. Do we need to pre-filter by any other way? Decided not to do it, but filter by protein-coding genes.
#tab=tab[!is.na(tab$Expr),]  

tab=merge(tab,ids_map_s1,all.x=T)
gene_tab$ID=rownames(gene_tab)

tab=merge(tab,gene_tab,all.x=T)

#Do transformation log2(x+1):
tab$log2_Expr=log2(tab$Expr+1)

write.table(tab, 'inputs/Gene_expr_perSample.20231113.tsv',sep='\t',quote=F,row.names=F)
tab=read.table('inputs/Gene_expr_perSample.20231113.tsv',sep='\t',header=T)

##################################
###Now, looking for the markers:
groups=c('UV','Smoking','Aging','MSI','HRD','POLE','APOBEC')
cancers=c('BRCA','HGSC','COAD','ccRCC','LUAD','UCEC','HNSCC','GBM','LSCC','PDAC')

all_test_stat=NULL
for (group in groups){

print(group)

categ=read.table(paste('DATA/MutSignatures/',group,
'.txt',sep=''),header=TRUE,sep="\t")

categ=as.data.frame(categ)
#categ=categ[categ$case_ID %in% samples_tab$Sample,]

categ$Status=NA
categ$Status=ifelse(categ$High_candidate=='Yes','High',categ$Status)
categ$Status=ifelse(categ$Low_candidate=='Yes','Low',categ$Status)
colnames(categ)[2]='Sample'
categ_1=categ[,c('Sample','Status')]

#We have category assigned for 1,062 out of 1,064 samples (for all with available mutational signatures)

tab_1=tab[tab$Sample %in% categ_1$Sample,]
tab_2=merge(tab_1,categ_1)
tab_3=tab_2[!is.na(tab_2$Status),]

library(doParallel)
registerDoParallel(cores=40)

#Estimate number that should be processed by a single core, when using 20 cores:
all_unique_ids=unique(tab_3$gene_id)
n_per_core=round(length(all_unique_ids)/40)+1


all_st=NULL
all_st<-foreach(ids_set_n=c(1:40)) %dopar% {
    first_id=(ids_set_n-1)*n_per_core+1
    last_id=min((ids_set_n)*n_per_core,length(all_unique_ids))
    sel_ids=all_unique_ids[first_id:last_id]
    tab_4=tab_3[tab_3$gene_id %in% sel_ids,]

all_st_per_core=NULL
    for (id in unique(tab_4$gene_id)){
    	tab_s=tab_4[tab_4$gene_id==id,]
    	gene=unique(tab_s$gene_name)
    	g1=tab_s$log2_Expr[tab_s$Status=='Low']
    	g2=tab_s$log2_Expr[tab_s$Status=='High']
#    	if (length(g1)>5 & length(g2)>5 & length(unique(tab_s$Cancer))>1){
#          test=glm(tab_s$Expr~tab_s$Cancer+tab_s$PC1+tab_s$PC2+tab_s$PC3+tab_s$Status_1)
       	   test=glm(tab_s$log2_Expr~tab_s$Cancer+tab_s$Status)
       	   st=summary(test)$coefficients
       	   st_1=st[(nrow(st)),4]
       	   f_ch=mean(g2)-mean(g1) #Fold change should be calculated as the difference, as the values are in log2
       	   stat=cbind(id,gene,mean(g2),mean(g1),f_ch,st_1)
       	   all_st_per_core=rbind(all_st_per_core,stat)
#       }
    }
return(all_st_per_core)
}

all_st_f=NULL
for (i in 1:40){
    all_st_1=as.data.frame(all_st[[i]])
    all_st_f=rbind(all_st_f,all_st_1)
}


if (nrow(all_st_f)>0){
all_st_f=as.data.frame(all_st_f)
colnames(all_st_f)=c('ID','Gene','Mean_High_gr','Mean_Low_gr','Log2_Fch','P_value')
all_st_f$FDR=p.adjust(all_st_f$P_value,method='fdr')
all_st_f$Group=group

#Combine for all groups together:
all_test_stat=rbind(all_test_stat,all_st_f)
}
}

all_test_stat=all_test_stat[order(all_test_stat$FDR),]
write.table(all_test_stat, 'out/Gene_expresion_markers.20231204.tsv',sep='\t',quote=F,row.names=F)

###END