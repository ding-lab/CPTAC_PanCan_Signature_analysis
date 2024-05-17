#2023-12-19: update using latest markers.

library(plyr)
library(dplyr)
library(tibble)
library(reshape)

library(ggplot2)
library(ggrepel)
library(grid)
library(ggrastr)
library(tidyverse)
options(ggrepel.max.overlaps = Inf)

date='2023-12-19'

options(ggrepel.max.overlaps = Inf)

setwd('/Users/nadezhdaterekhanova/')

############Palette:
palette=c('COAD'= '#FF8C69', 'HGSC'= '#CDB4DB', 'LUAD' = '#1A759F', 'ccRCC' = '#FFD966', 'UCEC' = '#5A189A', 'PDAC' = '#962A13', 'LSCC' = '#91BDFF', 'GBM' = '#52B788', 'HNSCC' = '#B7C47D', 'BRCA' = '#CD6090')


###Load list of cancer-associated genes:
c_genes=read.table('624_Cancer_genes.txt',sep='\t',header=F)

##############################################
####Load proteome/phosphoprosite markers:#####
##############################################

all_w_test_stat=read.table('Markers/Protein_markers.20231204.tsv',header=TRUE,sep="\t")
all_w_test_stat=all_w_test_stat[order(all_w_test_stat$FDR),]
all_w_test_stat$Log_10_FDR=log10(all_w_test_stat$FDR)

all_w_test_stat_phospho=read.table('Markers/Phosphosites_markers.20231204.tsv',header=TRUE,sep="\t")
all_w_test_stat_phospho$Log_10_FDR=log10(all_w_test_stat_phospho$FDR)


#############################################
######Making volcano-plots for the APOBEC:###
#############################################

group="APOBEC"

sel_genes=c('APOBEC3A','APOBEC3B','APOBEC3A_B','APOBEC3G','STAT1') 

results=all_w_test_stat[all_w_test_stat$Group==group,]

#Do capping:
#results$Log_10_FDR=ifelse(-results$Log_10_FDR>10, -10, results$Log_10_FDR)
#results$Log_10_FDR=ifelse(results$Log_10_FDR>10, 10, results$Log_10_FDR)


res_subs_fdr=results[results$FDR<0.05,]

res_subs_fdr=res_subs_fdr[order(-abs(res_subs_fdr$Log2_Fch)),]
res_subs_fdr_gr_1=res_subs_fdr[res_subs_fdr$Gene %in% c(c_genes$V1),][1:20,]
res_subs_fdr_gr_2=res_subs_fdr[res_subs_fdr$Gene %in% c(res_subs_fdr$Gene[1:10],sel_genes),]
res_subs_fdr_gr=rbind(res_subs_fdr_gr_1,res_subs_fdr_gr_2)
res_subs_fdr_gr=res_subs_fdr_gr[!duplicated(res_subs_fdr_gr),]
res_subs_fdr_gr$Direction=ifelse(res_subs_fdr_gr$Log2_Fch>0,"UP","DOWN")

#Remove empty rows with NAs:
res_subs_fdr_gr=res_subs_fdr_gr[!is.na(res_subs_fdr_gr$Gene),]

p = ggplot(data = results,aes(y=-Log_10_FDR,x= Log2_Fch)) + geom_point_rast(alpha=0.3, colour="grey") 

p = p + theme_minimal()+geom_hline(yintercept=-log10(0.05), alpha=0.5,size=0.2) +geom_vline(xintercept=0, alpha=0.5,size=0.2)

p = p + geom_point(data = res_subs_fdr_gr, aes(color=Direction), alpha=0.95)

p = p + ylab("-log10(FDR)")

p = p + xlab("log2(protein expr. fold change (High/Low))") 

p = p + geom_text_repel(data=res_subs_fdr_gr,aes(y=-Log_10_FDR,x= Log2_Fch,label=ifelse(!is.na(Gene), as.character(Gene),NA)),
alpha=1,size=4.3,segment.size  = 0.2,min.segment.length = 0, box.padding = 0.3, color = "black")

p = p + scale_color_manual("Direction",values=c('UP'='red','DOWN'='blue'))

p = p + theme(strip.background = element_blank(),strip.text.x = element_blank(),  strip.text.y = element_text(size = 12))

p = p + theme(legend.text=element_text(size=12),axis.title.x = element_text(size=12), axis.title.y = element_text(size=12), axis.text.x = element_text(colour="black", size=12), 
axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank(),legend.position="none") + xlim(-4,4)


pdf(paste("Volcano-plots/plots/",group,"_Volcanos_",date,"_Top_markers.pdf",sep=""), width=5, height=4,useDingbats=FALSE)
print(p)
dev.off()


#####################
####PHOSPHO - APOBEC#
#####################

group="APOBEC"
sel_genes=c('APOBEC3A','APOBEC3B','APOBEC3A_B','APOBEC3G','STAT1') 

results=all_w_test_stat_phospho[all_w_test_stat_phospho$Group==group,]

#Do capping:
#results$Log_10_FDR=ifelse(-results$Log_10_FDR>5, -5, results$Log_10_FDR)
#results$Log_10_FDR=ifelse(results$Log_10_FDR>5, 5, results$Log_10_FDR)

res_subs_fdr=results[results$FDR<0.05,]
res_subs_fdr=res_subs_fdr[!duplicated(res_subs_fdr$Gene),]


res_subs_fdr=res_subs_fdr[order(-abs(res_subs_fdr$Log2_Fch)),]
res_subs_fdr_gr_1=res_subs_fdr[res_subs_fdr$Gene %in% c(c_genes$V1),][1:10,]
res_subs_fdr_gr_2=res_subs_fdr[res_subs_fdr$Gene %in% c(res_subs_fdr$Gene[1:10],sel_genes),]
res_subs_fdr_gr=rbind(res_subs_fdr_gr_1,res_subs_fdr_gr_2)
res_subs_fdr_gr=res_subs_fdr_gr[!duplicated(res_subs_fdr_gr),]
res_subs_fdr_gr$Direction=ifelse(res_subs_fdr_gr$Log2_Fch>0,"UP","DOWN")

#Remove empty rows with NAs:
res_subs_fdr_gr=res_subs_fdr_gr[!is.na(res_subs_fdr_gr$Gene),]

p = ggplot(data = results,aes(y=-Log_10_FDR,x= Log2_Fch)) 

p = p + geom_point_rast(alpha=0.3, colour="grey") 

p = p + theme_minimal()+geom_hline(yintercept=-log10(0.05), alpha=0.5,size=0.2) 

p = p + geom_vline(xintercept=0, alpha=0.5,size=0.2)+ geom_point(data = res_subs_fdr_gr, aes(color=Direction), alpha=0.95)

p = p + ylab("-log10(FDR)")

p = p + xlab("log2(phosphosite expr. fold change (High/Low))") 

p = p + geom_text_repel(data=res_subs_fdr_gr,aes(y=-Log_10_FDR,x= Log2_Fch,label=ifelse(!is.na(Gene), as.character(Site),NA)),
alpha=1,size=4.3,segment.size  = 0.2,min.segment.length = 0, box.padding = 0.3, color = "black")

p = p + scale_color_manual("Direction",values=c('UP'='red','DOWN'='blue'))

p = p + theme(strip.background = element_blank(),strip.text.x = element_blank(),  strip.text.y = element_text(size = 12))

p = p + theme(legend.text=element_text(size=12),axis.title.x = element_text(size=12), axis.title.y = element_text(size=12), 
axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),
axis.ticks = element_blank(),legend.position="none") + xlim(-2.5,2.5)

pdf(paste("Volcano-plots/plots/Phosphosite_",group,"_Volcanos_",date,"_Top_markers.pdf",sep=""), width=5, height=4,useDingbats=FALSE)
print(p)
dev.off()

###################

################
###MSI##########
################
group="MSI"

sel_genes=c('RPL22L1','MLH1','MRE11','RAD50','NBN')
 
results=all_w_test_stat[all_w_test_stat$Group==group,]

#Do capping:
results$Log_10_FDR=ifelse(-results$Log_10_FDR>10, -10, results$Log_10_FDR)
results$Log_10_FDR=ifelse(results$Log_10_FDR>10, 10, results$Log_10_FDR)


res_subs_fdr=results[results$FDR<0.05,]

res_subs_fdr=res_subs_fdr[order(-abs(res_subs_fdr$Log2_Fch)),]
res_subs_fdr_gr_1=res_subs_fdr[res_subs_fdr$Gene %in% c(c_genes$V1),][1:20,]
res_subs_fdr_gr_2=res_subs_fdr[res_subs_fdr$Gene %in% c(res_subs_fdr$Gene[1:10],sel_genes),]
res_subs_fdr_gr=rbind(res_subs_fdr_gr_1,res_subs_fdr_gr_2)
res_subs_fdr_gr=res_subs_fdr_gr[!duplicated(res_subs_fdr_gr),]
res_subs_fdr_gr$Direction=ifelse(res_subs_fdr_gr$Log2_Fch>0,"UP","DOWN")

#Remove empty rows with NAs:
res_subs_fdr_gr=res_subs_fdr_gr[!is.na(res_subs_fdr_gr$Gene),]

p = ggplot(data = results,aes(y=-Log_10_FDR,x= Log2_Fch)) + geom_point_rast(alpha=0.3, colour="grey") 

p = p + theme_minimal()+geom_hline(yintercept=-log10(0.05), alpha=0.5,size=0.2) +geom_vline(xintercept=0, alpha=0.5,size=0.2)

p = p + geom_point(data = res_subs_fdr_gr, aes(color=Direction), alpha=0.95)

p = p + ylab("-log10(FDR)")

p = p + xlab("log2(protein expr. fold change (High/Low))") 

p = p + geom_text_repel(data=res_subs_fdr_gr,aes(y=-Log_10_FDR,x= Log2_Fch,label=ifelse(!is.na(Gene), as.character(Gene),NA)),
alpha=1,size=4.3,segment.size  = 0.2,min.segment.length = 0, box.padding = 0.3, color = "black")

p = p + scale_color_manual("Direction",values=c('UP'='red','DOWN'='blue'))

p = p + theme(strip.background = element_blank(),strip.text.x = element_blank(),  strip.text.y = element_text(size = 12))

p = p + theme(legend.text=element_text(size=12),axis.title.x = element_text(size=12), axis.title.y = element_text(size=12), axis.text.x = element_text(colour="black", size=12), 
axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank(),legend.position="none") + xlim(-3.2,3.2)


pdf(paste("Volcano-plots/plots/",group,"_Volcanos_",date,"_Top_markers.pdf",sep=""), width=5, height=4,useDingbats=FALSE)
print(p)
dev.off()


##################
#for MSI PHOSPHO:#
##################
group="MSI"

sel_genes=c('RPL22L1','MLH1','MRE11','RAD50','NBN')

results=all_w_test_stat_phospho[all_w_test_stat_phospho$Group==group,]

#Do capping:
results$Log_10_FDR=ifelse(-results$Log_10_FDR>5, -5, results$Log_10_FDR)
results$Log_10_FDR=ifelse(results$Log_10_FDR>5, 5, results$Log_10_FDR)

res_subs_fdr=results[results$FDR<0.05,]
res_subs_fdr=res_subs_fdr[!duplicated(res_subs_fdr$Gene),]


res_subs_fdr=res_subs_fdr[order(-abs(res_subs_fdr$Log2_Fch)),]
res_subs_fdr_gr_1=res_subs_fdr[res_subs_fdr$Gene %in% c(c_genes$V1),][1:10,]
res_subs_fdr_gr_2=res_subs_fdr[res_subs_fdr$Gene %in% c(res_subs_fdr$Gene[1:10],sel_genes),]
res_subs_fdr_gr=rbind(res_subs_fdr_gr_1,res_subs_fdr_gr_2)
res_subs_fdr_gr=res_subs_fdr_gr[!duplicated(res_subs_fdr_gr),]
res_subs_fdr_gr$Direction=ifelse(res_subs_fdr_gr$Log2_Fch>0,"UP","DOWN")

#Remove empty rows with NAs:
res_subs_fdr_gr=res_subs_fdr_gr[!is.na(res_subs_fdr_gr$Gene),]

p = ggplot(data = results,aes(y=-Log_10_FDR,x= Log2_Fch)) 

p = p + geom_point_rast(alpha=0.3, colour="grey") 

p = p + theme_minimal()+geom_hline(yintercept=-log10(0.05), alpha=0.5,size=0.2) 

p = p + geom_vline(xintercept=0, alpha=0.5,size=0.2)+ geom_point(data = res_subs_fdr_gr, aes(color=Direction), alpha=0.95)

p = p + ylab("-log10(FDR)")

p = p + xlab("log2(phosphosite expr. fold change (High/Low))") 

p = p + geom_text_repel(data=res_subs_fdr_gr,aes(y=-Log_10_FDR,x= Log2_Fch,label=ifelse(!is.na(Gene), as.character(Site),NA)),
alpha=1,size=4.3,segment.size  = 0.2,min.segment.length = 0, box.padding = 0.3, color = "black")

p = p + scale_color_manual("Direction",values=c('UP'='red','DOWN'='blue'))

p = p + theme(strip.background = element_blank(),strip.text.x = element_blank(),  strip.text.y = element_text(size = 12))

p = p + theme(legend.text=element_text(size=12),axis.title.x = element_text(size=12), axis.title.y = element_text(size=12), 
axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),
axis.ticks = element_blank(),legend.position="none") + xlim(-2.5,2.5)

pdf(paste("Volcano-plots/plots/Phosphosite_",group,"_Volcanos_",date,"_Top_markers.pdf",sep=""), width=5, height=4,useDingbats=FALSE)
print(p)
dev.off()



#################
#for HRD:########
#################
group="HRD"

sel_genes=c('PARP1')
 
results=all_w_test_stat[all_w_test_stat$Group==group,]

#Do capping:
#results$Log_10_FDR=ifelse(-results$Log_10_FDR>10, -10, results$Log_10_FDR)
#results$Log_10_FDR=ifelse(results$Log_10_FDR>10, 10, results$Log_10_FDR)


res_subs_fdr=results[results$FDR<0.05,]

res_subs_fdr=res_subs_fdr[order(-abs(res_subs_fdr$Log2_Fch)),]
res_subs_fdr_gr_1=res_subs_fdr[res_subs_fdr$Gene %in% c(c_genes$V1),][1:20,]
res_subs_fdr_gr_2=res_subs_fdr[res_subs_fdr$Gene %in% c(res_subs_fdr$Gene[1:10],sel_genes),]
res_subs_fdr_gr=rbind(res_subs_fdr_gr_1,res_subs_fdr_gr_2)
res_subs_fdr_gr=res_subs_fdr_gr[!duplicated(res_subs_fdr_gr),]
res_subs_fdr_gr$Direction=ifelse(res_subs_fdr_gr$Log2_Fch>0,"UP","DOWN")

#Remove empty rows with NAs:
res_subs_fdr_gr=res_subs_fdr_gr[!is.na(res_subs_fdr_gr$Gene),]

p = ggplot(data = results,aes(y=-Log_10_FDR,x= Log2_Fch)) + geom_point_rast(alpha=0.3, colour="grey") 

p = p + theme_minimal()+geom_hline(yintercept=-log10(0.05), alpha=0.5,size=0.2) +geom_vline(xintercept=0, alpha=0.5,size=0.2)

p = p + geom_point(data = res_subs_fdr_gr, aes(color=Direction), alpha=0.95)

p = p + ylab("-log10(FDR)")

p = p + xlab("log2(protein expr. fold change (High/Low))") 

p = p + geom_text_repel(data=res_subs_fdr_gr,aes(y=-Log_10_FDR,x= Log2_Fch,label=ifelse(!is.na(Gene), as.character(Gene),NA)),
alpha=1,size=4.3,segment.size  = 0.2,min.segment.length = 0, box.padding = 0.3, color = "black")

p = p + scale_color_manual("Direction",values=c('UP'='red','DOWN'='blue'))

p = p + theme(strip.background = element_blank(),strip.text.x = element_blank(),  strip.text.y = element_text(size = 12))

p = p + theme(legend.text=element_text(size=12),axis.title.x = element_text(size=12), axis.title.y = element_text(size=12), axis.text.x = element_text(colour="black", size=12), 
axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank(),legend.position="none") + xlim(-2.5,2.5)

pdf(paste("Volcano-plots/plots/",group,"_Volcanos_",date,"_Top_markers.pdf",sep=""), width=5, height=4,useDingbats=FALSE)
print(p)
dev.off()

################
#PHOSPHO - BRCA#
################
group="HRD"
sel_genes=c("PARP1")

results=all_w_test_stat_phospho[all_w_test_stat_phospho$Group==group,]

#Do capping:
#results$Log_10_FDR=ifelse(-results$Log_10_FDR>5, -5, results$Log_10_FDR)
#results$Log_10_FDR=ifelse(results$Log_10_FDR>5, 5, results$Log_10_FDR)

res_subs_fdr=results[results$FDR<0.05,]
res_subs_fdr=res_subs_fdr[!duplicated(res_subs_fdr$Gene),]


res_subs_fdr=res_subs_fdr[order(-abs(res_subs_fdr$Log2_Fch)),]
res_subs_fdr_gr_1=res_subs_fdr[res_subs_fdr$Gene %in% c(c_genes$V1),][1:10,]
res_subs_fdr_gr_2=res_subs_fdr[res_subs_fdr$Gene %in% c(res_subs_fdr$Gene[1:10],sel_genes),]
res_subs_fdr_gr=rbind(res_subs_fdr_gr_1,res_subs_fdr_gr_2)
res_subs_fdr_gr=res_subs_fdr_gr[!duplicated(res_subs_fdr_gr),]
res_subs_fdr_gr$Direction=ifelse(res_subs_fdr_gr$Log2_Fch>0,"UP","DOWN")

#Remove empty rows with NAs:
res_subs_fdr_gr=res_subs_fdr_gr[!is.na(res_subs_fdr_gr$Gene),]

p = ggplot(data = results,aes(y=-Log_10_FDR,x= Log2_Fch)) 

p = p + geom_point_rast(alpha=0.3, colour="grey") 

p = p + theme_minimal()+geom_hline(yintercept=-log10(0.05), alpha=0.5,size=0.2) 

p = p + geom_vline(xintercept=0, alpha=0.5,size=0.2)+ geom_point(data = res_subs_fdr_gr, aes(color=Direction), alpha=0.95)

p = p + ylab("-log10(FDR)")

p = p + xlab("log2(phosphosite expr. fold change (High/Low))") 

p = p + geom_text_repel(data=res_subs_fdr_gr,aes(y=-Log_10_FDR,x= Log2_Fch,label=ifelse(!is.na(Gene), as.character(Site),NA)),
alpha=1,size=4.3,segment.size  = 0.2,min.segment.length = 0, box.padding = 0.3, color = "black")

p = p + scale_color_manual("Direction",values=c('UP'='red','DOWN'='blue'))

p = p + theme(strip.background = element_blank(),strip.text.x = element_blank(),  strip.text.y = element_text(size = 12))

p = p + theme(legend.text=element_text(size=12),axis.title.x = element_text(size=12), axis.title.y = element_text(size=12), 
axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),
axis.ticks = element_blank(),legend.position="none") + xlim(-1.5,1.5)

pdf(paste("Volcano-plots/plots/Phosphosite_",group,"_Volcanos_",date,"_Top_markers.pdf",sep=""), width=5, height=4,useDingbats=FALSE)
print(p)
dev.off()


##############
#Smoking######
##############
group="Smoking"

#sel_genes=c('PRKDC','EGFR')
sel_genes=NULL

results=all_w_test_stat[all_w_test_stat$Group==group,]

#Do capping:
#results$Log_10_FDR=ifelse(-results$Log_10_FDR>10, -10, results$Log_10_FDR)
#results$Log_10_FDR=ifelse(results$Log_10_FDR>10, 10, results$Log_10_FDR)


res_subs_fdr=results[results$FDR<0.05,]

res_subs_fdr=res_subs_fdr[order(-abs(res_subs_fdr$Log2_Fch)),]
res_subs_fdr_gr_1=res_subs_fdr[res_subs_fdr$Gene %in% c(c_genes$V1),][1:20,]
res_subs_fdr_gr_2=res_subs_fdr[res_subs_fdr$Gene %in% c(res_subs_fdr$Gene[1:10],sel_genes),]
res_subs_fdr_gr=rbind(res_subs_fdr_gr_1,res_subs_fdr_gr_2)
res_subs_fdr_gr=res_subs_fdr_gr[!duplicated(res_subs_fdr_gr),]
res_subs_fdr_gr$Direction=ifelse(res_subs_fdr_gr$Log2_Fch>0,"UP","DOWN")

#Remove empty rows with NAs:
res_subs_fdr_gr=res_subs_fdr_gr[!is.na(res_subs_fdr_gr$Gene),]

p = ggplot(data = results,aes(y=-Log_10_FDR,x= Log2_Fch)) + geom_point_rast(alpha=0.3, colour="grey") 

p = p + theme_minimal()+geom_hline(yintercept=-log10(0.05), alpha=0.5,size=0.2) +geom_vline(xintercept=0, alpha=0.5,size=0.2)

p = p + geom_point(data = res_subs_fdr_gr, aes(color=Direction), alpha=0.95)

p = p + ylab("-log10(FDR)")

p = p + xlab("log2(protein expr. fold change (High/Low))") 

p = p + geom_text_repel(data=res_subs_fdr_gr,aes(y=-Log_10_FDR,x= Log2_Fch,label=ifelse(!is.na(Gene), as.character(Gene),NA)),
alpha=1,size=4.3,segment.size  = 0.2,min.segment.length = 0, box.padding = 0.3, color = "black")

p = p + scale_color_manual("Direction",values=c('UP'='red','DOWN'='blue'))

p = p + theme(strip.background = element_blank(),strip.text.x = element_blank(),  strip.text.y = element_text(size = 12))

p = p + theme(legend.text=element_text(size=12),axis.title.x = element_text(size=12), axis.title.y = element_text(size=12), axis.text.x = element_text(colour="black", size=12), 
axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank(),legend.position="none") + xlim(-3,3)


pdf(paste("Volcano-plots/plots/",group,"_Volcanos_",date,"_Top_markers.pdf",sep=""), width=5, height=4,useDingbats=FALSE)
print(p)
dev.off()

######################
#for Smoking PHOSPHO:#
######################

group="Smoking"

#sel_genes=c('PRKDC','EGFR')
sel_genes=NULL

results=all_w_test_stat_phospho[all_w_test_stat_phospho$Group==group,]

#Do capping:
#results$Log_10_FDR=ifelse(-results$Log_10_FDR>5, -5, results$Log_10_FDR)
#results$Log_10_FDR=ifelse(results$Log_10_FDR>5, 5, results$Log_10_FDR)

res_subs_fdr=results[results$FDR<0.05,]
res_subs_fdr=res_subs_fdr[!duplicated(res_subs_fdr$Gene),]


res_subs_fdr=res_subs_fdr[order(-abs(res_subs_fdr$Log2_Fch)),]
res_subs_fdr_gr_1=res_subs_fdr[res_subs_fdr$Gene %in% c(c_genes$V1),][1:10,]
res_subs_fdr_gr_2=res_subs_fdr[res_subs_fdr$Gene %in% c(res_subs_fdr$Gene[1:10],sel_genes),]
res_subs_fdr_gr=rbind(res_subs_fdr_gr_1,res_subs_fdr_gr_2)
res_subs_fdr_gr=res_subs_fdr_gr[!duplicated(res_subs_fdr_gr),]
res_subs_fdr_gr$Direction=ifelse(res_subs_fdr_gr$Log2_Fch>0,"UP","DOWN")

#Remove empty rows with NAs:
res_subs_fdr_gr=res_subs_fdr_gr[!is.na(res_subs_fdr_gr$Gene),]

p = ggplot(data = results,aes(y=-Log_10_FDR,x= Log2_Fch)) 

p = p + geom_point_rast(alpha=0.3, colour="grey") 

p = p + theme_minimal()+geom_hline(yintercept=-log10(0.05), alpha=0.5,size=0.2) 

p = p + geom_vline(xintercept=0, alpha=0.5,size=0.2)+ geom_point(data = res_subs_fdr_gr, aes(color=Direction), alpha=0.95)

p = p + ylab("-log10(FDR)")

p = p + xlab("log2(phosphosite expr. fold change (High/Low))") 

p = p + geom_text_repel(data=res_subs_fdr_gr,aes(y=-Log_10_FDR,x= Log2_Fch,label=ifelse(!is.na(Gene), as.character(Site),NA)),
alpha=1,size=4.3,segment.size  = 0.2,min.segment.length = 0, box.padding = 0.3, color = "black")

p = p + scale_color_manual("Direction",values=c('UP'='red','DOWN'='blue'))

p = p + theme(strip.background = element_blank(),strip.text.x = element_blank(),  strip.text.y = element_text(size = 12))

p = p + theme(legend.text=element_text(size=12),axis.title.x = element_text(size=12), axis.title.y = element_text(size=12), 
axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),
axis.ticks = element_blank(),legend.position="none") + xlim(-2,2)

pdf(paste("Volcano-plots/plots/Phosphosite_",group,"_Volcanos_",date,"_Top_markers.pdf",sep=""), width=5, height=4,useDingbats=FALSE)
print(p)
dev.off()
###############