library(plyr)
library(dplyr)
library(tibble)
library(reshape)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(grid)

theme_nogrid = function(...) theme_bw() + theme(axis.line = element_line(colour = "black"),
                            panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.background = element_blank())
                             


mypalette=c('COAD'= '#FF8C69', 'HGSC'= '#CDB4DB', 'LUAD' = '#1A759F', 'ccRCC' = '#FFD966', 'UCEC' = '#5A189A', 'PDAC' = '#962A13', 'LSCC' = '#91BDFF', 'GBM' = '#52B788', 'HNSCC' = '#B7C47D', 'BRCA' = '#CD6090')


groups=c('MSI','HRD','POLE','APOBEC','UV','Smoking','Aging')
cancers=c('BRCA','HGSC','COAD','ccRCC','UCEC','GBM','HNSCC','LUAD','PDAC','LSCC')


all_count=NULL
for (group in groups){
categ=read.table(paste('outs/',group,'.txt',sep=''),header=TRUE,sep="\t")
categ$Group=group
colnames(categ)[5:7]=c('Total_mutation','Count','Fraction')

all_count=rbind(all_count, categ)
}

all_count$High_candidate=ifelse(all_count$High_candidate=='Yes',1,0)
all_count$Low_candidate=ifelse(all_count$Low_candidate=='Yes',1,0)

samples_high=aggregate(all_count$High_candidate, by=list(all_count$Disease,all_count$Group),FUN='sum')
samples_high$Status='High'
samples_low=aggregate(all_count$Low_candidate, by=list(all_count$Disease,all_count$Group),FUN='sum')
samples_low$Status='Low'

both=rbind(samples_high, samples_low)
colnames(both)[1:3]=c('Cancer','Group','Count')

to_plot=both

to_plot=to_plot[(to_plot$Group=='HRD' & to_plot$Cancer %in% c('BRCA','HGSC')) | to_plot$Group!='HRD',]

to_plot$Count=as.numeric(as.character(unlist(to_plot$Count)))
to_plot=to_plot %>% mutate(Status=factor(to_plot$Status, levels=c('Low','High'))) %>% arrange(Status)
to_plot$group_status = paste(to_plot$Group,' ',to_plot$Status,sep="")
to_plot$group_status=factor(to_plot$group_status,levels=unique(to_plot$group_status))

colnames(to_plot)[3]=c("Count_cases")

to_plot$Cancer=factor(to_plot$Cancer,levels=rev(c('BRCA','ccRCC','COAD','UCEC','GBM','HNSCC','LSCC','LUAD','HGSC','PDAC')))

#Manually remove HRD classification for all except BRCA and OV:


p <- ggplot(to_plot, aes(x = group_status, y = Cancer)) 

p <- p + geom_point(aes(color = Cancer, size = Count_cases))
 
p <- p + theme(axis.text.x = element_text(colour="black", size=10, angle=45, vjust = 1,hjust=1), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank()) + labs(x="") 

p <- p + facet_grid(.~Group,scales = "free", space = "free") +theme_nogrid() + theme_minimal()

p <- p + theme(axis.title = element_blank(), axis.text.x = element_text(colour="black", size=10, angle=45, vjust = 1,hjust=1), axis.text.y = element_text(colour="black", size=10),axis.ticks = element_blank())

p <- p + scale_size_continuous(range = c(-1, 7), breaks = seq(25, 100, 25)) + scale_color_manual(values=mypalette, guide='none')


pdf("Summary_categories_extremeCasesVers.2023-12-04.pdf", width=9, height=4,useDingbats=FALSE)
p
dev.off()
