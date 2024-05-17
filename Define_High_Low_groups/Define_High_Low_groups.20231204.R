groups=c('MSI','HRD','APOBEC','Smoking','Aging','UV','POLE')

for (group in groups){
	tt1=read.table(paste('inputs/',group,'.txt',sep=""),header=TRUE,sep="\t")

	tt1=as.data.frame(tt1)
	tt1$CASE_ID<-gsub('CPT000814','604', tt1$CASE_ID)
	tt1$CASE_ID<-gsub('CPT001846','01BR044', tt1$CASE_ID)
	categ=tt1
	colnames(categ)[c(5:7)]=c('Total_mutation','Count','Fraction')
	categ=categ[!is.na(categ$Count) & !is.na(categ$Total_mutation),]
	categ$Fraction=categ$Count/categ$Total_mutation


	categ$High_candidate=ifelse(categ$Count>20 & categ$Fraction>0.2 & categ$Total_mutation>20,"Yes","")
	categ$Low_candidate=ifelse(categ$Count<10 & categ$Total_mutation>20,"Yes","")
	
	samples_low=as.character(categ$CASE_ID[categ$Low_candidate=="Yes"])
	samples_high=as.character(categ$CASE_ID[categ$High_candidate=="Yes"])
	tt1$High_candidate=ifelse(tt1$CASE_ID %in% samples_high, "Yes", "")
	tt1$Low_candidate=ifelse(tt1$CASE_ID %in% samples_low, "Yes", "")

	write.table(tt1,paste("outs/",group,".txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
}