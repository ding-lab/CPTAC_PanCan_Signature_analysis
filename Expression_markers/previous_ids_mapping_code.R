#colnames(ph) <- gsub('(C3[LN]).([0-9]+)','\\1-\\2',colnames(ph))
#colnames(ph) <- gsub('X(.*)','\\1',colnames(ph))
#ph=ph[,!grepl('\\.N',colnames(ph))]
#ph=ph[,!grepl('^N',colnames(ph))]
#colnames(ph) <- gsub('\\.T','',colnames(ph))
#colnames(ph) <- gsub('^C(.*[CO].*)','\\1',colnames(ph))
#colnames(ph) <- gsub('^C(.*[OV].*)','\\1',colnames(ph))


#colnames(ph) <- gsub('CPT000814','604',colnames(ph))
#colnames(ph) <- gsub('CPT001846','01BR044',colnames(ph))
