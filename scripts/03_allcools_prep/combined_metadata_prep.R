

adrc<-read.table('~/adrc_metadata_wPatientID_BCincluded.csv',header=TRUE, sep=',')

pdbp<-read.table('/data/pdbp_01302024/metadata/pdbp_metadata_abbrev_0411.csv',sep=',',quote='', header=TRUE)

dx<-c(adrc$DX,pdbp$DX)
sex<-c(adrc$SEX,pdbp$SEX)



dx[is.na(dx)]<-'None'

dx[dx=="Control"]<-"CTRL"
dx[dx=="Parkinson's Disease"]<-"PD"
dx[dx=="Dementia with Lewy Bodies"]<-"DLB"
dx[dx=="PDD-PD"]<-"PD"

study<-c(rep('ADRC',dim(adrc)[1]),rep('PDBP',dim(pdbp)[1]))

out_table<-data.frame(dx,sex,study,row.names=c(paste0('ADRC',adrc$Sample.ID),paste0('PDBP',pdbp$Sample.num)))

write.table(out_table,file='/data/combined_wgbs/combined_metadata.csv',sep=',',col.names=NA)