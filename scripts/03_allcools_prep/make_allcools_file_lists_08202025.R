library(stringr)

metadata<-read.csv('/data/combined_wgbs/combined_metadata.csv',header=TRUE)

#AD samples

ad_metadata<-metadata[(metadata$dx=='CTRL' | metadata$dx=='AD') & metadata$study=='ADRC',]

pd_metadata<-metadata[(metadata$dx=='CTRL' | metadata$dx=='PD') & metadata$study=='PDBP',]

lbd_metadata<-metadata[(metadata$dx=='CTRL' | metadata$dx=='DLB'),]

pdbp_filenames<-list.files('/data/pdbp_01302024/aligned_hisat/allc_files',pattern='*.tsv.gz$')
names_pdbp_filenames<-paste0('PDBP',str_replace(str_replace(pdbp_filenames,'allc_',''),'.tsv.gz',''))

pdbp_filenames<-paste0('/data/pdbp_01302024/aligned_hisat/allc_files/',pdbp_filenames)
names(pdbp_filenames)<-names_pdbp_filenames

these_pds<-pd_metadata$X[pd_metadata$X %in% names(pdbp_filenames)]
pd_fl_out<-data.frame(samp=these_pds,fn=pdbp_filenames[these_pds])
write.table(pd_fl_out, file='/data/combined_wgbs/pd_allc_list.tsv', quote=FALSE,sep='\t', row.names=FALSE, col.names=FALSE)


#check for missing allcs

adrc_fastq_samples_P116<-list.files('/data/ADRC_fastq/ADRC_P116',pattern='*_R1.fastq.gz$')
adrc_fastq_samples_P117<-list.files('/data/ADRC_fastq/ADRC_P117',pattern='*_R1.fastq.gz$')
adrc_fastq_samples_Pilot<-list.files('/data/ADRC_fastq/ADRC_Pilot',pattern='*_R1.fastq.gz$')

names_adrc_in<-paste0("ADRC",str_replace(c(adrc_fastq_samples_P116,adrc_fastq_samples_P117,adrc_fastq_samples_Pilot),"_R1.fastq.gz",""))

adrc_filenames<-list.files('/data/ADRC_fastq/ADRC_hisat',pattern='*.tsv.gz$')
names_adrc_filenames<-str_replace(paste0('ADRC',str_replace(str_replace(adrc_filenames,'allc_',''),'.tsv.gz','')),'_','-')
adrc_filenames<-paste0('/data/ADRC_fastq/ADRC_hisat/',adrc_filenames)
names(adrc_filenames)<-names_adrc_filenames

these_ads<-ad_metadata$X[ad_metadata$X %in% names(adrc_filenames)]

ad_fl_out<-data.frame(samp=these_ads,fn=adrc_filenames[these_ads])
write.table(ad_fl_out, file='/data/combined_wgbs/ad_allc_list.tsv', quote=FALSE,sep='\t', row.names=FALSE, col.names=FALSE)

names_adrc_in[!(names_adrc_in %in% metadata$X)]

#lbd cohort

these_lbds<-lbd_metadata$X[lbd_metadata$X %in% c(names(pdbp_filenames),names(adrc_filenames))]

joined_names<-c(pdbp_filenames,adrc_filenames)
names(joined_names)<-c(names(pdbp_filenames),names(adrc_filenames))

lbd_fl_out<-data.frame(samp=these_lbds,fn=joined_names[these_lbds])
write.table(lbd_fl_out, file='/data/combined_wgbs/dlb_allc_list.tsv', quote=FALSE,sep='\t', row.names=FALSE, col.names=FALSE)

