
setwd("/out/path/processed_data/")

### Read  Data
### This table relates the different IDs that are used
# Study accession: phs000276.v2.p1
# Table accession: pht002004.v2.p1
# Consent group: All
# Citation instructions: The study accession (phs000276.v2.p1) is used to cite the study and its data tables and documents. The data in this file should be cited using the accession pht002004.v2.p1.
# To cite columns of data within this file, please use the variable (phv#) accessions below:
#
# 1) the table name and the variable (phv#) accessions below; or
# 2) you may cite a variable as phv#.v2.p1

##			phv00129602.v2.p1	phv00129603.v2.p1	phv00196503.v1.p1
NFBC66_Sample_table <- read.table(gzfile("/path/to/files/phs000276.v2.pht002004.v2.p1.NFBC66_Sample.MULTI.txt.gz"),sep="\t",header = T)
head(NFBC66_Sample_table)
table(NFBC66_Sample_table$SAMPLE_USE)
length(unique(NFBC66_Sample_table$BioSample.Accession))
length(unique(NFBC66_Sample_table$SUBJID))
length(unique(NFBC66_Sample_table$SUBJID[NFBC66_Sample_table$SAMPLE_USE=="CTS_SRA"]))
length(unique(NFBC66_Sample_table$SUBJID[duplicated(NFBC66_Sample_table$SUBJID)]))
### Neet to count the number of subjects with both CTS and SNP_Array in this dataset:
NFBC66_Sample_table_CHIP_CTS <- subset(NFBC66_Sample_table,SAMPLE_USE%in%c('SNP_Array','CTS_SRA'))
aggregate(SUBJID~SAMPLE_USE,NFBC66_Sample_table_CHIP_CTS, FUN=function(x){length(unique(x))})
NFBC66_Sample_table_CHIP_CTS_uniqSUBJID <- aggregate(SAMPLE_USE~SUBJID,NFBC66_Sample_table_CHIP_CTS, FUN=length)
dim(NFBC66_Sample_table_CHIP_CTS_uniqSUBJID[NFBC66_Sample_table_CHIP_CTS_uniqSUBJID$SAMPLE_USE==2,])

## Load IDs in the chip data:
IDs_chipdata <- read.table("NFBC_dbGaP_20091127_hg19.fam",header = F)
head(IDs_chipdata)
### All the subjects in the chip are present in the Sample table
# IDs_chipdata$V1[!IDs_chipdata$V1%in%NFBC66_Sample_table$SUBJID]

## Load IDs in the sequence data (after post-processing)
IDs_seqdata <- read.table("PLP.fam",header = F)
IDs_seqdata$Sample_Name <- paste(IDs_seqdata$V1,IDs_seqdata$V2,sep = "_")

IDs_seqdata <- merge(IDs_seqdata,NFBC66_Sample_table,by.x="Sample_Name",by.y="SAMPID")
head(IDs_seqdata)
### 12 subjects appear with sequence data but were not present on the genotype. These will need to be removed from the analysis
# IDs_seqdata[!IDs_seqdata$SUBJID%in%IDs_chipdata$V1,]
# length(IDs_seqdata$SUBJID[!IDs_seqdata$SUBJID%in%IDs_chipdata$V1])

### write to a file to remove these
write.table(IDs_seqdata[!IDs_seqdata$SUBJID%in%IDs_chipdata$V1,c("V1","V2")],file="IDs_to_remove_seqdatan=12.txt", quote = F, row.names = F, col.names = F)


### Write the convertion file for the IDS that are common to both genotype and sequence data
write.table(IDs_seqdata[IDs_seqdata$SUBJID%in%IDs_chipdata$V1, c("V1","V2","SUBJID","SUBJID")],file="Convert_IDs_seq_to_geno_n=4511.txt", quote = F, row.names = F, col.names = F)
