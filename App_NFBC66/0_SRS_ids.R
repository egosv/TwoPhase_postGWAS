
### Read  Data
### relates some IDs that are needed for downlading the sequences (SRA ids)
### this file was obtaned from ncbi website (restricted accesss catalog)
sra_run_table <- read.table("SraRunTable.txt",sep="\t",header = T)
sra_run_table_uniq <- unique(sra_run_table)


# Study accession: phs000276.v2.p1
# Table accession: pht002004.v2.p1
NFBC66_Sample_table <- read.table(gzfile("path/to/files/phs000276.v2.pht002004.v2.p1.NFBC66_Sample.MULTI.txt.gz"),sep="\t",header = T)


#### merge based on BioSample and BioSample.Accession for CTS_SRA records in sample table
SRA_to_download <- merge(sra_run_table_uniq,subset(NFBC66_Sample_table,SAMPLE_USE=="CTS_SRA"),by.x='BioSample',by.y='BioSample.Accession')


#### Save file for later use ----
write.table(unique(SRA_to_download$SRA_Sample),file="SRS_unique.txt", quote = F, row.names = F, col.names = F)
