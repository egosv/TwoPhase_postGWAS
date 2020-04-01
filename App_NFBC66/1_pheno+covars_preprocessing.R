
### Read data ----
# Study accession: phs000276.v2.p1
# Table accession: pht002005.v2.p1
NFBC66_Subject_Phenotypes_table <- read.table(gzfile("path/to/files/phs000276.v2.pht002005.v1.p1.c1.NFBC66_Subject_Phenotypes.GRU.txt.gz"),sep="\t",header = T,na.strings = "X")
summary(NFBC66_Subject_Phenotypes_table)
# length(unique(NFBC66_Subject_Phenotypes_table$dbGaP_Subject_ID))
# length(unique(NFBC66_Subject_Phenotypes_table$SUBJID))


### Prepare phenotype file ----
# This script generates the plink file with the 9 QTs from Sabatti et al, 2009 
cols <- c("SUBJID","SUBJID","FS_TRIGL","FS_KOL_H","FS_KOL_L","crp3dec","FB_GLUK","FS_INS","BMI","Systolic.Blood.Pressure.MEAN","Diastolic.Blood.Pressure.MEAN")
Phenotypes_dat <- NFBC66_Subject_Phenotypes_table[,cols]
names(Phenotypes_dat) <- c("FID","IID","TG","HDL","LDL","CRP","GLU","INS","BMI","SBP","DBP")

## log-transform some of the traits as described in the paper (page 44, left column bottom)
for( colind in c("TG","BMI","INS","GLU","CRP") ){
  Phenotypes_dat[,paste0(colind,".orig")] <- Phenotypes_dat[,colind]
  Phenotypes_dat[,colind] <- log(Phenotypes_dat[,colind])
}
Phenotypes_dat[!is.finite(Phenotypes_dat[,"CRP"]),"CRP"] <- NA

## save the phenotypes to text file
write.table(Phenotypes_dat,file="NFBC66_phenotypes_raw.txt",row.names = F,quote = F)

### Prepare covariates file ----
cols1 <- c("SUBJID","SUBJID","SEX","Pills31","ZP4202U","ZT20")
# table(NFBC66_Subject_Phenotypes_table$SEX)
Covariates_dat <- NFBC66_Subject_Phenotypes_table[,cols1]
# summary(Covariates_dat)

Covariates_dat$Sex <- factor(ifelse(Covariates_dat$SEX==2,"F","M"))
# table(Covariates_dat$Pills31,Covariates_dat$Sex,useNA = "ifany")
# table(Covariates_dat$ZP4202U,Covariates_dat$Sex,useNA = "ifany")
# table(Covariates_dat$ZT20,Covariates_dat$Sex,useNA = "ifany")
# table(Covariates_dat$Pills31,Covariates_dat$ZP4202U,useNA = "ifany")

## The part below constructs the adjusting covariate in the GWAS

Covariates_dat$OC <- factor(mapply(function(o1,o2){
  if( is.na(o1) & is.na(o2) ) return("U")
  if( (o1==1 & o2==1) | (is.na(o1) & o2==1) ) return("Y")
  if( (o1==0 & o2==0) | (o1==0 & is.na(o2))) return("N")
  else(NA)
},o1=Covariates_dat$Pills31,o2=Covariates_dat$ZP4202U))
# table(Covariates_dat$OC,useNA = "ifany")

Covariates_dat$PG <- factor(mapply(function(p,s){
  if( (p==3 | is.na(p)) & s=="F" ) return("U")
  if( p==1 & !is.na(s) & s=="F") return("Y")
  if( p==2 | s=="M" ) return("N")
},p=Covariates_dat$ZT20,s=Covariates_dat$Sex))
# table(Covariates_dat$PG,useNA = "ifany")

# table(Covariates_dat$OC,Covariates_dat$PG,useNA = "ifany")

Covariates_dat$SexOCPG <- factor(mapply(function(o,p,s){
  if( s=='M' ) return( "1_M" )
  return( paste0(o,"-",p) )
},o=as.character(Covariates_dat$OC),p=as.character(Covariates_dat$PG),s=Covariates_dat$Sex))
# table(Covariates_dat$SexOCPG,useNA = "ifany")

## save the covariates to text file
### Note that the level U-U needs to be removed as it is always constant after removing missing values across phenotypes
write.table(cbind(FID=Covariates_dat$SUBJID,IID=Covariates_dat$SUBJID,model.matrix(~SexOCPG, Covariates_dat)[,-c(1,5)]),file="NFBC66_covar_SexOCPG.txt",row.names = F,quote = F)
