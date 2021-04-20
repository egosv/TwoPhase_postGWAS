args=(commandArgs(TRUE))
print(args)
for(k in 1:length(args)){
  eval(parse(text=args[[k]]))
}

setwd("/out/path/processed_data/")

### Packages needed ----
## install.packages(c('enrichwith', 'nloptr', 'dfoptim'))
## install.packages("twoPhaseGAS_1.07.tar.gz",repos = NULL, type="source")
library(twoPhaseGAS)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("snpStats")
library(snpStats)
library(iterators)
library(foreach)
library(doSNOW)
library(rlecuyer)

### Read data ----
# Study accession: phs000276.v2.p1
# Table accession: pht002005.v2.p1
NFBC66_Subject_Phenotypes_table <- read.table(gzfile("/path/to/files/phs000276.v2.pht002005.v1.p1.c1.NFBC66_Subject_Phenotypes.GRU.txt.gz"),sep="\t",header = T,na.strings = "X")
summary(NFBC66_Subject_Phenotypes_table)
# length(unique(NFBC66_Subject_Phenotypes_table$dbGaP_Subject_ID))
# length(unique(NFBC66_Subject_Phenotypes_table$SUBJID))


### Prepare phenotype file (same as before) ----
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

#### Load GWAS data ----
### Load results from the GWAS on TG
GWAS_resadj_TG <- read.table(gzfile("NFBC66_GWAS_adj_hg19.TG.assoc.linear"),header = T)
head(GWAS_resadj_TG)

### min P overall (signal in GCKR)
minp <- min(GWAS_resadj_TG$P, na.rm=T)
topSNP <- GWAS_resadj_TG[GWAS_resadj_TG$P==minp & !is.na(GWAS_resadj_TG$P),]

### min SNP in chr 8 (LPL region)
minp_chr8 <- min(GWAS_resadj_TG$P[GWAS_resadj_TG$CHR==8], na.rm=T)
topSNP_chr8 <- GWAS_resadj_TG[GWAS_resadj_TG$CHR==8 & GWAS_resadj_TG$P==minp_chr8 & !is.na(GWAS_resadj_TG$P),]

### min SNP in chr 11 (APOA1 region). 
minp_chr11 <- min(GWAS_resadj_TG$P[GWAS_resadj_TG$CHR==11], na.rm=T)
topSNP_chr11 <- GWAS_resadj_TG[GWAS_resadj_TG$CHR==11 & GWAS_resadj_TG$P==minp_chr11 & !is.na(GWAS_resadj_TG$P),]

### rbind(topSNP,topSNP_chr8,topSNP_chr11)
# CHR        SNP        BP (HG19) A1 TEST NMISS     BETA   STAT         P
#   2  rs1260326  27730940  A  ADD  5285  0.06135  6.566 5.671e-11
#   8 rs10096633  19830921  A  ADD  5300 -0.08974 -5.929 3.236e-09
#  11  rs2000571 116585533  A  ADD  5299  0.04967  4.320 1.591e-05

#### Load required GWAS-SNPs based on above data
GWAS_SNPs <- read.plink("NFBC_dbGaP_20091127_hg19", na.strings = c("-9"), sep = ".", select.snps = c(as.character(topSNP$SNP),as.character(topSNP_chr8$SNP),as.character(topSNP_chr11$SNP)))
# str(GWAS_SNPs)
GWAS_SNPs_ <- 2-as(GWAS_SNPs$genotypes,"numeric")
# colMeans(GWAS_SNPs_,na.rm = T)/2
GWAS_SNPs_df <- cbind(data.frame(FID=rownames(GWAS_SNPs_)),GWAS_SNPs_)

### Load subset of subjects that were sequenced (ID)
seq_subjects <- read.table("Convert_IDs_seq_to_geno_n=4511.txt", header = F)
# head(seq_subjects)

### Merge phenotype, covariate with top SNPs
Phase1_data <- merge(merge(Phenotypes_dat[,c("FID","IID","TG","TG.orig")],Covariates_dat[,c("SUBJID","Sex","OC","PG","SexOCPG")],by.x="FID",by.y="SUBJID"),GWAS_SNPs_df)

Phase1_data <- merge(Phase1_data,seq_subjects[,c("V3","V1","V2")],by.x="FID",by.y="V3",all.x=TRUE)
Phase1_data$Has_CTS <- factor(ifelse(is.na(Phase1_data$V1),"N","Y"))
Phase1_data$SampleID <- paste(Phase1_data$V1,Phase1_data$V2,sep = "_")
Phase1_data <- Phase1_data[,-which(names(Phase1_data)%in%c("V1","V2"))]

### summary(Phase1_data)

### Merge with principal components
evec <- read.table("NFBC_dbGaP_20091127_hg19_pca.eigenvec")
names(evec)[1:2] <- c("FID","IID")
names(evec)[3:ncol(evec)] <- paste0("PC",1:(ncol(evec)-2))

### use first 4 principal components in analysis
pcs <- paste0("PC",1:4)

Phase1_data <- merge(Phase1_data, evec, by=c("FID","IID"))


#### I'm going to focus on TG, note that this phenotype can be stratified as follows:
Phase1_data$S <- sapply(Phase1_data$TG.orig,function(x){
  if(is.na(x)) return(NA)
  if(x<1.8) return(1)
  if(x>=1.8 & x<2.3) return(2)
  if(x>=2.3) return(3)
})

### rbind(topSNP,topSNP_chr8,topSNP_chr11)
# CHR  SNP      BP(HG19) A1 TEST NMISS     BETA   STAT      P
#   2  rs1260326  27730940  A  ADD  5285  0.06135  6.566 5.671e-11
#   8 rs10096633  19830921  A  ADD  5300 -0.08974 -5.929 3.236e-09
#  11  rs2000571 116585533  A  ADD  5299  0.04967  4.320 1.591e-05

### Independent signals in GWAS of TG (Table 2 from Service et al)
# LDL rs268 19813529
# APOA5 rs2266788 116660686
# APOA5 rs3135506 116662407

### MAFs and distributions
colMeans(Phase1_data[,grep("^rs",names(Phase1_data))],na.rm = T)/2

# hist(Phase1_data$TG)

### Phase 2 selection 

### I'm using only the two genes that were actually identified by GWAS
snp = c("rs1260326", "rs10096633")


Phase1_data_nonmiss <- na.omit(Phase1_data[,c("FID","TG",snp,"S","Has_CTS","Sex", "OC", "PG", "SexOCPG", pcs)]) ## This is the GWAS-SNP for GCKR
strataZ <- survival:::strata(Phase1_data_nonmiss[,snp])
Phase1_data_nonmiss$StrataZ <- as.numeric(factor(strataZ,labels=1:length(levels(strataZ))))
names(Phase1_data_nonmiss)[2] <- c("Y")
head(Phase1_data_nonmiss)
dim(Phase1_data_nonmiss)

### select the phase 2 subjects from the ones who has CTS
Phase1_data_nonmiss_wCTS <- Phase1_data_nonmiss[Phase1_data_nonmiss$Has_CTS=="Y",c("FID","Y","S",snp,pcs,"StrataZ","SexOCPG")]
table(Phase1_data_nonmiss_wCTS[,snp[1]],Phase1_data_nonmiss_wCTS$S)
table(Phase1_data_nonmiss_wCTS[,snp[2]],Phase1_data_nonmiss_wCTS$S)


### run GWAS regression to select the effects described in the application (assuming additive effect for the snps)
regformula <- as.formula(paste0("Y~G+",paste0(c(snp,pcs),collapse="+"),"+SexOCPG"))
GWASfit <- glm(as.formula(paste0("Y~",paste0(c(snp,pcs),collapse="+"),"+SexOCPG")),data=Phase1_data_nonmiss_wCTS,family=gaussian)
resy = residuals(GWASfit)
beta0 <- c(coef(GWASfit)[1],0,coef(GWASfit)[-1])
disp0 <- summary(GWASfit)$dispersion


N = nrow(Phase1_data_nonmiss_wCTS)
### change the phase sampling fraction here
sf <- 3
samp.frac = c(0.1, 0.25, 0.5)[sf]
n2 = round(N*samp.frac) 

### I'm selecting only one design for all the regions based on the top 2 snps from GCKR and LPL
sel.prob.bal <- rep(1/9,9)
sel.prob.com <- c(1/6,1/6,1/6,0,0,0,1/6,1/6,1/6)
optMethod <- c("Par-spec","A-opt","D-opt")[1]
Kind <- 2

# RDS
order_resy = order(resy)
phase2_id_rds = c(order_resy[1:(n2/2)], order_resy[(N-(n2/2)+1):N])
R1 <- 1*( 1:N %in% phase2_id_rds )

# give different MAF/LD combinations for the grid search:
maf_range = seq(0.05, 0.35, length.out = 5)
LD_range = seq(-0.5, 0.5, length.out = 7)
Z.df <- data.frame(Z=snp,P_z=colMeans(Phase1_data_nonmiss_wCTS[,snp])/2)
maf_ld_comb <- expand.grid(P_z=Z.df$P_z,P_g=maf_range,LD.r=LD_range)
ind_comb <- rep(NA,nrow(maf_ld_comb))
for( i in 1:nrow(maf_ld_comb) ){
  p_gz0 <- tryCatch(p_gz_func(maf_ld_comb$P_z[i], maf_ld_comb$P_g[i], maf_ld_comb$LD.r[i], "G", "Z"), error=function(e){print(e); NULL}) 
  if( !is.null(p_gz0) ) ind_comb[i] <- i
}
maf_ld_comb <- maf_ld_comb[!is.na(ind_comb),]
maf_ld_comb <- merge(maf_ld_comb,Z.df)
numsce <- nrow(maf_ld_comb)


cl <- makeCluster(parallel:::detectCores(), type="SOCK", outfile="")

registerDoSNOW(cl)
clusterSetupRNGstream(cl, seed=rep(148399,6))

ptm <- Sys.time()

fitnessevals <- foreach(i=1:numsce, .combine='c', .inorder=F, .noexport=NULL, .verbose=T, .packages = c("twoPhaseGAS")) %dopar% {

  P_z <- maf_ld_comb[i,1]
  P_g <- maf_ld_comb[i,2]
  LD.r <- maf_ld_comb[i,3]
  Z_ <- as.character(maf_ld_comb[i,4])
  
  err <- 0
  ### calculate p_gz using MAF_g, LD and available data from Z
  p_gz0 <- tryCatch(p_gz_func(P_z, P_g, LD.r, "G", Z_), error=function(e){print(e); err <<- 1; return(NULL)}) 
  
  ## error objects for the optimal allocations
  R3.err <- 0; R4.err <-0
  
  auxformula <- as.formula(paste0("~",Z_))
  
  ### I'm looking at joint covariate and trait allocations for now (use the same allocation for null and alternative SNPs)    
  IM <- obsIM.R(formula=regformula,miscov=~G,auxvar=auxformula,family=gaussian,Phase1_data_nonmiss_wCTS,beta=beta0,p_gz=p_gz0,disp=disp0)
  
  # TZL 
  ## var G/Z
  p_z <- aggregate(as.formula(paste0("q~",Z_)), data = p_gz0, FUN=sum)
  names(p_z)[2] <- "Freq"
  p_g_z <- merge(p_gz0, p_z, by=Z_)
  p_g_z$qcond <- p_g_z$q/p_g_z$Freq
  p_g_z$E1 <- with(p_g_z, G*qcond)
  p_g_z$E2 <- with(p_g_z, G^2*qcond)
  varG_Z = do.call(rbind,by(p_g_z,INDICES=p_g_z[,Z_],function(x){
    c(Z=x[1,Z_],var=sum(x$E2)-sum(x$E1)^2)
  }))
  resyopt = resy*sqrt(varG_Z[match(Phase1_data_nonmiss_wCTS[,Z_], varG_Z[,"Z"]),"var"])
  order_resyopt = order(resyopt)
  phase2_id_opt = c(order_resyopt[1:(n2/2)], order_resyopt[(N-(n2/2)+1):N])
  R2 <- 1*( 1:N %in% phase2_id_opt )
  
  strataformula <- as.formula(paste0("~",Z_,"+S"))
  ## balanced
  pop1a <- sapply(1:100,function(x){
    sa <- which(BSS(samp.fracs=sel.prob.bal,n2,Phase1_data_nonmiss_wCTS,strataformula)$R==1)
    len <- length(sa)
    if( len==n2 ){ return(sa)
    }else if( len>n2 ){ return( sa[-sample(len, len-n2)] )
    }else return( c(sa,sample((1:N)[-sa],n2-len)) )
  })
  fit1 <- apply(pop1a,2,function(r){
    return(fitnessTP(IM,1*(1:N %in% r),optimMeasure=optMethod,K.idx=Kind))
  })
  
  ## combined
  pop2a <- sapply(1:100,function(x){
    sa <- which(BSS(samp.fracs=sel.prob.com,n2,Phase1_data_nonmiss_wCTS,strataformula)$R==1)
    len <- length(sa)
    if( len==n2 ){ return(sa)
    }else if( len>n2 ){ return( sa[-sample(len, len-n2)] )
    }else return( c(sa,sample((1:N)[-sa],n2-len)) )
  })
  fit2 <- apply(pop2a,2,function(r){
    return(fitnessTP(IM,1*(1:N %in% r),optimMeasure=optMethod,K.idx=Kind))
  })
  
  opt.prop <- tryCatch(optimTP.LM(formula=regformula,miscov=~G,auxvar=auxformula,strata=strataformula,family=gaussian,n=n2,Phase1_data_nonmiss_wCTS,beta=beta0,p_gz=p_gz0,disp=disp0,optimMeasure=optMethod,K.idx=Kind),error=function(e){print(e); cat("Error in optimTP.LM on iteration:",it,", using same sampling fraction for all strata instead.\n");
    R3.err <<- 1;
    stratadf <- data.frame(xtabs(strataformula,data=data))
    stratadf$prR_cond_optim <- n2/N; return(stratadf)})
  
  pop3a <- sapply(1:100,function(x){
    sa <- which(BSS(samp.fracs=opt.prop$prR_cond_optim,n2,Phase1_data_nonmiss_wCTS,strataformula)$R==1)
    len <- length(sa)
    if( len==n2 ){ return(sa)
    }else if( len>n2 ){ return( sa[-sample(len, len-n2)] )
    }else return( c(sa,sample((1:N)[-sa],n2-len)) )
  })
  fit3 <- apply(pop3a,2,function(r){
    return(fitnessTP(IM,1*(1:N %in% r),optimMeasure=optMethod,K.idx=Kind))
  })
  R3 <- 1*(1:N %in% pop3a[,which(fit3 == min(fit3) )])
  
  if( (lenrds <- length(phase2_id_rds))!=n2 ){
    if( lenrds>n2 ){
      phase2_id_rds1 <- phase2_id_rds[-sample(lenrds, lenrds-n2)]
    }else {
      phase2_id_rds1 <- c(phase2_id_rds,sample((1:N)[-phase2_id_rds],n2-lenrds))
    }
  }else phase2_id_rds1 <- phase2_id_rds
  
  ## instead of providing a random initialization, give an informed initialization based on combined, balanced and optimal (via Lagrange multipliers) allocations (a third each)
  ## balanced
  pop1 <- pop1a[,order(fit1)[1:20]]
  ## combined
  pop2 <- pop2a[,order(fit2)[1:20]]
  ## LM
  pop3 <- pop3a[,order(fit3)[1:19]]
  
  GA.sol <- tryCatch(optimTP.GA(ncores=1,formula=regformula,miscov=~G,auxvar=auxformula,family=gaussian,n=n2,Phase1_data_nonmiss_wCTS,beta=beta0,p_gz=p_gz0,disp=disp0,ga.popsize=60,ga.propelit=0.9,ga.proptourney=0.9,ga.ngen=500,ga.mutrate=0.001,ga.initpop=t(cbind(pop1,pop2,pop3,phase2_id_rds1)),optimMeasure=optMethod,K.idx=Kind,seed=1),error=function(e){print(e); cat("Error in optimTP.GA on iteration:",it,", using SRS instead.\n");
    R4.err <<- 1;
    return(list(bestsol=sample(N,n2)))})
  # library(kofnGA); plot(GA.sol)
  R4 <- rep(0,N); R4[GA.sol$bestsol] <- 1
  
  indx_i <- list(R2=R2,R3=R3,R4=R4)
  # indx_i <- list(R2=R2,R1=R1)
  FitnessArray <- matrix(NA,3,numsce)
  for( k in 1:numsce ){ # k <- 2
    P_z_b <- maf_ld_comb[k,1]
    P_g_b <- maf_ld_comb[k,2]
    LD.r_b <- maf_ld_comb[k,3]
    Z_b <- as.character(maf_ld_comb[i,4])
    
    p_gz0_b <- p_gz_func(P_z_b, P_g_b, LD.r_b, "G", Z_b)
    
    IM_b <- obsIM.R(formula=regformula,miscov=~G,auxvar=as.formula(paste0("~",Z_b)),family=gaussian,Phase1_data_nonmiss_wCTS,beta=beta0,p_gz=p_gz0_b,disp=disp0)
    
    FitnessArray[1,k] <- fitnessTP(IM_b,Rj=R2,optimMeasure=optMethod,K.idx=Kind)
    FitnessArray[2,k] <- fitnessTP(IM_b,Rj=R3,optimMeasure=optMethod,K.idx=Kind)
    FitnessArray[3,k] <- fitnessTP(IM_b,Rj=R4,optimMeasure=optMethod,K.idx=Kind)
  }
  return(list(list(indx_i,FitnessArray)))
}

stopCluster(cl)

FitnessArray <- array(as.numeric(unlist(lapply(fitnessevals,function(x){x[[2]]}))), dim=c(3, numsce, numsce))

indx <- apply(apply(FitnessArray,c(1,3),mean),1,function(x){
  x1 <- which(x==min(x))
  return(x1[1])
})

# R1 # RDS
R2 <- fitnessevals[[indx[1]]][[1]]$R2 # TZL
R3 <- fitnessevals[[indx[2]]][[1]]$R3 # LM
R4 <- fitnessevals[[indx[3]]][[1]]$R4 # GA

### save designs
save(Phase1_data_nonmiss_wCTS,R1,R2,R3,R4,fitnessevals,file=paste0(savedir,"Phase2_indicators_sf=",sf,".RData"))