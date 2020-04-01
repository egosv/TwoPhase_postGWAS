#  Check for problematic values
args=(commandArgs(TRUE))
print(args)
for(k in 1:length(args)){
  eval(parse(text=args[[k]]))
}

setwd("/out/path/processed_data/")

### change these parameters for each combination accordingly
GENE = c("GCKR", "LPL", "APOA5")[1]
### change the phase sampling fraction here
sf <- 3

library(twoPhaseGAS)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("snpStats")
library(snpStats)

load(file=paste0("Phase2_indicators_sf=",sf,".RData"))

data_phase2 <- cbind(Phase1_data_nonmiss_wCTS,R1=R1,R2=R2,R3=R3,R4=R4)
data_phase2$SexOCPG = factor(data_phase2$SexOCPG) 

### I'm going to use only the two genes that were actually identified by GWAS
snp = c("rs1260326", "rs10096633")
 

### load sequence data and select only the ones with R==1 
seq_SNPs <- read.plink(paste0(GENE,"_GotCloud_filteredPASS_flt"), na.strings = c("-9"), sep = ".")
# seq_SNPs$fam
# head(seq_SNPs$genotypes)

### Filter monomorphic an low frequency SNPs
SNPs_sum <- col.summary(seq_SNPs$genotypes)
SNPs_sum_flt <- SNPs_sum[!SNPs_sum$MAF%in%c(0,0.5),] ### non-monomorphic
seq_SNPs$genotypes = seq_SNPs$genotypes[,colnames(seq_SNPs$genotypes)%in%rownames(SNPs_sum_flt)]
seq_SNPs$map = seq_SNPs$map[seq_SNPs$map$snp.name%in%rownames(SNPs_sum_flt),]
## all.equal(seq_SNPs$map$snp.name,colnames(seq_SNPs$genotypes))
## summary(seq_SNPs$genotypes)
## head(seq_SNPs$map$snp.name)
## head(colnames(seq_SNPs$genotypes))

library(iterators)
library(foreach)
library(doSNOW)

cl <- makeCluster(detectCores()-1, type="SOCK", outfile="")

registerDoSNOW(cl)

ptm <- Sys.time()
# p<-101; G0 <- seq_SNPs$genotypes[,p,drop=FALSE]

formula_Ha <- as.formula(paste0("Y ~ G +", paste0(snp,collapse="+") ,"+ SexOCPG"))
formula_Ho <- as.formula(paste0("Y ~ ", paste0(snp,collapse="+") ,"+ SexOCPG"))
vars_1 = c("Y","G","StrataZ",snp,"S","SexOCPG")
vars_0 = c("Y","StrataZ",snp,"S","SexOCPG")

ptm <- Sys.time()
scanresults <- foreach(G0=iter(seq_SNPs$genotypes,by='column'), p=icount(), .combine='rbind', .inorder=F, .noexport=NULL, .verbose=T, .packages = c("statmod","twoPhaseGAS")) %dopar% {
  
  G = (2-as(G0,"numeric"))
  Gnam <- colnames(G) 
  colnames(G) <- "G"
  dat_sim <- merge(data_phase2,data.frame(FID=rownames(G),G),all.x=TRUE)
  # table(dat_sim$StrataZ,dat_sim$G,dat_sim$R1,useNA="ifany")
  # table(dat_sim$G,dat_sim$R1,useNA="ifany")
  # table(dat_sim$StrataZ,dat_sim$G,dat_sim$R2,useNA="ifany")
  # table(dat_sim$StrataZ,dat_sim$G,dat_sim$R3,useNA="ifany")
  
  # 1) Complete data analysis 
  dat_complete <- na.omit(dat_sim[,vars_1])
  lmfit.com <- glm(formula_Ha, data=dat_complete, family=gaussian) 
  lmfit.com.Ho <- glm(formula_Ho, data=dat_complete, family=gaussian)
  ### maf
  mafG_all = mean(dat_complete$G)/2
  
  # 2) SPML analyses
  for( i in c(1:4)){ # i=1
    cond = eval(parse(text=paste0("dat_sim$R",i,"==1")))
    data1 <- dat_sim[cond,vars_1] # phase 2 data
    data0 <- dat_sim[!cond,vars_0] # phase 2 data complement
    
    resHo <- tryCatch(twoPhaseSPML(formula=formula_Ho,miscov=~G,auxvar=~StrataZ,family=gaussian,data0,data1,verbose=FALSE), error=function(e){print(e); list(Sobs=NA,ll=NA) })
    resHa <- tryCatch(twoPhaseSPML(formula=formula_Ha,miscov=~G,auxvar=~StrataZ,family=gaussian,data0,data1,verbose=FALSE), error=function(e){print(e); list(theta=rep(NA,4),var_theta=rep(NA,4),Wobs=NA,ll=NA) })
    maf_Gexp = sum(resHa$qG$G*resHa$qG$q)/2
    
    assign(paste0("res",i), c( resHa$theta[1:4], resHa$var_theta[1:4], resHa$Wobs, resHo$Sobs, 2*(resHa$ll-resHo$ll), maf_Gexp ) )
  }
  
  ### put results together
  out_com0 <- c(coef(lmfit.com)[1:4], diag(vcov(lmfit.com))[1:4])
  Wcom <- out_com0[2]^2/ diag(vcov(lmfit.com))[2]
  LRcom <- as.numeric(2*(logLik(lmfit.com) - logLik(lmfit.com.Ho)))
  Scom <- glm.scoretest(lmfit.com.Ho,dat_complete[,"G"])^2
  
  out_com <- c( out_com0, Wcom, Scom, LRcom, mafG_all )
  
  names(out_com) <- c(paste0("beta",0:3), paste0("var_beta",0:3), "W","S","LR","MAF_G")
  
  res <- data.frame(p=p,SNP=Gnam,Alloc=c("Complete","RDS","TZL","LM","GA"),rbind(out_com,res1,res2,res3,res4))
  
  return(res) 
}
tim2 <- difftime(Sys.time(),ptm,units = "hours")

scanresults = merge(cbind(p=1:nrow(seq_SNPs$map),seq_SNPs$map[,-3]),scanresults,by=c("p"),sort=TRUE)


save(scanresults,file=paste0("NFBC66_TG_twophase_",GENE,"_sf=",sf,".RData"))

cat("Done with  gene=",GENE," and sf=",sf,"; Elapsed time: ", tim2, " hours\n",sep = "")
stopCluster(cl)
