args=(commandArgs(TRUE))
print(args)
for(k in 1:length(args)){
  eval(parse(text=args[[k]]))
}

setwd("/out/path/processed_data/")
  
library(twoPhaseGAS)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("snpStats")
library(snpStats)
library(statmod)
library(iterators)
library(foreach)
library(doSNOW) 

### change the phase sampling fraction here
sf <- 3
ss = c(449,1123,2246)[sf] 
s

load(file=paste0("Phase2_indicators_sf=",sf,".RData"))
# Phase1_data_nonmiss_wCTS,R1,R2,R3,R4,fitnessevals


N = nrow(Phase1_data_nonmiss_wCTS)

## table(R1,R2)
# all.equal(R1,R2)

data_phase2 <- cbind(Phase1_data_nonmiss_wCTS,R1=R1,R2=R2,R3=R3,R4=R4)
data_phase2$SexOCPG = factor(data_phase2$SexOCPG) 

### I'm going to use only the two genes that were actually identified by GWAS
snp = c("rs1260326", "rs10096633")
pcs <- paste0("PC",1:4)

#
formula_Ha <- as.formula(paste0("Y ~ RV +", paste0(c(snp,pcs),collapse="+") ,"+ SexOCPG"))
formula_Ho <- as.formula(paste0("Y ~ ", paste0(c(snp,pcs),collapse="+") ,"+ SexOCPG"))
vars_1 = c("Y","RV","StrataZ",snp,pcs,"S","SexOCPG")
vars_0 = c("Y","StrataZ",snp,pcs,"S","SexOCPG")



cl <- makeCluster(outfile="")
# cl <- makeCluster(parallel::detectCores(), type="SOCK", outfile="")

# library(doParallel); registerDoParallel(cl)
registerDoSNOW(cl)
#clusterSetupRNGstream(cl, seed=rep(148399,6))

ptm <- Sys.time()
RVcomb <- expand.grid(GENE=c("GCKR","LPL","APOA5"), RVthres=c(0.01, 0.001), wg=1:2) 

resRV <- foreach(indx=iter(RVcomb, by='row'), .combine='rbind', .inorder=F, .noexport=NULL, .verbose=T, .packages = c("statmod","twoPhaseGAS","snpStats")) %dopar% {

  GENE <- indx$GENE
  RVthres <- indx$RVthres
  wg <- indx$wg
  
  setwd("/out/path/processed_data/")
  seq_SNPs <- read.plink(paste0(GENE,"_GotCloud_filteredPASS_flt"), na.strings = c("-9"), sep = ".")
  
  ### Select rare variants an low frequency SNPs
  SNPs_sum <- col.summary(seq_SNPs$genotypes)
  SNPs_sum_flt <- SNPs_sum[SNPs_sum$MAF<RVthres,] 
  
  if( (nRV <- NROW(SNPs_sum_flt))==0 ) next
  
  seq_SNPs_flt <- seq_SNPs
  seq_SNPs_flt$genotypes = seq_SNPs_flt$genotypes[,colnames(seq_SNPs_flt$genotypes)%in%rownames(SNPs_sum_flt)]
  seq_SNPs_flt$map = seq_SNPs_flt$map[seq_SNPs_flt$map$snp.name%in%rownames(SNPs_sum_flt),]
  
  ### weighted RV 
  wj <- if( wg==1 ){ 
    ## Madsen and Browning  wj = 1 / [MAFj(1-MAFj)]1/2
    1/sqrt(SNPs_sum_flt$MAF*(1-SNPs_sum_flt$MAF))
  }else{
    ## Wu et al wj = beta(MAFj,a1,a2) with a1,a2=c(1,1) or c(1,25)
    dbeta(SNPs_sum_flt$MAF,1,25) 
  } 
  wjlabel <- c("M&B","WEA")[wg]
  
  RV <- round(rowSums(((2-as(seq_SNPs_flt$genotypes,"numeric"))>0) %*% diag(wj)),1)
  
  dat_sim <- merge(data_phase2,data.frame(FID=rownames(seq_SNPs_flt$genotypes),RV),all.x=TRUE)

  ## 
  corSNP1vsRV <- cor(dat_sim$rs1260326,dat_sim$RV)
  corSNP2vsRV <- cor(dat_sim$rs10096633,dat_sim$RV)
  
  # 1) Complete data analysis 
  dat_complete <- na.omit(dat_sim[,vars_1])
  lmfit.com <- glm(formula_Ha, data=dat_complete, family=gaussian) 
  lmfit.com.Ho <- glm(formula_Ho, data=dat_complete, family=gaussian)
  ### maf
  
  # 2) ML analyses
  for( i in c(1:4)){ # i=1
    cond = eval(parse(text=paste0("dat_sim$R",i,"==1")))
    data1 <- dat_sim[cond,vars_1] # phase 2 data
    data0 <- dat_sim[!cond,vars_0] # phase 2 data complement
    
    resHo <- tryCatch(twoPhaseSPML(formula=formula_Ho,miscov=~RV,auxvar=~StrataZ,family=gaussian,data0,data1,verbose=FALSE), error=function(e){print(e); list(Sobs=NA,ll=NA) })
    resHa <- tryCatch(twoPhaseSPML(formula=formula_Ha,miscov=~RV,auxvar=~StrataZ,family=gaussian,data0,data1,verbose=FALSE), error=function(e){print(e); list(theta=rep(NA,4),var_theta=rep(NA,4),Wobs=NA,ll=NA) })
    
    assign(paste0("res",i), c( resHa$theta[1:4], resHa$var_theta[1:4], resHa$Wobs, resHo$Sobs, 2*(resHa$ll-resHo$ll), NA, NA ) )
  }
  
  ### put results together
  out_com0 <- c(coef(lmfit.com)[1:4], diag(vcov(lmfit.com))[1:4])
  Wcom <- out_com0[2]^2/ diag(vcov(lmfit.com))[2]
  LRcom <- as.numeric(2*(logLik(lmfit.com) - logLik(lmfit.com.Ho)))
  Scom <- glm.scoretest(lmfit.com.Ho,dat_complete[,"RV"])^2
  
  out_com <- c( out_com0, Wcom, Scom, LRcom, corSNP1vsRV, corSNP2vsRV )
  
  names(out_com) <- c(paste0("beta",0:3), paste0("var_beta",0:3), "W","S","LR", "corSNP1vsRV", "corSNP2vsRV")
  
  res <- data.frame(GENE=GENE,RVthres=RVthres,Weight=wjlabel,nRV=nRV,Alloc=c("Complete","RDS","TZL","LM","GA"),rbind(out_com,res1,res2,res3,res4))
  
  return(res)
}
tim2 <- difftime(Sys.time(),ptm,units = "hours")

save(resRV,file=paste0("NFBC66_TG_twophase_RVanalysis_sf=",sf,".RData"))

cat("Done with RV analysis; Elapsed time: ", tim2, " hours\n",sep = "")
stopCluster(cl)
