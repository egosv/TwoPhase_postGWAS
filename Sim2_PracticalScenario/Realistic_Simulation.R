
### Packages ----
# library(sim1000G)  # optional for data generation 

## required:
library(MASS)
library(Matrix)
library(data.table)
library(statmod)
library(enrichwith)
library(nloptr)
library(dfoptim)
library(aod)
library(parallel)
library(doParallel)

### install and load package here
# install.packages("twoPhaseGAS_1.07.tar.gz", repos=NULL, type="source")
library(twoPhaseGAS)


#### packages for visualization results:
library(reshape)
library(ggplot2)
library(gridExtra)
library(RColorBrewer) 
library(scales) 
library(vcd)


## additionally for parallel computation:
library(iterators)
library(foreach)
library(rlecuyer)


###  Global parameters ----
N = 5000 ## phase 1 sample size
ncores = 4 ## number of cores available for parallel computing, assuming a single computer but can be extended for clusters with multiple nodes, via snow and doSNOW packages

optM <- 1 # Can vary from 1-3 for each optimality criterion. When optM=1, beta1 is a scalar and a two-dimensional vector when is not, so the formula needs to be adapted accordingly.
if( optM==1 ){
  Kind <- 2
  formula_Ha <- Y~G+fZ
  miscov_formula <- ~G
} else{
  Kind <- NULL
  formula_Ha <- Y~G1+G2+fZ
  miscov_formula <- ~G1+G2
} 
Gs <- all.vars(miscov_formula)
### The formula under the null is the same regardless
formula_Ho <- Y~fZ

optMethod <- c("Par-spec","A-opt","D-opt")[optM]
optMetout <- c("Param","Aopt","Dopt")[optM]

sampl <- c("tagSNP","trait","trait-tagSNP")[3] ## either marginal or joint sampling (only joint sampling, i.e."trait-tagSNP", was used in paper)
sfrac <- 0.5 ## determines phase 2 sample size (0.25/0.5 in the paper)
n2 <- round(N*sfrac)

## Stratification formula for the designs and other needed quantities
if ( sampl == "trait" ){
  strataformula <- ~S
  sel.prob.bal <- rep(1/3,3)
  sel.prob.com <- c(1/2,0,1/2)
}else if ( sampl == "tagSNP" ){
  strataformula <- ~Z
  sel.prob.bal <- sel.prob.com <- rep(1/3,3)
}else{
  strataformula <- ~Z+S
  sel.prob.bal <- rep(1/9,9)
  sel.prob.com <- c(1/6,1/6,1/6,0,0,0,1/6,1/6,1/6)
}

save_ind <- TRUE ## logical indicator for whether results for each step should be saved
savedir <- "./results/" ## path where files will be saved

### Data generation ( just for information purposes, this can be safely skipped ) ----
if(0) source('Step_optional_Data_generation.R') ## Produces a new set of genotypes and saves it into "./data_Realistic_R=1K_N=5K.RData". Note that generating new genotypes will not reproduce the results from the paper.


### Step 0 - Load realistic data and other parameters ---- 
load(file=paste0("data_Realistic_R=1K_N=5K.RData"))
dat_sim_all_sub <- subset(dat_sim_all,iteration<=ncores) # For quick run simply filter according to the number of cores, change to 500 for full reproducibility (as used in the paper)
Rep <- max(dat_sim_all_sub$iteration)

### Some data preparation for the design quantities
Zpos <- "56990716"
Z <- genotype[,Zpos]
genotype_noZ <- genotype[,-which(colnames(genotype)==Zpos)]
dist_P_g <- summary(colMeans(genotype_noZ)/2) ## MAF distribution for Gs
cor.mat <- round(cor(genotype),5)
dist_LD.r <- summary(cor.mat[rownames(cor.mat)==Zpos,-which(colnames(genotype)==Zpos)]) ## distribution of LD between each G in the region of interest and Z (removing Z of course)
p_Z <- data.frame(xtabs(~Z)/length(Z)) 

### check if design LD and MAF combinations are plausible.
if( optM==1 ){
  maf_ld_comb <- expand.grid(P_g=dist_P_g[c(2,3,5)],LD.r=dist_LD.r[c(2,3,5)])
  ind_comb <- rep(NA,nrow(maf_ld_comb))
  
  for( i in 1:nrow(maf_ld_comb) ){
    p_gz0 <- tryCatch(p_gz_func(p_Z, maf_ld_comb$P_g[i], maf_ld_comb$LD.r[i], "G", "Z"), error=function(e){print(e); return(tryCatch(p_gz_func(mean(Z)/2, maf_ld_comb$P_g[i], maf_ld_comb$LD.r[i], "G", "Z"),error=function(e){NULL})) })
    if( !is.null(p_gz0) ) ind_comb[i] <- i
  }
}else{
  maf_ld_comb <- expand.grid(P_g1=dist_P_g[c(2,5)],P_g2=dist_P_g[c(2,5)],LD.rZG1=dist_LD.r[c(2,5)],LD.rZG2=dist_LD.r[c(2,5)],LD.rG1G2=distG_LD.r[c(2,5)])
  
  for( i in 1:nrow(maf_ld_comb) ){
    pG1 <- 1-maf_ld_comb$P_g1[i]
    pG2 <- 1-maf_ld_comb$P_g2[i]
    LDi <- rbind(cbind(1,2, maf_ld_comb$LD.rZG1[i]), 
                 cbind(1,3, maf_ld_comb$LD.rZG2[i]), 
                 cbind(2,3, maf_ld_comb$LD.rG1G2[i]) )
    p_gz0 <- multilocus.Pgz(3, p=c(1-mean(Z)/2, pG1, pG2), LD = LDi)
    
    if( !is.null(p_gz0) ) ind_comb[i] <- i
  }
}
maf_ld_comb <- maf_ld_comb[!is.na(ind_comb),]

ptm <- Sys.time()
### Cluster prep. ----
cl <- makeCluster(ncores)
registerDoParallel(cl)
clusterSetRNGStream(cl, iseed=148399)

### Step 1 - Select phase 2 samples ----
source('Step_1.R') ## Produces dataframes dat_ids_all and dat_alloc which contain the ids for phase 2 sequencing per replicate and corresponding summaries

### Step 2 - Single-variant analysis ----
source('Step_2.R') ## Produces dataframe resOpt_single, which contains the results for the single-SNP analyses across replicates

### Step 3 - Conditional analysis ----
source('Step_3.R') ## Produces dataframe resOpt_conditional, which contains the results for the conditional analyses across replicates
stopCluster(cl)

### Step 4 - Generate tables and plots ----
### Note that in order to recreate the all the figures in the paper, all optimality criteria and phase 2 sample sizes need to be previously analyzed and corresponding datafreames saved. Here only one combination is displayed.
source('Step_4.R')

tim1 <- difftime(Sys.time(),ptm,units = "hours")
cat("Elapsed time: ", tim1, " hours\n",sep = "")
