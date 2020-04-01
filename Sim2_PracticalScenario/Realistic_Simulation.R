
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
# install.packages("twoPhaseGAS_1.05.tar.gz", repos=NULL, type="source")
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

optM <- 1 # Can vary from 1-3 for other optimality criterion (Note that these have no effect in this setting given that the beta_des is considered under the null only)
if(optM==1){
  Kind <- 2
} else Kind <- NULL
optMethod <- c("Par-spec","A-opt","D-opt")[optM]
optMetout <- c("Param","Aopt","Dopt")[optM]

sampl <- c("tagSNP","trait","trait-tagSNP")[3] ## either marginal or joint sampling (only joint sampling, i.e."trait-tagSNP", was used in paper)
sfrac <- 0.5 ## determines phase 2 sample size
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
if(0) source('Step_optional_Data_generation.R') ## Produces a new set of genotypes and saves it into "./data_Realistic_R=1K_N=5K.RData". Generating new genotypes will not reproduce the results from the paper.


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

### check if design LD and MAF combinations are plausible:
maf_ld_comb <- expand.grid(P_g=dist_P_g[c(2,3,5)],LD.r=dist_LD.r[c(2,3,5)])
ind_comb <- rep(NA,nrow(maf_ld_comb))
for( i in 1:nrow(maf_ld_comb) ){
  p_gz0 <- tryCatch(p_gz_func(p_Z, maf_ld_comb$P_g[i], maf_ld_comb$LD.r[i], "G", "Z"), error=function(e){print(e); return(tryCatch(p_gz_func(mean(Z)/2, maf_ld_comb$P_g[i], maf_ld_comb$LD.r[i], "G", "Z"),error=function(e){NULL})) })
  if( !is.null(p_gz0) ) ind_comb[i] <- i
}
maf_ld_comb <- maf_ld_comb[!is.na(ind_comb),]

ptm <- Sys.time()
### Cluster prep. ----
cl <- makeCluster(ncores)
registerDoParallel(cl)
clusterSetRNGStream(cl, iseed=148399)

### Step 1 - Select phase 2 samples ----
source('Step_1.R') ## Produces dataframe dat_ids_all, which contains the ids for phase 2 sequencing per replicate

### Step 2 - Single-SNP analysis ----
source('Step_2.R') ## Produces dataframe resOpt_single, which contains the results for the single-SNP analyses across replicates

### Step 3 - Conditional analysis ----
source('Step_3.R') ## Produces dataframe resOpt_conditional, which contains the results for the conditional analyses across replicates
stopCluster(cl)

### Step 4 - Generate tables and plots ----
source('Step_4.R')

tim1 <- difftime(Sys.time(),ptm,units = "hours")
cat("Elapsed time: ", tim1, " hours\n",sep = "")
