
### First, a region of interest needs to be extracted, I'm using data from 1KG.
### This is a vcftools command to extract the region of interest:
# file=ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
# panel=integrated_call_samples_v3.20130502.ALL.panel
# 
# vcftools --gzvcf $file --keep <(grep EUR $panel | cut -f1) --chr 16 --from-bp 56964792 --to-bp 57019949 --recode --stdout | gzip -c > EUR.chr16.HERPUD1-CETP.genotypes.vcf.gz


### genetic map in GRCh37 (from ftp://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37/) or rather https://github.com/adimitromanolakis/geneticMap-GRCh37/ (downloaded by default)

vcfpath <- "."
vcf_file = file.path(vcfpath, "EUR.chr16.HERPUD1-CETP.genotypes.vcf.gz")

### Note that the functions below are non-deterministic

vcf = readVCF( vcf_file , maxNumberOfVariants = 43 , min_maf = 0.05 , max_maf = NA, region_start = 56990716-5000, region_end = 56990716 + 5000 ) 
readGeneticMap( chromosome = 16, dir=vcfpath ) 

# c("56993324","56995236","56989830","56994990","56990716") %in% vcf$vcf$POS

set.seed(123)
startSimulation(vcf, totalNumberOfIndividuals = N)
biminfo <- vcf$vcf
ids = generateUnrelatedIndividuals(N)
genotype = retrieveGenotypes(ids)

# Compare MAF of simulataed data and vcf
# plot( apply(genotype,2,mean)/2 ,  apply(vcf$gt1+vcf$gt2,1,mean)/2 )
# abline(0,1,lty=1,lwd=9,col=rgb(0,0,1,0.3))

### code all SNPs as number of copies of the minor allele
genotype <- apply(genotype,2,function(x){ if((mean(x)/2)>0.5){return(2-x)}else return(x) })
# dim(genotype)

# gplots::heatmap.2(genotype,col=c("white","orange","red"),Colv=F, trace="none")
# gplots::heatmap.2(cor(genotype)^2, col=rev(heat.colors(10)) , trace="none",Rowv=F,Colv=F)

cor.lowertri <- (cor(genotype)^2)[lower.tri(cor(genotype), diag = FALSE)]
# summary(cor.lowertri, na.rm=TRUE)
# summary((cor(genotype))[lower.tri(cor(genotype), diag = FALSE)], na.rm=TRUE)
# head(genotype)
colnames(genotype) <- vcf$vcf$POS
# 56990716
cor.mat <- round(cor(genotype),3)
rownames(cor.mat) <- colnames(cor.mat) <- vcf$vcf$POS

### summaries
# rbind(cor.mat[rownames(cor.mat)=="56990716",]^2, ### LD (r^2)
#       cor.mat[rownames(cor.mat)=="56990716",], 
#       colMeans(genotype)/2) ## MAFs
# genotype[,which(colnames(cor.mat)%in%c(56990716,56989590))]

Gsnps <- c("56993324","56995236","56989830","56994990")

rbind(cor.mat[rownames(cor.mat)=="56990716",]^2, ### LD (r^2)
      cor.mat[rownames(cor.mat)=="56990716",], 
      colMeans(genotype)/2)[,colnames(cor.mat)%in%Gsnps]

tao <- 2/5
Sigma2 <- 1
Beta <- c(2, 0.125,-0.15,-0.20,0.25)

### if regenerating data for this, check that the stats are rougly like this
## G SNPs are 56993324, 56995236, 56989830, 56994990 (which roughly have MAFs of 0.294000, 0.4858, 0.478100, 0.126200 and r(r^2) with Z of 0.942(0.887364), 0.54(0.2916), -0.414(0.171396), -0.147(0.021609)  ). The rest LD structure is:

# cor.mat <- round(cor(genotype),3)
# cor.mat[Gsnps,Gsnps]^2
# cor.mat[c("56990716",Gsnps),c("56990716",Gsnps)]
## Z SNP is 56990716 (MAF=0.2926)
X <- cbind(1, genotype[,Gsnps])
Z <- genotype[,"56990716"]

library(iterators)
library(foreach)
library(rlecuyer)
library(doParallel)

cl <- makeCluster(1)
registerDoParallel(cl)
clusterSetRNGStream(cl, iseed=148399)

dat_sim_all <- foreach(it=icount(1000), .combine='rbind', .inorder=F, .noexport=NULL, .verbose=F) %dopar% {
  pval<-1
  k<-0
  while( pval>=1e-05 ){
    #Simulate the data with Gs, the seq SNPs, and Z the GWAS SNP
    Y <- X %*% Beta + rnorm(N,sd=sqrt(Sigma2))
    lmfit.comp <- lm(Y~Z)   
    pval <- 1-pchisq(as.numeric(coef(lmfit.comp)[2]^2/diag(vcov(lmfit.comp))[2]),1)
    k <- k+1
    if( pval<1e-05 ){       
      Yst <- sapply(Y,function(x){
        if( x<qnorm(tao,mean=Beta[1],sd=sqrt(Sigma2)) ){return(1)
        }else if( x>qnorm(1-tao,mean=Beta[1],sd=sqrt(Sigma2)) ){return(3)
        }else return(2)})
      dat_sim <- data.frame(wait_it=k,Y=Y,S=Yst)
    }
  }
  return(cbind(iteration=it,dat_sim))
}

stopCluster(cl)

save(dat_sim_all,genotype,biminfo,file="data_Realistic_R=1K_N=5K.RData")