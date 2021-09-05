### This is an example on how to obtain the GA design described in Espin-Garcia, Craiu, and Bull (2021) under Poisson regression.
### It can be run by simply sourcing it, i.e. source('example_Poisson.R'), once the following two packages are installed.

# install.packages("../twoPhaseGAS_1.07.tar.gz", type="source", repos = NULL)
library(twoPhaseGAS)
library(statmod)

### Simulating data ----
Beta0 <- 1 # intercept
Beta1 <- 0.1 # effect size of the sequencing variant
LD.r <- 0.75 # LD between the sequencing variant and the GWAS variant
P.g <- 0.2 # minor allele frequency of the sequencing variant
P.z <- 0.3 # minor allele frequency of the  GWAS variant
N <- 5000 # phase 1  sample size
n2 <- 2000 # phase 2 sample size

### simulate data under Poisson distribution
DataGeneration_TPD <- function (Beta0, Beta1, N, LD.r = 0.75, P_g = 0.2, P_z = 0.3){
  cat("Beta0 =", Beta0, " Beta1 =", Beta1, "\n")
  P_G <- 1 - P_g
  P_Z <- 1 - P_z
  Rsq <- LD.r^2
  LD <- LD.r * sqrt(P_G * P_g * P_Z * P_z)
  Dpr <- LD/((LD < 0) * max(-P_G * P_Z, -P_g * P_z) + (LD >= 
                                                         0) * min(P_G * P_z, P_g * P_Z))
  h.freqs <- rep(0, 4)
  h.freqs[1] <- LD + (1 - P_g) * (1 - P_z)
  h.freqs[2] <- 1 - P_g - h.freqs[1]
  h.freqs[3] <- 1 - P_z - h.freqs[1]
  h.freqs[4] <- P_g - h.freqs[3]
  if( any(h.freqs<0) ){
    h.freqs <- pmax(h.freqs, 1e-05)
    h.freqs <- h.freqs/sum(h.freqs)
  } 
  pval <- 1
  k <- 0
  while (pval >= 1e-05) {
    h1 <- sample(1:4, size = N, prob = h.freqs, replace = TRUE)
    h2 <- sample(1:4, size = N, prob = h.freqs, replace = TRUE)
    G <- as.numeric(h1 == 3 | h1 == 4) + as.numeric(h2 == 
                                                      3 | h2 == 4)
    Z <- as.numeric(h1 == 2 | h1 == 4) + as.numeric(h2 == 
                                                      2 | h2 == 4)
    
    Y <- rpois(N, lambda = exp(Beta0 + G * Beta1))
    lmfit.comp <- glm(Y ~ Z, family = poisson)
    pval <- pchisq(as.numeric(coef(lmfit.comp)[2]^2/diag(vcov(lmfit.comp))[2]), 1, lower.tail = FALSE)
    k <- k + 1
    if (pval < 1e-05) {
      G0 <- sample(0:2, size = N, prob = c(P_G^2, 2 * P_G * 
                                             P_g, P_g^2), replace = TRUE)
      dat_sim_it <- data.frame(wait_it = k, Y = Y, G1 = G, 
                               Z = Z, G0 = G0)
    }
  }
  return(dat_sim_it)
}

set.seed(1279)
dat_sim <- DataGeneration_TPD(N = N, Beta0 = Beta0, Beta1=Beta1, LD.r=LD.r, P_g=P.g, P_z=P.z)
# summary(glm(Y~G1+factor(Z), data=dat_sim, family=poisson))

### Data prep for GA ----
dat_sim$fZ <- factor(dat_sim$Z)
dat = dat_sim[ , c("Y","Z","fZ")]
fit = glm(Y~fZ, data=dat, family = poisson)
disp0 = summary(fit)$dispersion 
pZ <- 1-mean(dat$Z)/2 ##
### summary(fit)

optM <- 1 # One can vary this value from 1-3 for each optimality criterion. When optM==1, beta1 is a scalar and a two-dimensional vector when is not, so the values below need to be adapted accordingly.
if( optM==1 ){
  Kind <- 1
  formula_Ha <- Y~G+fZ
  miscov_formula <- ~G
  beta0 = c(fit$coef[1],0,fit$coef[-1])
  pG <- 1-P.g
  LDj <- cbind(1,2, LD.r)
  p_gz0 <- suppressMessages(multilocus.Pgz(2, p=c(pZ, pG), LD = LDj))
  
} else{
  Kind <- NULL
  formula_Ha <- Y~G0+G1+fZ
  miscov_formula <- ~G0+G1
  beta0 = c(fit$coef[1],0,0,fit$coef[-1])
  pG0 <- 1-P.g
  pG1 <- 1-P.g
  LDj <- rbind(cbind(1,2, LD.r), 
               cbind(1,3, LD.r), 
               cbind(2,3, 0.1) )
  p_gz0 <- suppressMessages(multilocus.Pgz(3, p=c(pZ, pG0, pG1), LD = LDj))
  colnames(p_gz0)[2:3] <- c("G0","G1")
} 
### The formula under the null is the same regardless
formula_Ho <- Y~fZ
auxvar = ~Z
strataformula = ~Z ## Note this will perform covariate-dependent sampling for LM
fam = poisson
optMethod <- c("Par-spec","A-opt","D-opt")[optM]
  
### Obtain GA design ----

#### Inform GA with a good starting population, e.g. by using LM
IM <- obsIM.R(formula=formula_Ha, 
              miscov=miscov_formula, auxvar=~Z,
              family=poisson, dat,
              beta=beta0,p_gz=p_gz0,disp=disp0)

opt.prop <- tryCatch(optimTP.LM(formula=formula_Ha,
                                miscov=miscov_formula,
                                auxvar=~Z,
                                strata=strataformula,
                                family=poisson,
                                n=n2,dat,
                                beta=beta0,p_gz=p_gz0,
                                disp=disp0,
                                optimMeasure=optMethod,K.idx=Kind),error=function(e){print(e); cat("Error in optimTP.LM on iteration:",it,", using same sampling fraction for all strata instead.\n");
  R3.err <<- 1;
  stratadf <- data.frame(xtabs(strataformula,data=data))
  stratadf$prR_cond_optim <- n2/N; return(stratadf)})

pop1a <- sapply(1:120,function(x){
  sa <- which(BSS(samp.fracs=opt.prop$prR_cond_optim,n2,dat,strataformula)$R==1)
  len <- length(sa)
  if( len==n2 ){ return(sa)
  }else if( len>n2 ){ return( sa[-sample(len, len-n2)] )
  }else return( c(sa,sample((1:N)[-sa],n2-len)) )
})
fit1 <- apply(pop1a,2,function(r){
  return(fitnessTP(IM,1*(1:N %in% r),optimMeasure=optMethod,K.idx=Kind))
})

pop2a <- sapply(1:120,function(x){
  sa <- sample(N,n2)
})
fit2 <- apply(pop2a,2,function(r){
  return(fitnessTP(IM,1*(1:N %in% r),optimMeasure=optMethod,K.idx=Kind))
})

### Put initial population (mix between LM and best candidates from SRS)
pop1 <- pop1a[,order(fit1)[1:30]]
pop2 <- pop2a[,order(fit2)[1:30]]

### It may take a few minutes to run
GA.sol1 <- optimTP.GA(ncores=1,
                     formula=formula_Ha,
                     miscov=miscov_formula,
                     auxvar=~Z,
                     family=fam,
                     n=n2,
                     dat,
                     beta=beta0,p_gz=p_gz0,disp=disp0,
                     ga.popsize=60,
                     ga.propelit=0.9,
                     ga.proptourney=0.9,
                     ga.ngen=200,
                     ga.mutrate=0.001,
                     ga.initpop=t(cbind(pop1,pop2)),
                     optimMeasure=optMethod,K.idx=Kind,seed=1)

### Obtain the indicators for simple random sampling (SRS) -best fitness from the ones drawn before- and GA
R1 <- 1*(1:N %in% pop2a[,which(fit2 == min(fit2) )])
R2 <- rep(0,N); R2[GA.sol1$bestsol] <- 1

### Perform analysis ----
## Here, we are analyzing the data in four ways: 1) Complete data, 2) ML with SRS, 3) ML with GA, and 4) Naive (i.e. phase 2 data alone selected using SRS).

res <- data.frame()
for (gcol in c("G0","G1")){ # gcol <- "G0"; gcol <- "G1"
  names(dat_sim)[names(dat_sim)==gcol] <- "G" 
  start.Ho <- list(betas=c(Beta0,0,0))
  if( gcol=="G1" ){
    start.Ha <- list(betas=c(Beta0,Beta1,0,0))
  }else start.Ha <- list(betas=c(Beta0,0,0,0))
    
## Complete data analysis 
lmfit.com <- glm(Y~G+fZ, data=dat_sim, family=poisson)
lmfit.com.Ho <- glm(Y~fZ, data=dat_sim,family=poisson)
### summary(lmfit.com)

### extract test info for complete data
out_com0 <- c(coef(lmfit.com)[2], diag(vcov(lmfit.com))[2])
Wcom <- out_com0[1]^2/out_com0[2]
LRcom <- as.numeric(2*(logLik(lmfit.com) - logLik(lmfit.com.Ho)))
Scom <- glm.scoretest(lmfit.com.Ho,dat_sim[,"G"])^2

out_com <- c( out_com0, Wcom, Scom, LRcom)

names(out_com)<-c("beta1","var_beta1","W","S","LR")

for( i in c(1:2)){ # i=1
  cond = eval(parse(text=paste0("R",i,"==1")))
  data1 <- dat_sim[cond,c("Y","G","Z","fZ")] # phase 2 data
  data0 <- dat_sim[!cond,c("Y","Z","fZ")] # phase 2 data complement
  
  resHo <- twoPhaseSPML(formula=Y~fZ,miscov=~G,auxvar=~Z,family=poisson,data0,data1,start.values=start.Ho,verbose=FALSE)
  resHa <- twoPhaseSPML(formula=Y~G+fZ,miscov=~G,auxvar=~Z,family=poisson,data0,data1,start.values=start.Ha,verbose=FALSE)
  
  assign(paste0("res",i), c( resHa$theta[2], resHa$var_theta[2], resHa$Wobs, resHo$Sobs, 2*(resHa$ll-resHo$ll) ) )
}

## Naive analysis 
dat_nai <- dat_sim[R1==1,c("Y","G","Z","fZ")] # phase 2 data alone
lmfit.nai <- glm(Y~G+fZ, data=dat_nai, family=poisson)
lmfit.nai.Ho <- glm(Y~fZ, data=dat_nai,family=poisson)
### summary(lmfit.nai)

### extract test info for complete data
out_nai0 <- c(coef(lmfit.nai)[2], diag(vcov(lmfit.nai))[2])
Wnai <- out_nai0[1]^2/out_nai0[2]
LRnai <- as.numeric(2*(logLik(lmfit.nai) - logLik(lmfit.nai.Ho)))
Snai <- glm.scoretest(lmfit.nai.Ho,dat_nai[,"G"])^2

out_nai <- c( out_nai0, Wnai, Snai, LRnai)

res <- rbind(res,data.frame(Beta1=Beta1,snp=gcol,Analysis=c("Complete","ML-SRS","ML-GA","Naive"),rbind(out_com,res1,res2,out_nai)))

names(dat_sim)[names(dat_sim)=="G"] <- gcol # put back to the original value
}    

print(res)
### Expected output:
#          Beta1 snp Analysis     beta1       var_beta1          W            S           LR
# out_com    0.1  G0 Complete -0.0099435819 0.0002213691  0.446651419  0.446655131  0.447548135
# res1       0.1  G0   ML-SRS  0.0011314560 0.0005107608  0.002506443  0.002510957  0.002511346
# res2       0.1  G0    ML-GA -0.0007389052 0.0004921429  0.001109395  0.001108488  0.001107784
# out_nai    0.1  G0    Naive  0.0011862183 0.0005349665  0.002630284  0.002630284  0.002629682
# out_com1   0.1  G1 Complete  0.1051296161 0.0004948109  22.336283169 22.353003255 22.573748809
# res11      0.1  G1   ML-SRS  0.0678366207 0.0011091880  4.148807181  4.132324066  4.158262909
# res21      0.1  G1    ML-GA  0.1170936479 0.0005939658  23.083689002 23.050909698 23.322514207
# out_nai1   0.1  G1    Naive  0.0742660816 0.0012214881  4.515353653  4.516995817  4.546878257

### In this example, ML under GA reaches almost the same power as the complete data analysis and is orders of magnitude more significant than both ML under SRS and the Naive analyses.
