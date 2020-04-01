# Packages and functions ----
library(reshape)
library(ggplot2)
library(vcd)
library(gridExtra)
library(RColorBrewer) 
library(data.table)

### Prepare auxiliary information ----

### Set up the true values for Beta (as designed in Step_optional_Data_generation.R)
Betas_df <- data.frame(Gpos=c("56993324","56995236","56989830","56994990"),Beta1=c(0.125,-0.15,-0.20,0.25))

formt <- function(x,ndec){format(round(x, ndec), nsmall = ndec)} # auxiliary function

### SNPs classification
Zsnp <- "56990716"
Gsnps <- c("56993324","56995236","56989830","56994990" )
Hitchsnps <- c("56995038","56994192","56993935","56990803","56988958","56986762","56986045", "56991143")

MAFs <- colMeans(genotype)/2
mean_LD.r <- mean(cor.mat[rownames(cor.mat)==Zpos,-which(colnames(genotype)==Zpos)])

Dprime <- function(LD.r,P_z,P_g){
  P_G <- 1 - P_g; P_Z <- 1 - P_z ;
  Dgz <- LD.r * sqrt(P_G*P_g*P_Z*P_z)
  Dpr <- Dgz/((Dgz<0)*min(P_G*P_Z,P_g*P_z)+(Dgz>=0)*min(P_G*P_z,P_g*P_Z)) # Lies between -1 and 1
  return(Dpr)
}
#### pairwise D'
p <- length(MAFs)
Dprime.mat <- matrix(NA,ncol=p,nrow=p)
for(i in 1:p){
  for(j in 1:p){
    Dprime.mat[i,j] <- Dprime(cor.mat[i,j],MAFs[i],MAFs[j])
  }
}
colnames(Dprime.mat) <- rownames(Dprime.mat) <- colnames(cor.mat)

reg_sum <- data.frame(Gpos=colnames(genotype_noZ),MAF=colMeans(genotype_noZ)/2,r=cor.mat[c(Zsnp),colnames(genotype_noZ)],Dpr=Dprime.mat[c(Zsnp),colnames(genotype_noZ)])
reg_sum$r2 <- reg_sum$r^2
reg_sum$ID0 <- substr(reg_sum$Gpos,4,nchar(as.character(reg_sum$Gpos)))
reg_sum$ID <- sapply(reg_sum$ID0,function(s){
  if(s %in% substr(Gsnps,4,nchar(Gsnps))) return(paste0(s,"$^*$"))
  if(s %in% substr(Hitchsnps,4,nchar(Hitchsnps))) return(paste0(s,"$^\\dagger$"))
  return(s)
})

#* Causal SNPs characteristics ----
{   
  subset(reg_sum,grepl("[*]",ID))[,c("ID0","MAF","r","r2","Dpr")]
} 

### Allocation plots ----
#* Data prep ----
{
  ### dat_alloc is one output of Step_1.R 
  indx_tot_sim2 <- subset(dat_alloc, error==0 | is.na(error))
  indx_tot_sim2$N <- N
  indx_tot_sim2$OptCriterion <- optMethod
  indx_tot1_sim2_mosaic <- melt(indx_tot_sim2,id.vars=c("N","samp","ss","alloc","it", "OptCriterion"),measure.vars=c("Z0.S1", "Z1.S1", "Z2.S1", "Z0.S2","Z1.S2", "Z2.S2", "Z0.S3", "Z1.S3", "Z2.S3"))

  strata <- strsplit(as.character(indx_tot1_sim2_mosaic$variable),split = "[.]")
  indx_tot1_sim2_mosaic$Z <- gsub("Z","",sapply(strata,`[[`,1))
  indx_tot1_sim2_mosaic$Yst <- sapply(strata,`[[`,2)

  indx_tot1_sim2_mosaic$alloc <- factor(indx_tot1_sim2_mosaic$alloc,levels = c("Opt.Lagr", "Opt.GA", "Combined", "RDS", "TZL","Complete"))
  indx_tot1_sim2_mosaic$alloc <- factor(indx_tot1_sim2_mosaic$alloc)
}
#* Mosaic plots ----
{
  PaletteAlloc <- grey.colors(3, start = 0.4, end = 0.8, gamma = 2.2, alpha = NULL)
  
  vnames <- list(set_varnames = c(ss="",OptCriterion="",alloc=""))
  vnames1 <- list(set_varnames = c(ss="",OptCriterion="",alloc="", Z=""))
  lnames <- list(ss = c("540","1.25K","2.5K"), Yst=c(expression(T["1"]),expression(T["2"]),expression(T["3"])), OptCriterion = c("A-opt","D-opt","Par-spec"),alloc=c("LM","GA","Comb","RDS","TZL","Comp"))
  lnames1 <- list(ss = c("540","1.25K","2.5K"), Yst=c(expression(T["1"]),expression(T["2"]),expression(T["3"])), OptCriterion = c("A-opt","D-opt","Par-spec"),alloc=rep("",5))
  
  n1 <- grid.grabExpr(mosaic(xtabs(value ~ alloc + Yst + Z, droplevels(subset(indx_tot1_sim2_mosaic,ss==540 & alloc!="Complete"))), zero_size = 0, main="540", labeling_args=vnames1, set_labels=lnames, keep_aspect_ratio=FALSE, margins=unit(c(2.5,2.5,0.5,2.5), "lines"), highlighting=3, highlighting_fill = PaletteAlloc, main_gp = gpar(fontsize = 14)))
  
  n2 <- grid.grabExpr(mosaic(xtabs(value ~ alloc + Yst + Z, droplevels(subset(indx_tot1_sim2_mosaic,ss==1250 & alloc!="Complete"))), zero_size = 0, main="1250", labeling_args=vnames1, set_labels=lnames1, keep_aspect_ratio=FALSE, margins=unit(c(2.5,2.5,0.5,2.5), "lines"), highlighting=3, highlighting_fill = PaletteAlloc, main_gp = gpar(fontsize = 14)))
  
  n3 <- grid.grabExpr(mosaic(xtabs(value ~ alloc + Yst + Z, droplevels(subset(indx_tot1_sim2_mosaic,ss==2500 & alloc!="Complete"))), zero_size = 0, main="2500", labeling_args=vnames, set_labels=lnames1, keep_aspect_ratio=FALSE, margins=unit(c(2.5,2.5,0.5,2.5), "lines"), highlighting=3, highlighting_fill = PaletteAlloc, main_gp = gpar(fontsize = 14)))
  
  
  obj <- grid.arrange(n1, n2, n3, ncol=3)
  if( save_ind ) png(paste0(savedir,"/Mosaics_Sim2.png"),width = 7, height = 4, units = "in", bg="transparent", res=300)
  grid.draw(obj)
  if( save_ind ) dev.off()
  
}
### Single-SNP analysis ----
#* Data prep ----
{
  resOpt_sim2 <- resOpt_single
  resOpt_sim2$N <- N
  resOpt_sim2$OptCriterion <- optMethod
  resOpt_sim2$Alloc <- relevel(resOpt_sim2$Alloc,ref="Complete")
  resOpt_sim2 <- merge(resOpt_single,Betas_df,by="Gpos",all.x=TRUE)
  resOpt_sim2$Beta1[is.na(resOpt_sim2$Beta1)] <- 0
  resOpt_sim2$Alloc <- factor(resOpt_sim2$Alloc, levels = c("Complete","Opt.Lagr","Opt.GA","Combined","RDS","TZL"))
  resOpt_sim2$Alloc_expr <- resOpt_sim2$Alloc
  levels(resOpt_sim2$Alloc_expr) <- expression("bold(Complete)","bold(LM)","bold(GA)","bold(Comb)","bold(RDS)","bold(TZL)")
}
#* Summarize across replicates ----
{
  sig_Ha<-0.05/29; sig_Ho<-0.05/29
  tab_sim2 <- do.call( rbind,by(resOpt_sim2, INDICES = list(resOpt_sim2$samp,resOpt_sim2$Alloc,resOpt_sim2$Gpos,resOpt_sim2$ss,resOpt_sim2$Beta1), function(x){
    data.frame(ss=x$ss[1],Gpos=x$Gpos[1],Beta1=x$Beta1[1],Alloc=x$Alloc[1],samp=x$samp[1],
               Mean_beta1=mean(x$beta1,na.rm=T),Median_beta1=median(x$beta1,na.rm=T),
               SD_beta1=sd(x$beta1,na.rm=T),Mean_SEbeta1=mean(sqrt(x$var_beta1),na.rm=T),
               Bias_beta1=mean(x$beta1-x$Beta1,na.rm=T),MAB_beta1=mean(abs(x$beta1-x$Beta1),na.rm=T),
               RMSE_beta1=sqrt(mean((x$beta1-x$Beta1)^2,na.rm=T)),
               CP_beta1=sum(x$Beta1>x$beta1-qnorm(0.975)*sqrt(x$var_beta1) & x$Beta1<x$beta1+qnorm(0.975)*sqrt(x$var_beta1),na.rm=T)/sum(!is.na(sqrt(x$var_beta1))), 
               t1e_power_Wald=sum( ((x$Beta1==0)*(pchisq(x$W,1,lower.tail = FALSE)<sig_Ho)+(x$Beta1!=0)*(pchisq(x$W,1,lower.tail = FALSE)<sig_Ha)),na.rm=T )/ sum(!is.na(pchisq(x$W,1,lower.tail = FALSE))),
               t1e_power_LRS=sum( ((x$Beta1==0)*(pchisq(x$LR,1,lower.tail = FALSE)<sig_Ho)+(x$Beta1!=0)*(pchisq(x$LR,1,lower.tail = FALSE)<sig_Ha)),na.rm=T )/ sum(!is.na(pchisq(x$LR,1,lower.tail = FALSE))),
               t1e_power_Score=sum( ((x$Beta1==0)*(pchisq(x$S,1,lower.tail = FALSE)<sig_Ho)+(x$Beta1!=0)*(pchisq(x$S,1,lower.tail = FALSE)<sig_Ha)),na.rm=T )/ sum(!is.na(pchisq(x$S,1,lower.tail = FALSE))),
               N_Wald = sum(!is.na(pchisq(x$W,1,lower.tail = FALSE))),
               N_LRS = sum(!is.na(pchisq(x$LR,1,lower.tail = FALSE))),
               N_Score = sum(!is.na(pchisq(x$S,1,lower.tail = FALSE)))
    )
  } ) )
  
  Tpow_sim2 <- rbind(cbind(Test="Wald",cast(tab_sim2,ss+Gpos+Beta1+Alloc~.,value="t1e_power_Wald")),
                   cbind(Test="LR",cast(tab_sim2,ss+Gpos+Beta1+Alloc~.,value="t1e_power_LRS")),
                   cbind(Test="Score",cast(tab_sim2,ss+Gpos+Beta1+Alloc~.,value="t1e_power_Score")))
  names(Tpow_sim2)[which(names(Tpow_sim2)=="(all)")] <- "Power"

  TpowScore_sim2 <- subset(Tpow_sim2, Test=="Score") 
}
#* Estimation and Power ----
{
  summreps <- tab_sim2
  summreps$Sum <- paste0(formt(summreps$Mean_beta1,3)," (",formt(summreps$Mean_SEbeta1,3),")") 
  summreps$Pow <- formt(summreps$t1e_power_Score*100,1)
  #summreps$Pow <- formt(summreps$t1e_power_Wald*100,1)
  #summreps$Pow <- formt(summreps$t1e_power_LRS*100,1)
  summreps <- merge(reg_sum,summreps,by="Gpos")
  summreps$MAF <- formt(summreps$MAF*100,1)
  summreps <- summreps[order(as.numeric(as.character(summreps$Gpos))),]
  
  cast(subset(summreps,ss==540),ID+MAF+r+r2+Dpr+Beta1~Alloc,value="Sum")
  cast(subset(summreps,ss==540),ID+Beta1~Alloc,value="Pow")
  
  ## sumarize other sample sizes but only causal variants:
  cast(subset(summreps,ss!=540 & grepl("[*]",ID)),ss+ID0+MAF+r+r2+Dpr+Beta1~Alloc,value="Sum")
  cast(subset(summreps,ss!=540 & grepl("[*]",ID)),ss+ID0+Beta1~Alloc,value="Pow")
  
}

### Conditional analysis ----
#* Data prep ----
{
  resOpt_sim2_cond <- resOpt_conditional
  resOpt_sim2_cond$N <- N
  resOpt_sim2_cond$OptCriterion <- optMethod
  resOpt_sim2_cond$Alloc <- relevel(resOpt_sim2_cond$Alloc,ref="Complete")
  
  resOpt_sim2_cond <- merge(resOpt_sim2_cond,Betas_df,by="Gpos",all.x=TRUE)
  resOpt_sim2_cond <- merge(resOpt_sim2_cond,Betas_df,by.x="Gcond",by.y="Gpos",all.x=TRUE)
  names(resOpt_sim2_cond)[which(names(resOpt_sim2_cond)%in%c("Beta1.x","Beta1.y"))] <- c("Beta1","Beta2_cond")
  resOpt_sim2_cond$Beta1[is.na(resOpt_sim2_cond$Beta1)] <- 0
  resOpt_sim2_cond$Beta2_cond[is.na(resOpt_sim2_cond$Beta2_cond)] <- 0
  resOpt_sim2_cond$Alloc <- factor(resOpt_sim2$Alloc, levels = c("Complete","Opt.Lagr","Opt.GA","Combined","RDS","TZL"))
  resOpt_sim2_cond$Alloc_expr <- resOpt_sim2_cond$Alloc
  levels(resOpt_sim2_cond$Alloc_expr) <- expression("bold(Complete)","bold(LM)","bold(GA)","bold(Comb)","bold(RDS)","bold(TZL)")
  levels(resOpt_sim2_cond$OptCriterion) <- c("Par-spec")
}
#* Summarize across replicates----
{
  #** Beta1 ----
  sig_Ha<-0.05/29; sig_Ho<-0.05/29
  tab_sim2cond <- do.call( rbind,by( resOpt_sim2_cond, INDICES = list(resOpt_sim2_cond$OptCriterion,resOpt_sim2_cond$samp,resOpt_sim2_cond$Alloc,resOpt_sim2_cond$Gpos,resOpt_sim2_cond$ss,resOpt_sim2_cond$Beta1), function(x){
    data.frame(OptCriterion=x$OptCriterion[1],ss=x$ss[1],Gpos=x$Gpos[1],Beta1=x$Beta1[1],Alloc=x$Alloc[1],samp=x$samp[1],
               Mean_beta1=mean(x$beta1,na.rm=T),Median_beta1=median(x$beta1,na.rm=T),
               SD_beta1=sd(x$beta1,na.rm=T),Mean_SEbeta1=mean(sqrt(x$var_beta1),na.rm=T),
               Bias_beta1=mean(x$beta1-x$Beta1,na.rm=T),MAB_beta1=mean(abs(x$beta1-x$Beta1),na.rm=T),
               RMSE_beta1=sqrt(mean((x$beta1-x$Beta1)^2,na.rm=T)),
               CP_beta1=sum(x$Beta1>x$beta1-qnorm(0.975)*sqrt(x$var_beta1) & x$Beta1<x$beta1+qnorm(0.975)*sqrt(x$var_beta1),na.rm=T)/sum(!is.na(sqrt(x$var_beta1))), 
               
               t1e_power_Wald=sum( ((x$Beta1==0)*(pchisq(x$W_G,1,lower.tail = FALSE)<sig_Ho)+(x$Beta1!=0)*(pchisq(x$W_G,1,lower.tail = FALSE)<sig_Ha)),na.rm=T )/ sum(!is.na(pchisq(x$W_G,1,lower.tail = FALSE))),
               t1e_power_LRS=sum( ((x$Beta1==0)*(pchisq(x$LR_G,1,lower.tail = FALSE)<sig_Ho)+(x$Beta1!=0)*(pchisq(x$LR_G,1,lower.tail = FALSE)<sig_Ha)),na.rm=T )/ sum(!is.na(pchisq(x$LR_G,1,lower.tail = FALSE))),
               t1e_power_Score=sum( ((x$Beta1==0)*(pchisq(x$S_G,1,lower.tail = FALSE)<sig_Ho)+(x$Beta1!=0)*(pchisq(x$S_G,1,lower.tail = FALSE)<sig_Ha)),na.rm=T )/ sum(!is.na(pchisq(x$S_G,1,lower.tail = FALSE))),
               
               N_Wald = sum(!is.na(pchisq(x$W_G,1,lower.tail = FALSE))),
               N_LRS = sum(!is.na(pchisq(x$LR_G,1,lower.tail = FALSE))),
               N_Score = sum(!is.na(pchisq(x$S_G,1,lower.tail = FALSE)))
    )
  } ) )
  # cast(tab_sim2cond,ss+OptCriterion+Beta1+Gpos~Alloc,value="N_Score")
  formula <- Gpos+OptCriterion+Beta1+ss+Alloc~.
  Tpow_sim2cond <- rbind(cbind(Test="Wald",cast(tab_sim2cond,formula,value="t1e_power_Wald")),
                         cbind(Test="LR",cast(tab_sim2cond,formula,value="t1e_power_LRS")),
                         cbind(Test="Score",cast(tab_sim2cond,formula,value="t1e_power_Score")))
  names(Tpow_sim2cond)[which(names(Tpow_sim2cond)=="(all)")] <- "Power"
  
  # subset(Tpow_sim1, snp=="G0" & Alloc=="Complete")
  ## only score test
  TpowScore_sim2cond <- subset(Tpow_sim2cond, Test=="Score")
  
  #** Beta2 ----
  sig_Ha<-0.05/29; sig_Ho<-0.05/29
  tab_sim2cond1 <- do.call( rbind,by( resOpt_sim2_cond,INDICES = list(resOpt_sim2_cond$OptCriterion,resOpt_sim2_cond$samp,resOpt_sim2_cond$Alloc,resOpt_sim2_cond$Gcond,resOpt_sim2_cond$ss,resOpt_sim2_cond$Beta2_cond), function(x){
    data.frame(OptCriterion=x$OptCriterion[1],ss=x$ss[1],Gcond=x$Gcond[1],Beta2_cond=x$Beta2_cond[1],Alloc=x$Alloc[1],samp=x$samp[1],
               Mean_beta2=mean(x$beta2_cond,na.rm=T),Median_beta2=median(x$beta2_cond,na.rm=T),
               SD_beta2=sd(x$beta2_cond,na.rm=T),Mean_SEbeta2=mean(sqrt(x$var_beta2_cond),na.rm=T),
               Bias_beta2=mean(x$beta2_cond-x$Beta2_cond,na.rm=T),MAB_beta2=mean(abs(x$beta2_cond-x$Beta2_cond),na.rm=T),
               RMSE_beta2=sqrt(mean((x$beta2_cond-x$Beta2_cond)^2,na.rm=T)),
               CP_beta2=sum(x$Beta2_cond>x$beta2_cond-qnorm(0.975)*sqrt(x$var_beta2_cond) & x$Beta2_cond<x$beta2_cond+qnorm(0.975)*sqrt(x$var_beta2_cond),na.rm=T)/sum(!is.na(sqrt(x$var_beta2_cond))), 
               
               t1e_power_Wald=sum( ((x$Beta2_cond==0)*(pchisq(x$W_Gcond,1,lower.tail = FALSE)<sig_Ho)+(x$Beta2_cond!=0)*(pchisq(x$W_Gcond,1,lower.tail = FALSE)<sig_Ha)),na.rm=T )/ sum(!is.na(pchisq(x$W_Gcond,1,lower.tail = FALSE))),
               t1e_power_LRS=sum( ((x$Beta2_cond==0)*(pchisq(x$LR_Gcond,1,lower.tail = FALSE)<sig_Ho)+(x$Beta2_cond!=0)*(pchisq(x$LR_Gcond,1,lower.tail = FALSE)<sig_Ha)),na.rm=T )/ sum(!is.na(pchisq(x$LR_Gcond,1,lower.tail = FALSE))),
               t1e_power_Score=sum( ((x$Beta2_cond==0)*(pchisq(x$S_Gcond,1,lower.tail = FALSE)<sig_Ho)+(x$Beta2_cond!=0)*(pchisq(x$S_Gcond,1,lower.tail = FALSE)<sig_Ha)),na.rm=T )/ sum(!is.na(pchisq(x$S_Gcond,1,lower.tail = FALSE))),
               
               N_Wald = sum(!is.na(pchisq(x$W_Gcond,1,lower.tail = FALSE))),
               N_LRS = sum(!is.na(pchisq(x$LR_Gcond,1,lower.tail = FALSE))),
               N_Score = sum(!is.na(pchisq(x$S_Gcond,1,lower.tail = FALSE)))
    )
  } ) )
}
#* Estimation and Power ----
{   

  #** Beta1 ----
  summreps <- subset(tab_sim2cond, OptCriterion=="Par-spec")
  summreps$Sum <- paste0(formt(summreps$Mean_beta1,3)," (",formt(summreps$Mean_SEbeta1,3),")") 
  
  summreps$Pow <- formt(summreps$t1e_power_Score*100,1)
  #summreps$Pow <- formt(summreps$t1e_power_Wald*100,1)
  #summreps$Pow <- formt(summreps$t1e_power_LRS*100,1)
  summreps <- merge(reg_sum,summreps,by="Gpos")
  summreps$MAF <- formt(summreps$MAF*100,1)
  summreps <- summreps[order(as.numeric(as.character(summreps$Gpos))),]
  
  ## Results
  cast(subset(summreps,!Gpos%in%c("56989830") & ss==540),ID+Beta1~Alloc,value="Sum")
  cast(subset(summreps,!Gpos%in%c("56989830") & ss==540),ID+Beta1~Alloc,value="Pow")
  
  ## sumarize other sample sizes but only causal variants:
  cast(subset(summreps,!Gpos%in%c("56989830") & ss!=540 & grepl("[*]",ID)),ss+ID0+Beta1~Alloc,value="Sum")
  cast(subset(summreps,!Gpos%in%c("56989830") & ss!=540 & grepl("[*]",ID)),ss+ID0+Beta1~Alloc,value="Pow")
  
  #** Beta2 ----
  summreps1 <- subset(tab_sim2cond1, OptCriterion=="Par-spec")
  summreps1$Sum <- paste0(formt(summreps1$Mean_beta2,3)," (",formt(summreps1$Mean_SEbeta2,3),")") 
  
  summreps1$Pow <- formt(summreps1$t1e_power_Score*100,1)
  #summreps$Pow <- formt(summreps$t1e_power_Wald*100,1)
  #summreps$Pow <- formt(summreps$t1e_power_LRS*100,1)
  summreps2 <- merge(reg_sum,summreps1,by.x="Gpos",by.y="Gcond")
  summreps2$MAF <- formt(summreps2$MAF*100,1)
  summreps2 <- summreps2[order(as.numeric(as.character(summreps2$Gpos))),]
  
  ## Results
  cast(subset(summreps2,Gpos%in%c("56986762",Gsnps)),ss+ID+Beta2_cond~Alloc,value="Sum")
  cast(subset(summreps2,Gpos%in%c("56986762",Gsnps)),ss+ID+Beta2_cond~Alloc,value="Pow")
  
} 

# Combined boxplot ----
{
  Hitchsnps <- c("56995038","56994192","56993935","56990803","56988958","56986762","56986045", "56991143")
  
  resOpt_sim2$mlog10p <- -log10(pchisq(resOpt_sim2$S,1,lower.tail = FALSE))
  resOpt_sim2$Kbp <- round(as.numeric(as.character(resOpt_sim2$Gpos))/1000,2)
  resOpt_sim2$causal <- ifelse(resOpt_sim2$Beta1==0,"Non-causal","Causal")
  resOpt_sim2$causal[resOpt_sim2$Gpos%in%Hitchsnps] <- "Hitchhiker"
  
  resOpt_sim2_cond$mlog10p <- -log10(pchisq(resOpt_sim2_cond$S_G,1,lower.tail = FALSE))
  resOpt_sim2_cond$Kbp <- round(as.numeric(as.character(resOpt_sim2_cond$Gpos))/1000,2)
  resOpt_sim2_cond$causal <- ifelse(resOpt_sim2_cond$Beta1==0,"Non-causal","Causal")
  resOpt_sim2_cond$causal[resOpt_sim2_cond$Gpos%in%Hitchsnps] <- "Hitchhiker"
  
  cols.sel <- c("ss","mlog10p","Kbp","Alloc_expr","causal")
  
  df <- rbind(cbind(Analysis="Conditional",subset(resOpt_sim2_cond,OptCriterion=="Par-spec")[,cols.sel]),
              cbind(Analysis="Single-variant",subset(resOpt_sim2)[,cols.sel]))
  df$Analysis <- relevel(df$Analysis, ref="Single-variant")
  
  for(ph2ss in c(540,1250,2500)){
    boxplts_combined <- ggplot(subset(df,ss==ph2ss), aes(y=mlog10p,x=Kbp,group=Kbp)) + 
      geom_abline(intercept=-log10(0.05/29), slope=0, colour="red", size=.9, linetype=2) + 
      geom_boxplot(aes(fill=causal),color="grey25",position=position_dodge(width = 1.5),width=.5,outlier.size=0.4,lwd=0.2) + 
      facet_grid(Analysis ~ Alloc_expr, labeller=labeller(Alloc_expr=label_parsed), scales = "free_y") + 
      scale_fill_manual(values=cbPalette) + guides(fill=guide_legend(title=""))  + 
      theme_bw(base_size = 11) + ylab(expression(-log["10"]*(p))) + xlab("Chromosome 16 position (Kbp)") + theme(axis.text.x = element_text(vjust=0.5,angle = 0, size=8), axis.text.y  = element_text(vjust=0.5, size=10), axis.title.x = element_text(face="bold",size=15,vjust=-0.5), axis.title.y = element_text(face="bold",size=15), strip.text.x = element_text(size=12, face="bold",color="#FFFFFF"), strip.text.y = element_text(size=17, face="bold",color="#FFFFFF"), strip.background = element_rect(fill="black"), legend.position="bottom", legend.title=element_text(size=17), legend.text=element_text(size=13),legend.key = element_rect(color="black",fill="transparent"),panel.background=element_blank(), panel.grid.minor=element_blank(), panel.grid.major=element_blank(), plot.background=element_blank(),legend.background=element_blank(),plot.title = element_text(lineheight=.8, face="bold",size=17))
    
    if( save_ind ) ggsave(boxplts_combined, file=paste0(savedir,"/Boxplots_sing-cond_analysis_ss",ph2ss,".png"), width=11, height=7, units="in", bg="transparent", dpi=500)
  }
}







