
formula_Ha <- Y~G+Gcond+fZ
formula_Ho <- Y~fZ
formula_Ho1 <- Y~Gcond+fZ
formula_Ho2 <- Y~G+fZ


### define function to perform multiple df. score test (extending glm.scoretest from statmod)
glm.scoretest.multidf <- function(fit, x2, dispersion = NULL) {
  w <- fit$weights
  r <- fit$residuals
  if (any(w <= 0)) {
    r <- r[w > 0]
    x2 <- x2[w > 0]
    w <- w[w > 0]
  }
  if (is.null(dispersion)) {
    fixed.dispersion <- (fit$family$family %in% c("poisson", 
                                                  "binomial"))
    if (fixed.dispersion) 
      dispersion <- 1
    else if (fit$df.residual > 0) {
      dispersion <- sum(w * r^2)/fit$df.residual
    }
    else {
      stop("No residual df available to estimate dispersion")
    }
  }
  # ws <- sqrt(w)
  # x2.1w <- qr.resid(fit$qr, ws * x2)
  # zw <- ws * r
  # colSums(as.matrix(x2.1w * zw))/sqrt(colSums(as.matrix(x2.1w * x2.1w)))/sqrt(dispersion)
  
  X1 <- cbind(model.matrix(fit$formula,fit$data),x2) 
  
  eta <- fit$linear.predictors ### this only works if the null is for zero
  W <- fit$family$mu.eta(eta)^2/fit$family$variance(fit$family$linkinv(eta))
  
  F1 <- crossprod(qr.R(qr(X1*sqrt(as.numeric(W)))))/dispersion
  
  Ind <- rep(FALSE,ncol(X1))
  Ind[which(colnames(X1)%in%colnames(x2))] <- TRUE
  V1 <- F1[Ind,Ind] - F1[Ind,!Ind] %*% solve(F1[!Ind,!Ind]) %*% F1[!Ind,Ind]
  
  mu.eta.val <- fit$family$mu.eta(eta)
  S <- X1 * W * (fit$y - fit$family$linkinv(eta))/(dispersion*mu.eta.val)
  U <- colSums(S)[Ind]
  return( as.numeric(t(U) %*% solve( V1 ) %*% U) )
}

### identify the SNP with minimum p-value (or equivalently the maximum test statistic) for each iteration and design
max_S <- aggregate(resOpt_single$S,by=list(it=resOpt_single$it,Alloc=resOpt_single$Alloc), FUN=max)
Gmax <- merge(resOpt_single[,c(c("it","Alloc","S","p","Gpos"))], max_S, by.x=c("it","Alloc","S"),by.y=c("it","Alloc","x"))
Gmax <- Gmax[order(Gmax$it),]

namscols <- c("beta1","beta2_cond","var_beta1","var_beta2_cond","W2df","S2df","LR2df", "W_G","S_G","LR_G", "W_Gcond","S_Gcond","LR_Gcond")
res_NULL <- rep(NA,13); names(res_NULL) <- namscols

resOpt_conditional <- foreach(G=iter(genotype_noZ,by='column'), p=icount(), .combine='rbind', .verbose=T) %:%
  foreach(data_iter=isplit( dat_sim_all, list(it=dat_sim_all$iteration) ), ids_iter=isplit( dat_ids_all, list(it=dat_ids_all$it) ), .combine='rbind', .inorder=F, .verbose=T, .packages = c("MASS","statmod","Matrix","data.table","enrichwith","aod","twoPhaseGAS")) %dopar% {
    
    it <- paste(data_iter$key,collapse="_")
    dat_sim0 <- data_iter$value
    
    it_idx <- paste(ids_iter$key,collapse="_")
    indx_it <- ids_iter$value
    
    condGs <- Gmax[Gmax$it==it,]
    condGs <- condGs[match(c("Complete","Combined","RDS","TZL","Opt.Lagr","Opt.GA"),condGs$Alloc),]
    condGs$Gpos <- as.character(condGs$Gpos)
    uniqcondGs <- unique(condGs$Gpos)
    condGs_mat <- genotype_noZ[,uniqcondGs,drop=FALSE]
    
    Gnam <- colnames(G) 
    if ( all(uniqcondGs==Gnam) ) return( data.frame(p=p,Gpos=Gnam,Gcond=condGs$Gpos,samp=sampl,ss=n2,it=as.integer(it),itidx=as.integer(it_idx),Alloc=c("Complete","Combined","RDS","TZL","Opt.Lagr","Opt.GA"),rbind(res_NULL,res_NULL,res_NULL,res_NULL,res_NULL,res_NULL)) )
    
    colnames(G) <- "G"
    dat_sim <- cbind(dat_sim0,Z=Z,G,condGs_mat)
    dat_sim$fZ <- factor(dat_sim$Z)
    
    ## Complete data analysis 
    Gpos.Complete <- condGs$Gpos[condGs$Alloc=="Complete"]
    if( Gnam != Gpos.Complete ){
      names(dat_sim)[which(names(dat_sim)==Gpos.Complete)] <- "Gcond"
      lmfit.com <- glm(formula_Ha, data=dat_sim, family=gaussian) 
      lmfit.com.Ho <- glm(formula_Ho, data=dat_sim, family=gaussian)
      lmfit.com.Ho1 <- glm(formula_Ho1, data=dat_sim, family=gaussian)
      lmfit.com.Ho2 <- glm(formula_Ho2, data=dat_sim, family=gaussian)
      ### 2df test
      Wcom <- unname(wald.test(b=coef(lmfit.com),Sigma = vcov(lmfit.com),Terms = 2:3)$result$chi2[1])
      LRcom <- as.numeric(2*(logLik(lmfit.com) - logLik(lmfit.com.Ho)))
      Scom <- glm.scoretest.multidf(lmfit.com.Ho,as.matrix(dat_sim[,c("G","Gcond")]))
      ### 1df test for G (adjusted for Gcond)
      Wcom1 <- unname(wald.test(b=coef(lmfit.com),Sigma = vcov(lmfit.com),Terms = 2)$result$chi2[1])
      LRcom1 <- as.numeric(2*(logLik(lmfit.com) - logLik(lmfit.com.Ho1)))
      Scom1 <- glm.scoretest.multidf(lmfit.com.Ho1,as.matrix(dat_sim[,c("G"),drop=FALSE]))
      
      ### 1df test for Gcond (adjusted for G)
      Wcom2 <- unname(wald.test(b=coef(lmfit.com),Sigma = vcov(lmfit.com),Terms = 3)$result$chi2[1])
      LRcom2 <- as.numeric(2*(logLik(lmfit.com) - logLik(lmfit.com.Ho2)))
      Scom2 <- glm.scoretest.multidf(lmfit.com.Ho2,as.matrix(dat_sim[,c("Gcond"),drop=FALSE]))
      
      out_com <- c( coef(lmfit.com)[2:3], diag(vcov(lmfit.com))[2:3], Wcom, Scom, LRcom, Wcom1,Scom1,LRcom1, Wcom2,Scom2,LRcom2 )
      names(out_com) <- namscols
      names(dat_sim)[which(names(dat_sim)=="Gcond")] <- Gpos.Complete
    }else out_com <- res_NULL
    
    
    ### Extract ids
    R0 <- unlist(indx_it$R0)
    R1 <- unlist(indx_it$R1)
    R2 <- unlist(indx_it$R2)
    R3 <- unlist(indx_it$R3)
    R4 <- unlist(indx_it$R4)
    
    ### SPML analyses
    for( i in c(0:4)){ # i=1
      des = c("Combined","RDS","TZL","Opt.Lagr","Opt.GA")[i+1]
      Gpos.des <- condGs$Gpos[condGs$Alloc==des]
      
      if( Gnam != Gpos.des ){
        names(dat_sim)[which(names(dat_sim)==Gpos.des)] <- "Gcond"
        cond = eval(parse(text=paste0("R",i,"==1")))
        data1 <- dat_sim[cond,c("Y","G","Gcond","Z","fZ","S")] # phase 2 data
        data0 <- dat_sim[!cond,c("Y","Z","fZ","S")] # phase 2 data
        
        names(dat_sim)[which(names(dat_sim)=="Gcond")] <- Gpos.des
        
        resHo_a <- tryCatch(twoPhaseSPML(formula=formula_Ho,miscov=~G+Gcond,auxvar=~Z,family=gaussian,data0,data1,verbose=FALSE), error=function(e){print(e); list(Sexp=NA, Sobs=NA,ll=NA) })
        resHo1_a <- tryCatch(twoPhaseSPML(formula=formula_Ho1,miscov=~G+Gcond,auxvar=~Z,family=gaussian,data0,data1,verbose=FALSE), error=function(e){print(e); list(Sexp=NA, Sobs=NA,ll=NA) })
        resHo2_a <- tryCatch(twoPhaseSPML(formula=formula_Ho2,miscov=~G+Gcond,auxvar=~Z,family=gaussian,data0,data1,verbose=FALSE), error=function(e){print(e); list(Sexp=NA, Sobs=NA,ll=NA) })
        resHa_a <- tryCatch(twoPhaseSPML(formula=formula_Ha,miscov=~G+Gcond,auxvar=~Z,family=gaussian,data0,data1,verbose=FALSE), error=function(e){print(e); list(theta=c(NA,NA),var_theta=c(NA,NA),Wexp=NA,Wobs=NA,ll=NA) })
        
        ### 1df for G (adjusted for Gcond)
        W1a <- resHa_a$theta[2]^2/resHa_a$var_theta[2]
        LR1a <- 2*(resHa_a$ll-resHo1_a$ll)
        S1a <- resHo1_a$Sobs
        
        ### 1df for Gcond (adjusted for G)
        W2a <- resHa_a$theta[3]^2/resHa_a$var_theta[3]
        LR2a <- 2*(resHa_a$ll-resHo2_a$ll)
        S2a <- resHo2_a$Sobs
        
        resi <- c( resHa_a$theta[2:3], resHa_a$var_theta[2:3], resHa_a$Wobs, resHo_a$Sobs, 2*(resHa_a$ll-resHo_a$ll), W1a,S1a,LR1a, W2a,S2a,LR2a  )
        names(resi) <- namscols
        assign(paste0("res",i), resi )
        
      }else assign(paste0("res",i), res_NULL)
      
    }
    
    res <- data.frame(p=p,Gpos=Gnam,Gcond=condGs$Gpos,samp=sampl,ss=n2,it=as.integer(it),itidx=as.integer(it_idx),Alloc=c("Complete","Combined","RDS","TZL","Opt.Lagr","Opt.GA"),rbind(out_com,res0,res1,res2,res3,res4))
    
    return(res)    
  } 


if( save_ind ) save(resOpt_conditional,file=paste0(savedir,"RealisticMultiDesBeta1_ss",n2,"_samp=",sampl,"_N=5K_R=",ifelse(Rep>=1000,paste0(Rep/1000,"K"),Rep),"_",optMetout,"_CondAnalysis.RData"))