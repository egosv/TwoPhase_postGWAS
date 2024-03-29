
### because we don't need to run the optimal designs for each G but rather for each Y, Z (fixed) combination, we split in two separate loops.
if(0){ # For testing only
  it <- sample(Rep,1)
  dat_sim0 <- subset(dat_sim_all, iteration==it)
}
des_ids_all <- foreach(data_iter=isplit( dat_sim_all_sub, list(it=dat_sim_all_sub$iteration) ), .combine='c', .inorder=F, .noexport=NULL, .verbose=T, .packages = c("twoPhaseGAS")) %dopar% {
  
  it <- paste(data_iter$key,collapse="_")
  dat_sim0 <- data_iter$value
  
  dat_sim <- cbind(dat_sim0,Z=Z)
  dat_sim$fZ <- factor(dat_sim$Z)
  
  olsfit = lm(Y~fZ,dat_sim)
  disp_ols = summary(olsfit)$sigma
  resy = residuals(olsfit)
  
  beta0 <- if( optM==1 ){
    c(coef(olsfit)[1], 0, coef(olsfit)[-1]) 
  }else{
    c(coef(olsfit)[1], 0, 0, coef(olsfit)[-1]) 
  }
  
  # RDS
  order_resy = order(resy)
  phase2_id_rds = c(order_resy[1:(n2/2)], order_resy[(N-(n2/2)+1):N])
  R1 <- 1*( 1:N %in% phase2_id_rds )
  
  data <- dat_sim[,c("Y","Z","fZ","S")] # all data without G (phase 1 data so to speak)
  numsce <- nrow(maf_ld_comb)
  
  ### Precompute IMs (and p_gz0)
  IMlist <- lapply(1:numsce, function(j){
    
    ### calculate p_gz using AF and LD and available data from Z
    if( optM==1 ){
      pZ <- 1-mean(Z)/2
      pG <- 1-maf_ld_comb$P_g[j]
      LDj <- cbind(1,2, maf_ld_comb$LD.r[j])
      p_gz0 <- suppressMessages(multilocus.Pgz(2, p=c(pZ, pG), LD = LDj))

    }else{
      pZ <- 1-mean(Z)/2
      pG1 <- 1-maf_ld_comb$P_g1[j]
      pG2 <- 1-maf_ld_comb$P_g2[j]
      LDj <- rbind(cbind(1,2, maf_ld_comb$LD.rZG1[j]), 
                   cbind(1,3, maf_ld_comb$LD.rZG2[j]), 
                   cbind(2,3, maf_ld_comb$LD.rG1G2[j]) )
      p_gz0 <- suppressMessages(multilocus.Pgz(3, p=c(pZ, pG1, pG2), LD = LDj))
      
    }
    
    ## error objects for the optimal allocations
    R3.err <- 0; R4.err <-0
    
    
    IM <- obsIM.R(formula=formula_Ha, miscov=miscov_formula, auxvar=~Z, family=gaussian, data, beta=beta0, p_gz=p_gz0, disp=disp_ols)
    
    list(IM=IM,p_gz0=p_gz0)
  })
  
  ### Compute the designs 
  indx_it <- lapply(1:numsce, function(j){
    # j <- 1
    IM <- IMlist[[j]]$IM
    p_gz0 <- IMlist[[j]]$p_gz0
    
    ## error objects for the optimal allocations
    R3.err <- 0; R4.err <-0
    
    # data.frame(xtabs(strataformula,dat_sim))
    # TZL - OPT
    ## var G/Z
    p_z <- aggregate(q~Z, data = p_gz0, FUN=sum)
    names(p_z)[2] <- "Freq"
    p_g_z <- merge(p_gz0, p_z, by="Z")
    p_g_z$qcond <- p_g_z$q/p_g_z$Freq
    ### compute var(G|Z) [or the corresponding optimality criterion if a matrix]
    varG_Z = do.call(rbind, by(p_g_z, INDICES=p_g_z[,"Z"],function(x){
      Xmat <- x[, Gs, drop=F]*sqrt(x$qcond)
      EXmat <- colSums(x[, Gs, drop=F]*x$qcond)
      vcovmat <- t(Xmat) %*% as.matrix(Xmat) - EXmat %*% t(EXmat)
      
      sumval <- twoPhaseGAS:::.optMsure(optMethod)(vcovmat, Kind)
      c(Z=x[1,"Z"],var=sumval)
    }))
    ### scale the residuals
    resyopt <- if( optMethod=="Par-spec" ){
      resy*sqrt(varG_Z[match(dat_sim$Z, varG_Z[,"Z"]),"var"])
    }else{
      resy*varG_Z[match(dat_sim$Z, varG_Z[,"Z"]),"var"]
    }
    order_resyopt = order(resyopt)
    phase2_id_opt = c(order_resyopt[1:(n2/2)], order_resyopt[(N-(n2/2)+1):N])
    R2 <- 1*( 1:N %in% phase2_id_opt )
    
    ## balanced
    pop1a <- sapply(1:100,function(x){
      sa <- which(BSS(samp.fracs=sel.prob.bal,n2,data,strataformula)$R==1)
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
      sa <- which(BSS(samp.fracs=sel.prob.com,n2,data,strataformula)$R==1)
      len <- length(sa)
      if( len==n2 ){ return(sa)
      }else if( len>n2 ){ return( sa[-sample(len, len-n2)] )
      }else return( c(sa,sample((1:N)[-sa],n2-len)) )
    })
    fit2 <- apply(pop2a,2,function(r){
      return(fitnessTP(IM,1*(1:N %in% r),optimMeasure=optMethod,K.idx=Kind))
    })
    R0 <- 1*(1:N %in% pop2a[,which(fit2 == min(fit2) )])
    
    opt.prop <- tryCatch(optimTP.LM(formula=formula_Ha,miscov=miscov_formula,auxvar=~Z,strata=strataformula,family=gaussian,n=n2,data,beta=beta0,p_gz=p_gz0,disp=disp_ols,optimMeasure=optMethod,K.idx=Kind),error=function(e){print(e); cat("Error in optimTP.LM on iteration:",it,", using same sampling fraction for all strata instead.\n");
      R3.err <<- 1;
      stratadf <- data.frame(xtabs(strataformula,data=data))
      stratadf$prR_cond_optim <- n2/N; return(stratadf)})
    
    pop3a <- sapply(1:100,function(x){
      sa <- which(BSS(samp.fracs=opt.prop$prR_cond_optim,n2,data,strataformula)$R==1)
      len <- length(sa)
      if( len==n2 ){ return(sa)
      }else if( len>n2 ){ return( sa[-sample(len, len-n2)] )
      }else return( c(sa,sample((1:N)[-sa],n2-len)) )
    })
    fit3 <- apply(pop3a,2,function(r){
      return(fitnessTP(IM,1*(1:N %in% r),optimMeasure=optMethod,K.idx=Kind))
    })
    R3 <- 1*(1:N %in% pop3a[,which(fit3 == min(fit3) )])
    
    ## instead of providing a random initialization, give an informed initialization based on combined, balanced and optimal (via Lagrange multipliers) allocations (a third each)
    ## balanced
    pop1 <- pop1a[,order(fit1)[1:20]]
    ## combined
    pop2 <- pop2a[,order(fit2)[1:20]]
    ## LM
    pop3 <- pop3a[,order(fit3)[1:19]]
    
    GA.sol <- tryCatch(optimTP.GA(ncores=1,formula=formula_Ha,miscov=miscov_formula,auxvar=~Z,family=gaussian,n=n2,data,beta=beta0,p_gz=p_gz0,disp=NULL,ga.popsize=60,ga.propelit=0.9,ga.proptourney=0.9,ga.ngen=500,ga.mutrate=0.001,ga.initpop=t(cbind(pop1,pop2,pop3,phase2_id_rds)),optimMeasure=optMethod,K.idx=Kind,seed=1),error=function(e){print(e); cat("Error in optimJTC_GA on iteration:",it,", using SRS instead.\n");
      R4.err <<- 1;
      return(list(bestsol=sample(N,n2)))})
    # library(kofnGA); plot(GA.sol)
    R4 <- rep(0,nrow(data)); R4[GA.sol$bestsol] <- 1

    list(R0=R0,R2=R2,R3=R3,R4=R4)
  })
  
  ### Evaluate the designs
  FitnessArray <- array(NA, dim=c(4,numsce,numsce))
  dimnames(FitnessArray)[[1]] <- c("Combined","TZL","Opt.Lagr","Opt.GA")
  dimnames(FitnessArray)[[2]] <- paste0("Comb",1:numsce)
  dimnames(FitnessArray)[[3]] <- paste0("Val",1:numsce)
  
  for ( j in 1:numsce ){ # j <- 1
    R0 <- indx_it[[j]]$R0
    R2 <- indx_it[[j]]$R2
    R3 <- indx_it[[j]]$R3
    R4 <- indx_it[[j]]$R4
    for( k in 1:numsce ){ # k <- 2
      IM_b <- IMlist[[k]]$IM
      
      FitnessArray[1,j,k] <- fitnessTP(IM_b,Rj=R0,optimMeasure=optMethod,K.idx=Kind)
      FitnessArray[2,j,k] <- fitnessTP(IM_b,Rj=R2,optimMeasure=optMethod,K.idx=Kind)
      FitnessArray[3,j,k] <- fitnessTP(IM_b,Rj=R3,optimMeasure=optMethod,K.idx=Kind)
      FitnessArray[4,j,k] <- fitnessTP(IM_b,Rj=R4,optimMeasure=optMethod,K.idx=Kind)
    }
  }
  
  indx <- apply(apply(FitnessArray,c(1,2), median),1,function(x){
    x1 <- which(x==min(x))
    return(x1[1])
  })
  
  R0 <- indx_it[[indx[1]]]$R0
  R2 <- indx_it[[indx[2]]]$R2
  R3 <- indx_it[[indx[3]]]$R3
  R4 <- indx_it[[indx[4]]]$R4
  
  data$R0 <- R0
  data$R1 <- R1
  data$R2 <- R2
  data$R3 <- R3
  data$R4 <- R4
  
  df <- data.frame(xtabs(~Z+S,data=data))
  xn <- setNames(df$Freq, paste0("Z",df$Z,"-","S",df$S)) ## alldata
  
  for( i in c(0:4) ){ # i=0
    assign(paste0("opt",i), data.frame(xtabs(as.formula(paste0("~Z+S+R",i)),data=data)))
    assign(paste0("opt",i), eval(parse(text=paste0("subset(opt",i,", R",i,"==1)"))))
    assign(paste0("x",i), eval(parse(text=paste0("setNames(opt",i,"$Freq, paste0('Z',opt",i,"$Z,'-','S',opt",i,"$S))"))))
  }
  
  
  alloc.res <- rbind(c("N_n2"=nrow(dat_sim), error=NA, xn), 
                     c("N_n2"=sum(R0), error=NA, x0),
                     c("N_n2"=sum(R1), error=NA, x1),
                     c("N_n2"=sum(R2), error=NA, x2),
                     c("N_n2"=sum(R3), error=NA, x3),
                     c("N_n2"=sum(R4), error=NA, x4)
  )
  alloc.res <- data.frame(it=it,samp=sampl,ss=n2,
                          alloc=c("Complete","Combined","RDS","TZL","LM","GA"),alloc.res)
  
  alloc_ids <- data.frame(it=as.numeric(it))
  alloc_ids$R0 = list(R0)
  alloc_ids$R1 = list(R1)
  alloc_ids$R2 = list(R2)
  alloc_ids$R3 = list(R3)
  alloc_ids$R4 = list(R4)
  
  return(list(list(alloc_df=alloc.res, alloc_ids=alloc_ids)))  
} 

### sort output
dat_alloc <- rbindlist( lapply(des_ids_all, `[[`, 1) )
dat_ids_all <- rbindlist( lapply(des_ids_all, `[[`, 2) )
setorder(dat_ids_all, it)
dat_ids_all[,1]
class(dat_alloc) <- "data.frame"
class(dat_ids_all) <- "data.frame"

### save design ids
if(save_ind) save(dat_alloc,dat_ids_all,file=paste0(savedir,"RealisticMultiDesBeta1_ss",n2,"_samp=",sampl,"_N=5K_R=",ifelse(Rep>=1000,paste0(Rep/1000,"K"),Rep),"_",optMetout,"_desigids.RData")) 

