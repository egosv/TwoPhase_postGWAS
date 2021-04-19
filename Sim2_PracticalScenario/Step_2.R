
resOpt_single <- foreach(G=iter(genotype_noZ,by='column'), p=icount(), .combine='rbind', .verbose=T) %:%
  foreach(data_iter=isplit( dat_sim_all, list(it=dat_sim_all$iteration) ), ids_iter=isplit( dat_ids_all, list(it=dat_ids_all$it) ), .combine='rbind', .inorder=F, .verbose=T, .packages = c("statmod","twoPhaseGAS")) %dopar% {
    
    it <- paste(data_iter$key,collapse="_")
    dat_sim0 <- data_iter$value
    
    it_idx <- paste(ids_iter$key,collapse="_")
    indx_it <- ids_iter$value
    
    Gnam <- colnames(G) 
    colnames(G) <- "G"
    dat_sim <- cbind(dat_sim0,Z=Z,G)
    dat_sim$fZ <- factor(dat_sim$Z)
    
    ### Complete data analysis 
    lmfit.com <- glm(formula_Ha, data=dat_sim, family=gaussian); lmfit.com.Ho <- glm(formula_Ho, data=dat_sim, family=gaussian)
    
    ### Extract ids
    R0 <- unlist(indx_it$R0)
    R1 <- unlist(indx_it$R1)
    R2 <- unlist(indx_it$R2)
    R3 <- unlist(indx_it$R3)
    R4 <- unlist(indx_it$R4)
  
    ### SMPL Analyses
    for( i in c(0:4)){ # i=1
      cond = eval(parse(text=paste0("R",i,"==1")))
      data1 <- dat_sim[cond,c("Y","G","Z","fZ","S")] # phase 2 data
      data0 <- dat_sim[!cond,c("Y","Z","fZ","S")] # phase 2 data complement
      
      resHo <- tryCatch(twoPhaseSPML(formula=Y~fZ,miscov=~G,auxvar=~Z,family=gaussian,data0,data1,verbose=FALSE), error=function(e){print(e); list(Sobs=NA,ll=NA) })
      resHa <- tryCatch(twoPhaseSPML(formula=Y~G+fZ,miscov=~G,auxvar=~Z,family=gaussian,data0,data1,verbose=FALSE), error=function(e){print(e); list(theta=rep(NA,4),var_theta=rep(NA,4),Wobs=NA,ll=NA) })
      
      assign(paste0("res",i), c( resHa$theta[2], resHa$var_theta[2], resHa$Wobs, resHo$Sobs, 2*(resHa$ll-resHo$ll) ) )
    }
    
    
    out_com0 <- c(coef(lmfit.com)[2], diag(vcov(lmfit.com))[2])
    Wcom <- out_com0[1]^2/out_com0[2]
    LRcom <- as.numeric(2*(logLik(lmfit.com) - logLik(lmfit.com.Ho)))
    Scom <- glm.scoretest(lmfit.com.Ho,dat_sim[,"G"])^2
    
    out_com <- c( out_com0, Wcom, Scom, LRcom)
    
    names(out_com)<-c("beta1","var_beta1","W","S","LR")
    
    
    res <- data.frame(p=p,Gpos=Gnam,samp=sampl,ss=n2,it=as.integer(it),itidx=as.integer(it_idx),Alloc=c("Complete","Combined","RDS","TZL","LM","GA"),rbind(out_com,res0,res1,res2,res3,res4))
    
    return(res)    
  } 


if( save_ind ) save(resOpt_single,file=paste0(savedir,"RealisticMultiDesBeta1_ss",n2,"_samp=",sampl,"_N=5K_R=",ifelse(Rep>=1000,paste0(Rep/1000,"K"),Rep),"_",optMetout,"_SingleSNP.RData")) 
