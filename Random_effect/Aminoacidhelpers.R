#dpriior is CRP
#log cpr probability
log.pr.crp<-function(S,alpha,n) {
  #log CRP cluster probability
  nS=unlist(lapply(S,length))
  K=length(S)
  pr=K*log(alpha)+sum(lgamma(nS))-lgamma(alpha+n)+lgamma(alpha) #dont need if alpha fixed
  return(pr)
}

str2lst <- function(strinput){
  return(strsplit(sub(pattern = "(", replacement = "", 
                      sub(pattern = ",(", replacement = "", 
                          strsplit(strinput, ")" )[[1]], 
                          fixed = T), fixed = T),","))
}



#log likelihood for all partitions in anchor
#for K fixed partition, compute tilde(p)(y|s), compare with tilde(p)(y_obs|s)
logpr.joint <- function(ysynth, anchor, B){
  N = length(anchor)
  logjoint = matrix(NA, ncol = 2, nrow = N)
  ii = 1
  for (s in anchor){
    gm = str2lst(s)
    Sf = lapply(gm, function(i){return(sub(pattern = "A", replacement = "a", x = i))})
    B$Sf = as.factor(gm2mg(Sf, B$Name))
    B$Y = ysynth
    b.m=lmer(formula = Y ~ Day * Stage + (-1 + Stage | Sf) + (-1 + Stage | Name),
             data = B, REML = T, control = DAVE)
    ps=BIC(b.m)/(-2)
    logpcrp = log.pr.crp(Sf,1,12)
    logjoint[ii,] = c(logpcrp, ps) 
    #print(c(ps, logpcrp))
    ii = ii+1
  }
  return(logjoint)
}


##simulating y from prior given partition (list form)
priorysimulate <- function(B, listpart){
  aalable = levels(B$Name)
  #listpart = mg2gm(part, aalable)
  B$Sf = as.factor(gm2mg(listpart, B$Name))
  b.m=lmer(formula = Y ~ Day * Stage + (-1 + Stage | Sf) + (-1 + Stage | Name),
           data = B, REML = T)
  modelmtrx = cbind(getME(b.m, "X"), as.matrix(getME(b.m, "Z")))
  npart = length(listpart)
  #fixed effect: normal(0,1) prior for all 6 variables
  beta = rnorm(6)
  G1cov = rIW(V = diag(0.25,3), nu = 5)
  G2cov = rIW(V = diag(0.25,3), nu = 5)
  G1effect = rmvnorm(12, rep(0,3), G1cov)
  G2effect = rmvnorm(npart, rep(0,3), G2cov)
  finalbeta = c(beta, as.vector(t(G1effect)), as.vector(t(G2effect)))
  sig = 1/rgamma(1,8,1) #equivalent to IG
  res = rnorm(dim(B)[1], mean = 0, sd = sig)
  priory = modelmtrx %*% finalbeta + res
  return(priory)
}


#simulate pair of (S,y_i)
priorfunc <- function(B){
  part = DP.sim(1,12)
  while(max(part) == 1) {part = DP.sim(1,12)}
  ysim = priorysimulate(B, mg2gm(part, levels(B$Name)))
  return(list(part,ysim))
}


#sampling from the approx.post for partition and 
#find the approximate credible set for synthetic data y_i
runMCMCmodified<-function(B,EST,is.VE,F,J,SS,WTFI,alpha=1,VARY.ALPHA=FALSE,w.alpha=NA,APR=NA,SAVE_TO_FILE=FALSE,outfile=NA,DEBUG=FALSE,PART) {
  
  n=dim(B)[1]
  
  #levels of Amino Acids
  AA=levels(B$Name); nl=length(AA)  
  Sf= PART #initial groups for FE or RE clustering
  Kf=length(Sf)
  B$Sf=as.factor(gm2mg(Sf,B$Name))
  #gm2mg: list of string representation to numeric label
  #mg2gm: numeric label to list of string
  
  if (is.VE) {
    Sv=as.list(AA) #initial groups for variance effect (VE) clusters
    Kv=length(Sv)  #initial number of VE clusters
    B$Sv=as.factor(gm2mg(Sv,B$Name))
  }
  
  op=NA
  if (EST=="LME") {
    b.m=lme(fixed=F$FE,random=F$RE,weights=F$VF,data=B,method=F$MT,control=BOB) 
    op=BIC(b.m)/(-2)
  }
  if (EST=="LMER") {
    b.m=lmer(formula=F$FRE.LMER,data=B,REML=F$MT.LMER,control=DAVE)
    op=BIC(b.m)/(-2)
  }
  if (EST=="LMBF") {
    b.m=lmBF(formula=F$EF.LMBF, data=B, whichRandom=F$RE.LMBF, progress=FALSE, method=F$METHOD.LMBF, iterations=F$IS_SAMPLES.LMBF)
    op=attributes(b.m)[[3]]$bf
  }
  if (is.na(op)) {error("ML estimator EST value non recognised. Use one of LME, LMER or LMBF.")}
  
  hashing=array(NA,1); row.names(hashing)<-prettyPart(Sf,AA); hashing[1]=op;
  
  #old log cluster prior
  o.crp=log.pr.crp(Sf,alpha,nl)
  if (is.VE) {o.crp=o.crp+log.pr.crp(Sv,alpha,nl)}
  
  #Clf: group lablel representation of the partition
  THf=list(); CLf=matrix(NA,J/SS,nl); 
  
  #nSf: Cardnality of each set in the partition
  #Af: the group lable for each acid
  #kf: num of partitions
  nSf=unlist(lapply(Sf,length)); Af=gm2mg(Sf,AA); Kf=length(Sf)
  THf[[1]]=list(nSf); CLf[1,]=Af; 
  
  if (!is.VE) {
    PL=matrix(NA,J/SS,4); 
    PL[1,]=c(alpha,op,o.crp,Kf) 
    THv=list(); CLv=matrix(NA,J/SS,nl); #not used just saves a test when we write
  } else {
    THv=list(); CLv=matrix(NA,J/SS,nl);
    nSv=unlist(lapply(Sv,length)); Av=gm2mg(Sv,AA); Kv=length(Sv)
    THv[[1]]=list(nSv); CLv[1,]=Av; 
    PL=matrix(NA,J/SS,5); 
    PL[1,]=c(alpha,op,o.crp,Kf,Kv) 
  }
  
  #MCMC J steps each step sweeps all 2xnl cluster memberships
  for (j in 1:J) {
    
    #write info about current state
    o1=paste(sprintf('%g',c(j,alpha,op,o.crp)),collapse=' ')
    o2=paste(sprintf('%g',c(Kf,nSf)),collapse=' ')
    if (!is.VE) {
      print(paste(c(o1,o2),collapse=' # '))
    } else {
      o3=paste(sprintf('%g',c(Kv,nSv)),collapse=' ')
      print(paste(c(o1,o2,o3),collapse=' # '))
    }
    
    #MCMC SRW update for alpha the CRP parameter
    if (VARY.ALPHA) {
      ap=abs(alpha+w.alpha*(2*runif(1)-1))
      n.crp=log.pr.crp(Sf,ap,nl)
      MHR = n.crp - o.crp + dexp(ap,rate=APR,log=TRUE) - dexp(alpha,rate=APR,log=TRUE)
      
      if (log(runif(1))<MHR) {
        alpha=ap
        o.crp=n.crp
      }
    }
    
    #go through each AA and update its cluster membership in f and v
    for (i in 1:nl) {
      
      Spf=GenerateCandidate(AA,Sf,i)
      
      #numeric lable for membership
      Afp=gm2mg(Spf,AA)
      NEW=any(Afp!=Af)
      
      if ( (NEW && length(Spf)>1) || DEBUG) {
        if (is.VE) {B$Sv=as.factor(gm2mg(Sv,B$Name))}
        
        if (DEBUG) {
          np=op; b.mp=b.m
        } else {
          
          
          
          
          np=NA
          if (EST=="LME") {
            B$Sf=as.factor(gm2mg(Spf,B$Name))
            b.mp=lme(fixed=F$FE,random=F$RE,weights=F$VF,data=B,method=F$MT,control=BOB) 
            np=BIC(b.mp)/(-2)
          }
          if (EST=="LMER") {
            #if ppSpf is different from the previous iteration (encoded in rowname), 
            #then refit the model, compute BIC
            #np is na if hashing does not store the partition before
            #nt that the first line assign values to np if np is not NA
            if (is.na(np<-hashing[ppSpf<-prettyPart(Spf,AA)])) {
              B$Sf=as.factor(gm2mg(Spf,B$Name))
              #  print(prettyPart(Spf,AA))
              b.mp=lmer(formula=F$FRE.LMER,data=B,REML=F$MT.LMER,control=DAVE)
              np=BIC(b.mp)/(-2)
              # print(np)
              #construct new_hashing, using string partition as rowname
              new_hashing=array(NA,1); row.names(new_hashing)<-ppSpf; new_hashing[1]=np;
              #record all old rownames(partition strings), add the new one.
              oldnames=row.names(hashing); hashing=rbind(hashing,new_hashing); row.names(hashing)<-c(oldnames,ppSpf);
            } else {print(hashing[ppSpf])}
          }
          if (EST=="LMBF") {
            #this is Approximate penalty method - refresh both estimates
            #recal old
            B$Sf=as.factor(gm2mg(Sf,B$Name))
            b.m=lmBF(formula=F$EF.LMBF, data=B, whichRandom=F$RE.LMBF, progress=FALSE, method=F$METHOD.LMBF, iterations=F$IS_SAMPLES.LMBF)
            op=attributes(b.m)[[3]]$bf
            #cal new
            B$Sf=as.factor(gm2mg(Spf,B$Name))
            b.mp=lmBF(formula=F$EF.LMBF, data=B, whichRandom=F$RE.LMBF, progress=FALSE, method=F$METHOD.LMBF, iterations=F$IS_SAMPLES.LMBF)
            np=attributes(b.mp)[[3]]$bf
          }
          if (is.na(np)) {error("ML estimator EST value non recognised. Use one of LME, LMER or LMBF.")}
        }
        
        
        
        
        
        
        
        #todo - do the cancelation in the MHR
        n.crp=log.pr.crp(Spf,alpha,nl)
        if (is.VE) {n.crp=n.crp+log.pr.crp(Sv,alpha,nl)} 
        
        MHR = np - op + n.crp - o.crp
        #MHR = n.crp - o.crp #CRP prior target for debugging
        
        if (log(runif(1))<MHR) {
          Sf=Spf
          o.crp=n.crp
          op=np
          
          #todo - in todo's above, nSp and Kp to be computed above
          nSf=unlist(lapply(Sf,length))
          Kf=length(Sf)
          Af=gm2mg(Sf,AA)
          
          b.m=b.mp           #not really needed
        }
      }  
      
      ####not needed########
      
      if (is.VE) {
        
        Spv=GenerateCandidate(AA,Sv,i)
        
        if (length(Spv)>1) {
          
          B$Sf=as.factor(gm2mg(Sf,B$Name))
          
          np=NA
          if (EST=="LME") {
            B$Sv=as.factor(gm2mg(Spv,B$Name))
            b.mp=lme(fixed=F$FE,random=RE,weights=F$VF,data=B,method=F$MT,control=BOB) 
            np=BIC(b.mp)/(-2)
          }
          if (EST=="LMER") {
            B$Sv=as.factor(gm2mg(Spv,B$Name))
            b.mp=lmer(formula=F$FRE.LMER,data=B,REML=F$MT.LMER,control=DAVE)
            np=BIC(b.mp)/(-2)
          }
          if (EST=="LMBF") {
            #this is Approximate penalty method - refresh both estimates
            B$Sv=as.factor(gm2mg(Sv,B$Name))
            b.m=lmBF(formula=F$EF.LMBF, data=B, whichRandom=F$RE.LMBF, progress=FALSE, method=F$METHOD.LMBF, iterations=F$IS_SAMPLES.LMBF)
            op=attributes(b.m)[[3]]$bf
            B$Sv=as.factor(gm2mg(Spv,B$Name))
            b.mp=lmBF(formula=F$EF.LMBF, data=B, whichRandom=F$RE.LMBF, progress=FALSE, method=F$METHOD.LMBF, iterations=F$IS_SAMPLES.LMBF)
            np=attributes(b.mp)[[3]]$bf
          }
          
          if (is.na(np)) {error("ML estimator EST value non recognised. Use one of LME, LMER or LMBF.")}
          
          #todo - do the cancelation in the MHR
          n.crp=log.pr.crp(Sf,alpha,nl)+log.pr.crp(Spv,alpha,nl)
          
          MHR = np - op + n.crp - o.crp
          #MHR = n.crp - o.crp #CRP prior target for debugging
          
          if (log(runif(1))<MHR) {
            Sv=Spv
            o.crp=n.crp
            op=np
            
            #todo - in todo's above, nSp and Kp to be computed above
            nSv=unlist(lapply(Sv,length))
            Kv=length(Sv)
            Av=gm2mg(Sv,AA)
            
            b.m=b.mp           #not really needed
          }
        }  
      }
    }
    
    #collect samples from the MCMC every SS steps
    if (j%%SS==0) {
      #conditional on the number in each cluster we know the dbn of the cluster weights w 
      THf[[j/SS]]=list(nSf)
      CLf[j/SS,]=Af; 
      if (!is.VE) {
        PL[j/SS,]=c(alpha,op,o.crp,Kf)
      } else {
        THv[[j/SS]]=list(nSv)
        CLv[j/SS,]=Av; 
        PL[j/SS,]=c(alpha,op,o.crp,Kf,Kv)
      }   
      if (SAVE_TO_FILE && ((j%%(WTFI*SS))==0)) save(CLf,PL,THf,CLv,THv,file=outfile)
    }
  }
  
  #END MCMC
  
  return(list(CLf=CLf,PL=PL,THf=THf,CLv=CLv,THv=THv))
  
}



GenerateCandidate2 <- function(AA,S,i){
  part = GenerateCandidate(AA,S,i)
  #print(length(part))
  while (length(part) == 1) {part = GenerateCandidate(AA,S,i)}
  return(part)
}



#get the 17 dimensional summary statistics and (c_i,s_i) pair from raw outputs
partitionreg17 <- function(tryregression){
  ci = rep(NA, length(tryregression))
  Si = matrix(NA, nrow = length(tryregression), ncol = 17)
  for (i in 1:length(tryregression)){
    #print(i)
    if (!is.null(tryregression[[i]][[1]])){
      ci[i] = tryregression[[i]][[1]]
    }
    
    if (!is.na(ci[i])){
      Si[i,] = tryregression[[i]][[2]]
    }
  }
  
  #remove the NAs
  naindica = which(is.na(ci))
  
  ci = ci[-naindica]
  Si = Si[-naindica,]
  return(cbind(ci,Si))
}



#the 72 dimensional summary statistics and (c_i,s_i) pair from raw outputs
partitionreg72 <- function(tryregression,B){
  ci = rep(NA, length(tryregression))
  Si = matrix(NA, nrow = length(tryregression), ncol = 72)
  for (i in 1:length(tryregression)){
    #print(i)
    if (!is.null(tryregression[[i]][[1]])){
      ci[i] = tryregression[[i]][[1]]
    }
    
    if (!is.na(ci[i])){
      Bsim = B
      Bsim$Y = tryregression[[i]][[3]]
      Si[i,] = aggregate(Bsim$Y, list(Bsim$Stage, Bsim$Day, Bsim$Name),mean)[,4]
    }
  }
  
  naindica = which(is.na(ci))
  
  ci = ci[-naindica]
  Si = Si[-naindica,]
  return(cbind(ci,Si))
}


#legal permutation of y corresponds to relabling the data
#since the design is balanced and complete, do not change the result
quadrandseq <- function(acidobs,B){
  quadrant = list()
  for (i in c("Day", "Stage", "Name")){
    factorvec = B[,i]
    quadrantseq = rep(NA, length(levels(factorvec)))
    l = 1
    for (j in levels(factorvec)){
      quadrantseq[l] = mean(acidobs[which(factorvec == j)])
      l = l+1
    }
    quadrant[[i]] = quadrantseq
  }
  
  return(quadrant)
}
