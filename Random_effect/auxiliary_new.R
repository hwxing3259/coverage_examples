ooCORE<-function(PartSim,F,bn,sb,ZZ) {

CLf=PartSim$CLf
PL=PartSim$PL

print(EST)

if (EST=="LME") {
  print(c(F$FE,F$RE)); attributes(F$VF);
}
if (EST=="LMER") {
  print (F$FRE.LMER)
}
if (EST=="LMBF") {
  print(c(F$EF.LMBF,F$RE.LMBF,F$METHOD.LMBF,F$IS_SAMPLES.LMBF))
}

#if we interrupt the MCMC remove the unwanted NA's 
lst=which(!is.na(PL[,1])); CLf=CLf[lst,]; PL=PL[lst,]
if (is.VE) {CLv=CLv[lst,];}

#burnin and subsampling
mx=dim(PL)[1]; sl=seq(from=bn,to=mx,by=sb)
CLf=CLf[sl,]; PL=PL[sl,]
if (is.VE) {CLv=CLv[sl,];}

#mean(PL[,3]) #these two lines check MCMC targeting CRP has correct mean
#EK=sum(alpha/(alpha+0:(nl-1))); EK

###
#Superficial Output Analysis
print(effectiveSize(as.mcmc(PL)))
print(summary(as.mcmc(PL)))

###
#visualise co-occurence matrix p(i,j) = prob i,j in same cluster 
comf=CooM(CLf,ZZ)
if (is.VE) {comv=CooM(CLv,ZZ)}

print(sprintf("%s",apply(hpd(CL=CLf,lev=0.5,ZZ=ZZinclude)$set,1,paste,collapse=' ')))

###
#MAP cluster (max in run)
#mi=which.max(PL[,1]+PL[,2])
#prettyPart(mg2gm(CLf[mi,],ZZ))

#ML cluster (max in run)
#mj=which.max(PL[,1])
#prettyPart(mg2gm(CLf[mj,],ZZ))

}

#################################
#produce the desired design matirix, maybe even synthetic data
makeSdata<-function(dataset,thin,synth=FALSE,base=NA,glmm=NA,prior=NA,St=NA,uniq.cov=FALSE,cut.cov=NA,keep.cov=NA) {

  #load the data
  if (dataset=="PPP") {
    load(".RData"); 
    all=PPP_all  
  } else {
    load(".RData")
    all=HHH_all
  }

  #must remove "total" which is an aggregate over ZZ not an ZZ itself
  all=all[all$NNNN!="total",]

  b=all[is.element(all$NNNN,thin),]
  b$NNNN=factor(b$NNNN)
  n=dim(b)[1]

  #seems strange but SSSS is uneven offsets 0,120,300 and DDDD is binary 0,7 anyway
  #and pressure has just 5 levels and isnt linear so might as well go categorical there
  b$SSSS=as.factor(b$SSSS)
  if (dataset=="PPP") {
    b$DDDD=as.factor(b$DDDD)  
  } else {
    b$PPPS=as.factor(b$PPPS) 
  }
  
  #create batch index marker to handle repeated measurements analysis showed that a batch 
  #is the set of 3 observations with the same SSSS, DDDD, PPPC/PPPS and MMMM value 
  #(and any NNNN value) these show clear correlation

  #b has a factor indicating batch
  if (dataset=="PPP") {
    X=model.matrix(Value~SSSS+PPPC+DDDD+MMMM,data=b)  
  } else {
    X=model.matrix(Value~SSSS+PPPS+MMMM,data=b)  
  }
  p=dim(X)[2]
  v=2^(0:(p-1))
  b$batch=as.factor(X%*%v)                  #no fine resolution batches ... *10+(0:(n-1))%%3)
  #b$batch=as.factor(X%*%v*10+(0:(n-1))%%3) #if you want fine resolution batches

  #print(b[1:30,]) #check the batch labels are right

  #br has one sample from each batch - the first in each batch
  bi=seq(from=1,to=n,by=3)
  br=b[bi,]
  nr=dim(br)[1]

  #bd takes the average within batches - EDA showed they were very consistent 
  bd=br                                  #take the covariates from each batch
  bd$Value=apply(matrix(b$Value,nr,3,byrow=TRUE),1,mean)  #overwrite averages
  bd$Y=scale(log(bd$Value))
  n=dim(bd)[1]

  #take subset of the data - just include one sample 
  #for each unique compination of covariates.
  if (uniq.cov) {
    bd=bd[,!is.element(names(bd),cut.cov)]
    row.names(bd)<-as.character(1:n)
    bd=bd[as.numeric(row.names(unique(bd[,keep.cov]))),]
  }
  if (synth) {

    Bt=bd
    Bt$St=as.factor(gm2mg(St,Bt$NNNN)) #St defined in calling level
    b.gm<-MCMCglmm(GLMM$FX.SYN, random=GLMM$RE.SYN, rcov=GLMM$RC.SYN, data=Bt, prior=prior, pr=TRUE, verbose=FALSE, nitt=2000, burnin=0)

    if (base=="APPROXPOST") {
      lls=b.gm$Deviance[dim(b.gm$Sol)[1]]/(-2)
    } else {
      lls=NA
      mev=0*Bt$Y+10000
      b.gm.prior<-MCMCglmm(GLMM$FX.SYN, random=GLMM$RE.SYN, rcov=GLMM$RC.SYN, data=Bt, mev=mev, prior=prior, pr=TRUE, verbose=FALSE, nitt=2000, burnin=0)
      loc=c(); for (k in colnames(b.gm$Sol)) loc=c(loc,which(k==colnames(b.gm.prior$Sol)))
      b.gm$Sol=b.gm.prior$Sol[,loc]
      loc=c(); for (k in colnames(b.gm$VCV)) loc=c(loc,which(k==colnames(b.gm.prior$VCV)))
      b.gm$VCV=b.gm.prior$VCV[,loc]
    }
    Bt$Value=exp(simulate(b.gm,marginal=NULL,type="response",it=dim(b.gm$Sol)[1],seed = as.integer(314159)))
    Bt$Y=scale(log(Bt$Value))
    Bt$Sf=Bt$St
    bd=Bt   
  }
 
  return(bd)
}



############################################################

#CRP simulation
#returns string group label representation

PriorPartSimS<-function(n,N,ZZ,alpha) {

  CL=matrix(NA,N,length(ZZ))
  
  for (k in 1:N) {
    g2=DP.sim(alpha,n=length(ZZ))
    while (max(g2)<2) { #CAN NOT BE IN THE SAME CLUSTER
      g2=DP.sim(alpha,n=length(ZZ))
    }
    CL[k,]=g2
  }
  return(list(CLf=CL,PL=NA))

}

#thin=ZZ;synth=TRUE;St=Strue
#simulate from the given partition?
LkdYSim<-function(CL,n,N,ZZ,Y.dat,dataset,base,glmm,prior,uniq.cov=FALSE,cut.cov=NA,keep.cov=NA) {
  Y.sim=matrix(NA,n,N)
  d=rep(NA,N)
  for (k in 1:N) {
    Strue=mg2gm(CL[k,],ZZ)
    bd=makeSdata(dataset,thin=ZZ,synth=TRUE,base,glmm,prior,St=Strue,uniq.cov,cut.cov,keep.cov)
    Y.sim[,k]=bd$Y
    d[k]=norm(Y.sim[,k]-Y.dat)
    print(c(k,d[k]))
  }
  return(list(Y.sim=Y.sim,d=d))
}


######################################################################

runMCMC<-function(B,EST,is.VE,F,J,SS,WTFI,alpha=1,VARY.ALPHA=FALSE,w.alpha=NA,APR=NA,SAVE_TO_FILE=FALSE,outfile=NA,DEBUG=FALSE) {

n=dim(B)[1]

#levels of Amino Acids
ZZ=levels(B$NNNN); nl=length(ZZ)  

#Start state
Sf=as.list(ZZ) #initial groups for FE or RE clustering
Kf=length(Sf)
B$Sf=as.factor(gm2mg(Sf,B$NNNN))

#gm2mg: list of string representation to numeric label
#mg2gm: numeric label to list of string


if (is.VE) {
  Sv=as.list(ZZ) #initial groups for variance effect (VE) clusters
  Kv=length(Sv)  #initial number of VE clusters
  B$Sv=as.factor(gm2mg(Sv,B$NNNN))
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

hashing=array(NA,1); row.names(hashing)<-prettyPart(Sf,ZZ); hashing[1]=op;

#old log cluster prior
o.crp=log.pr.crp(Sf,alpha,nl)
if (is.VE) {o.crp=o.crp+log.pr.crp(Sv,alpha,nl)}

#Clf: group lablel representation of the partition
THf=list(); CLf=matrix(NA,J/SS,nl); 

#nSf: Cardnality of each set in the partition
#Af: the group lable for each acid
#kf: num of partitions
nSf=unlist(lapply(Sf,length)); Af=gm2mg(Sf,ZZ); Kf=length(Sf)
THf[[1]]=list(nSf); CLf[1,]=Af; 

if (!is.VE) {
  PL=matrix(NA,J/SS,4); 
  PL[1,]=c(alpha,op,o.crp,Kf) 
  THv=list(); CLv=matrix(NA,J/SS,nl); #not used just saves a test when we write
} else {
  THv=list(); CLv=matrix(NA,J/SS,nl);
  nSv=unlist(lapply(Sv,length)); Av=gm2mg(Sv,ZZ); Kv=length(Sv)
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

  #go through each ZZ and update its cluster membership in f and v
  for (i in 1:nl) {

    Spf=GenerateCandidate(ZZ,Sf,i)
    
    #numeric lable for membership
    Afp=gm2mg(Spf,ZZ)
    NEW=any(Afp!=Af)

    if ( (NEW && length(Spf)>1) || DEBUG) {
      if (is.VE) {B$Sv=as.factor(gm2mg(Sv,B$NNNN))}

      if (DEBUG) {
        np=op; b.mp=b.m
      } else {
        
        
        
        
      np=NA
	if (EST=="LME") {
        B$Sf=as.factor(gm2mg(Spf,B$NNNN))
  	  b.mp=lme(fixed=F$FE,random=F$RE,weights=F$VF,data=B,method=F$MT,control=BOB) 
  	  np=BIC(b.mp)/(-2)
	}
	if (EST=="LMER") {
	  #if ppSpf is different from the previous iteration (encoded in rowname), 
	  #then refit the model, compute BIC
	  #np is na if hashing does not store the partition before
	  #nt that the first line assign values to np if np is not NA
        if (is.na(np<-hashing[ppSpf<-prettyPart(Spf,ZZ)])) {
          B$Sf=as.factor(gm2mg(Spf,B$NNNN))
        #  print(prettyPart(Spf,ZZ))
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
        B$Sf=as.factor(gm2mg(Sf,B$NNNN))
        b.m=lmBF(formula=F$EF.LMBF, data=B, whichRandom=F$RE.LMBF, progress=FALSE, method=F$METHOD.LMBF, iterations=F$IS_SAMPLES.LMBF)
        op=attributes(b.m)[[3]]$bf
        #cal new
        B$Sf=as.factor(gm2mg(Spf,B$NNNN))
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
        Af=gm2mg(Sf,ZZ)

        b.m=b.mp           #not really needed
      }
    }  

    ####not needed########
    
    if (is.VE) {

      Spv=GenerateCandidate(ZZ,Sv,i)

      if (length(Spv)>1) {

        B$Sf=as.factor(gm2mg(Sf,B$NNNN))

        np=NA
	  if (EST=="LME") {
          B$Sv=as.factor(gm2mg(Spv,B$NNNN))
  	    b.mp=lme(fixed=F$FE,random=RE,weights=F$VF,data=B,method=F$MT,control=BOB) 
  	    np=BIC(b.mp)/(-2)
	  }
	  if (EST=="LMER") {
          B$Sv=as.factor(gm2mg(Spv,B$NNNN))
  	    b.mp=lmer(formula=F$FRE.LMER,data=B,REML=F$MT.LMER,control=DAVE)
  	    np=BIC(b.mp)/(-2)
	  }
	  if (EST=="LMBF") {
          #this is Approximate penalty method - refresh both estimates
          B$Sv=as.factor(gm2mg(Sv,B$NNNN))
          b.m=lmBF(formula=F$EF.LMBF, data=B, whichRandom=F$RE.LMBF, progress=FALSE, method=F$METHOD.LMBF, iterations=F$IS_SAMPLES.LMBF)
          op=attributes(b.m)[[3]]$bf
          B$Sv=as.factor(gm2mg(Spv,B$NNNN))
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
          Av=gm2mg(Sv,ZZ)

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
