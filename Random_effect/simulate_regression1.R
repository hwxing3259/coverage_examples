
source("setup.R")
library(MCMCglmm)
library(mvtnorm)
library(optimx)
#local helper functions
source("auxiliary.R")
source("auxiliary_new.R")
source("Aminoacidhelpers.R")

#Generate samples from the base distribution over partitions
#(prior or approx prior)
BASE.DBN="PRIOR" #"PRIOR"  #"APPROXPOST" #basedbn is prior

#In the fitting stage are we using BIC (in lmer) or lmBF approx?
EST="LMER" #"LMER" # or "LMBF"

#where is all this going?
job=paste0(BASE.DBN,"-",EST,"regression")
gather.dir=paste(job,"/",sep="")
if (!dir.exists(gather.dir)) {
  dir.create(gather.dir)
}

#do we thin the data?
#everything
kill=c()
AAinclude=c("PHE","LEU","TRP","TYR","MET","LYS","ASP","GLU",
            "ASN","GLN","HIS","SER","THR","AIA","GLY","ILE","PRO","VAL","ORN")
#just 12
kill=c("total","ASN","ASP","HIS","GLN","SER","ORN","PRO")
AAinclude=setdiff(AAinclude,kill)
#just 6
#kill=c(kill,"AIA","TYR","MET","LYS","GLU","THR",
#AAinclude=setdiff(AAinclude,kill)

#do we thin the data for unique covariate patterns?
UNIQ.COV=FALSE; COV.REMOVE=c("Value","batch","Process","Meat"); COV.KEEP=c("Stage","Day","Name")

#choice of alpha based on mean #clusters 3-4
alpha=1; VARY.ALPHA=FALSE; w.alpha=NA; APR=NA


#in approx Posterior analysis need number steps for phi-side simulation
#give run parameters, output storage, hard disk write frequency
Jphi=2000; SS=1; WTFI=Jphi/10; DEBUG=FALSE

#where do we write? 
outfile.set=paste(gather.dir,"settings-",BASE.DBN,"-",EST,".RData",sep="")            #settings
outfile.phi=paste(gather.dir,"phiside-",BASE.DBN,"-",EST,".RData",sep="")            #after the run
outfile.phi.BAC=paste(gather.dir,"phiside-",BASE.DBN,"-",EST,"-BAC",".RData",sep="") #during the run

#which data set is it
EXP="PEF" #or "HPP" tho that isnt fully implemented in this file it is fine elsewhere

#Model for data given partition St
GLMM=list()
GLMM$FX.SYN=as.formula("Y~Day*Stage")
GLMM$RE.SYN=as.formula("~us(Stage):St+us(Stage):Name")
GLMM$RC.SYN=as.formula("~units")

#GLMM=list()
#GLMM$FX.SYN=as.formula("Y~1")
#GLMM$RE.SYN=as.formula("~us(Stage):St+us(Stage):Name")
#GLMM$RC.SYN=as.formula("~units")


F=list()
if (EST=="LMBF") {
  #set the model formula - using lmBF
  F$EF.LMBF=as.formula("Y ~ Stage*Day+Stage*Sf+Stage*Name") #just PEF
  #F$EF.LMBF=as.formula("Y ~ Stage*Sf+Stage*Name") #just PEF
  F$RE.LMBF=c("Name","Sf")
  F$METHOD.LMBF="auto"
  F$IS_SAMPLES.LMBF=1000
  is.VE=FALSE
} else {
  #set the model formula - using lme4/lmer()
  F$FRE.LMER=as.formula("Y~Day*Stage+(-1+Stage|Sf)+(-1+Stage|Name)") #just PEF
  F$MT.LMER=TRUE #if we vary the RE's only, REML seems fine for the BIC
  is.VE=FALSE
}

#prior on non-S parameters collectively psi - respects above structure
PRIOR<-list(R=list(V=1, nu=1), G=list(G1=list(V=diag(3)/4, nu=5),
                                      G2=list(V=diag(3)/4, nu=5)), B=list(mu=rep(0,6),V=1*diag(6)))

#PRIOR<-list(R=list(V=1, nu=1), G=list(G1=list(V=diag(3)/4, nu=5),
#            G2=list(V=diag(3)/4, nu=5)), B=list(mu=rep(0,1),V=1*diag(1)))



#
bd.real=makeSdata(dataset=EXP,synth=FALSE,thin=AAinclude,uniq.cov=UNIQ.COV,cut.cov=COV.REMOVE,keep.cov=COV.KEEP)
str(bd.real)
n=dim(bd.real)[1]
levels(bd.real$Name)<-sub(pattern="A", replacement="a", x=levels(bd.real$Name))

B = bd.real
AAlable = sub(pattern="a", replacement="A", x=levels(bd.real$Name))

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
  sig = 1/rgamma(1,10,1)
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



#fine the approximate credible set for synthetic data y_i
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

#leagle permutation of y corresponds to relabling the data
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


library(doParallel)
getDoParWorkers()
registerDoParallel()
getDoParWorkers()



lmeregression <- function(B,N) {
  F = list()
  F$FRE.LMER=as.formula("Y~Day*Stage+(-1+Stage|Sf)+(-1+Stage|Name)") #just PEF
  F$MT.LMER=TRUE #if we vary the RE's only, REML seems fine for the BIC
  simulateprior = priorfunc(B)
  Bsim = B
  part = prettyPart(mg2gm(simulateprior[[1]],AAlable),AAlable)
  Bsim$Y = simulateprior[[2]]
  partlistform = mg2gm(simulateprior[[1]],levels(Bsim$Name))
  print(part)
  ab = runMCMCmodified(B=Bsim,"LMER",FALSE,F,N,1,N/10,1,FALSE,NA,APR,SAVE_TO_FILE=FALSE,"", PART = partlistform)
  hpd95 = names(hpd(ab$CLf,0.95,NA,AAlable)$pmf)
  indica = is.element(part, hpd95)
  summarystat = unlist(lapply(quadrandseq(Bsim$Y,Bsim),sort))
  return(list(indica, summarystat, Bsim$Y, hpd95, part))
}


lmeregressioncatch <- function(B,N){
  result = tryCatch(lmeregression(B,N), error = function(err){NA})
  return(result)
}


tryregression <- foreach(k.dat = 1:160,.export = ls()) %dopar% {lmeregressioncatch(B,2000)}
outfile.theta=paste(gather.dir,"sim1.RData",sep="")

save(tryregression, file=outfile.theta) 