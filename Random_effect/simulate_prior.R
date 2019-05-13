
source("setup.R")
#local helper functions
source("auxiliary.R")
source("auxiliary_new.R")
source("partitionhelpers.R")

#####################################################################
################PROCESS THE RAW PPP DATA#############################
######################################################################

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
ZZinclude=c("1","2","3","4","5","6","7","8",
            "9","10","11","12","13","14","15","16","17","18","19")
#just 12
kill=c("total","9","7","11","10","12","19","17")
ZZinclude=setdiff(ZZinclude,kill)

#do we thin the data for unique covariate patterns?
UNIQ.COV=FALSE; COV.REMOVE=c("JJJJ","KKKK","LLLL","MMMM"); COV.KEEP=c("SSSS","DDDD","NNNN")

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
EXP="PPP" #or "HHH" tho that isnt fully implemented in this file it is fine elsewhere

#Model for data given partition St
GLMM=list()
GLMM$FX.SYN=as.formula("Y~DDDD*SSSS")
GLMM$RE.SYN=as.formula("~us(SSSS):St+us(SSSS):NNNN")
GLMM$RC.SYN=as.formula("~units")

#GLMM=list()
#GLMM$FX.SYN=as.formula("Y~1")
#GLMM$RE.SYN=as.formula("~us(SSSS):St+us(SSSS):NNNN")
#GLMM$RC.SYN=as.formula("~units")


F=list()
if (EST=="LMBF") {
  #set the model formula - using lmBF
  F$EF.LMBF=as.formula("Y ~ SSSS*DDDD+SSSS*Sf+SSSS*NNNN") #just PPP
  #F$EF.LMBF=as.formula("Y ~ SSSS*Sf+SSSS*NNNN") #just PPP
  F$RE.LMBF=c("NNNN","Sf")
  F$METHOD.LMBF="auto"
  F$IS_SAMPLES.LMBF=1000
  is.VE=FALSE
} else {
  #set the model formula - using lme4/lmer()
  F$FRE.LMER=as.formula("Y~DDDD*SSSS+(-1+SSSS|Sf)+(-1+SSSS|NNNN)") #just PPP
  F$MT.LMER=TRUE #if we vary the RE's only, REML seems fine for the BIC
  is.VE=FALSE
}

#prior on non-S parameters collectively psi - respects above structure
PRIOR<-list(R=list(V=1, nu=1), G=list(G1=list(V=diag(3)/4, nu=5),
                                      G2=list(V=diag(3)/4, nu=5)), B=list(mu=rep(0,6),V=1*diag(6)))


#
bd.real=makeSdata(dataset=EXP,synth=FALSE,thin=ZZinclude,uniq.cov=UNIQ.COV,cut.cov=COV.REMOVE,keep.cov=COV.KEEP)
str(bd.real)
n=dim(bd.real)[1]
levels(bd.real$NNNN)<-sub(pattern="A", replacement="a", x=levels(bd.real$NNNN))

B = bd.real
ZZlable = sub(pattern="a", replacement="A", x=levels(bd.real$NNNN))



######################################################################################
#############simulating synthetic data from the joint prior##########################
####################################################################################
library(doParallel)
getDoParWorkers()
registerDoParallel()
getDoParWorkers()



##sampling synthetic data from hte joint prior
lmeregression <- function(B,N) {
  F = list()
  F$FRE.LMER=as.formula("Y~DDDD*SSSS+(-1+SSSS|Sf)+(-1+SSSS|NNNN)") #just PPP
  F$MT.LMER=TRUE #if we vary the RE's only, REML seems fine for the BIC
  simulateprior = priorfunc(B)
  Bsim = B
  part = prettyPart(mg2gm(simulateprior[[1]],ZZlable),ZZlable)
  Bsim$Y = simulateprior[[2]]
  partlistform = mg2gm(simulateprior[[1]],levels(Bsim$NNNN))
  print(part)
  ab = runMCMCmodified(B=Bsim,"LMER",FALSE,F,N,1,N/10,1,FALSE,NA,APR,SAVE_TO_FILE=FALSE,"", PART = partlistform)
  hpd95 = names(hpd(ab$CLf,0.95,NA,ZZlable)$pmf)
  indica = is.element(part, hpd95)
  summarystat = unlist(lapply(quadrandseq(Bsim$Y,Bsim),sort))
  return(list(indica, summarystat, Bsim$Y, hpd95, part))
}


#optimx is not always stable, catch errors
lmeregressioncatch <- function(B,N){
  result = tryCatch(lmeregression(B,N), error = function(err){NA})
  return(result)
}



#simulating K synthetic data points
K = 300
tryregression <- foreach(k.dat = 1:K,.export = ls()) %dopar% {lmeregressioncatch(B,2000)}
outfile.theta=paste(gather.dir,"sim1.RData",sep="")

save(tryregression, file=outfile.theta) 
