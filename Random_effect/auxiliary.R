#some auxiliary functions and optimisation settings for lme/lmer


#gm:list of partition
#M: vector of factors
gm2mg<-function(GM,M) {
  #convert groups listing members to list of members indicating groups
  group.indices=lapply(GM,function(x){which(is.element(as.matrix(M),x))})
  MG=rep(NA,length(M))
  group.sizes=unlist(lapply(group.indices,length))
  MG[unlist(group.indices)]=rep(1:length(GM),group.sizes)
  return(MG)
}

#mg: vector of group labels
#M: vector of factors
mg2gm<-function(MG,M) {
  #convert list of members indicating groups to groups listing members
  GM=lapply(as.list(unique(MG)), function(x) {M[MG==x]})
  return(GM)
}


#DP.sim returns N lables of group assignment
DP.sim<-function(alpha,n) {
  #one person at table 1
  S=1
  K=1
  #bring in the customers 1 by 1
  for (i in 2:n) {
    #n_k=nS[k], table() counts the number at each table k=1..K
    nS=table(S)
    #the probability vector is propto [n_1,...,n_K, alpha]
    pr=c(nS,alpha); pr=pr/sum(pr)
    #pick a table
    k=sample(1:(K+1),1,prob=pr)
    if (k==(K+1)) { #new cluster
      K=K+1
      S=c(S,K)
    } else {
      S=c(S,k)
    }
  }
  return(S)
}


#given a list of partition, returns a string
#xx: list of partition
#M: vector of string of factors
prettyPart<-function(xx,ZZ) {
  #partition to string, pretty format
  ts=mg2gm(gm2mg(xx,ZZ),ZZ); 
  paste(lapply(ts,function(x){sprintf("(%s)",paste(x,collapse=","))}),collapse=",")
}





hpd<-function(CL,lev=0.95,test=NA,ZZ) {
  #calculate hpd set of partitions from an array CL in MG (vector of integer set indices) format
  #also check if given partition "test" (presented in GM/list of sets format) is in the output and report its post prob
  
  #maps numeric label to its order, so get unique representation
  #e.g. label = 4, appear 3rd, then it maps label 4 to 3, as it is the 3rd to appear
  eg1=t(apply(CL,1,function (x) {K=max(x); n=length(x); m=order(unique(x)); for (i in 1:n) {x[i]=m[x[i]]}; return(x)} ));
  eg2=apply(eg1,1,mg2gm,ZZ)
  eg3=sprintf("%s",lapply(eg2,function(xx) {paste(lapply(xx,function(x){sprintf("(%s)",paste(x,collapse=","))}),collapse=",")}))
  eg4=table(eg3)
  eg5=eg4[order(eg4,decreasing=TRUE)]/sum(eg4)
  eg6=names(eg5)
  eg7=cumsum(eg5)
  eg8=cbind(round(eg7[eg7<lev],2),eg6[eg7<lev])
  TEST=NA; plc=c()
  if (is.list(test)) {
    jts=prettyPart(test,ZZ)
    plc=which(eg6==jts); 
    if (length(plc)>0) {
      TEST=(eg7[plc]<lev)
      writeLines(sprintf("\n%s\nis in the level %.2g HPD set. Its estimated posterior probability is %.2g.\n",jts,eg7[plc],eg5[plc]))
    } else {
      TEST=0
      writeLines(sprintf("\n%s\ndoes not appear in the sample output.\n",jts))
    }
  }
  if (length(plc)==0) {plc=length(eg7)}
  #print(sprintf("%s",apply(eg8,1,paste,collapse=' ')))
  return(list(test=TEST,level=eg7[plc],set=eg8,pmf=eg5[eg7<lev]))
}



#log cpr probability
log.pr.crp<-function(S,alpha,n) {
  #log CRP cluster probability
  nS=unlist(lapply(S,length))
  K=length(S)
  pr=K*log(alpha)+sum(lgamma(nS))-lgamma(alpha+n)+lgamma(alpha) #dont need if alpha fixed
  return(pr)
}




#ZZ: vector of fectors
#S: list of partitions
#i: index of factor to change
GenerateCandidate<-function(ZZ,S,i) {

  #generate a new partition by moving element i to a randomly chosen set
  #this may be the same set it came from (TODO fix this waste of time) 
  #or a new set

  a=ZZ[i]            #the ZZ name
  A=gm2mg(S,ZZ)
  g=A[i]             #which group is it in?
  
  Sp=S
  #remove ZZ "a" from group "g"
  Sp[[g]]=setdiff(Sp[[g]],a)
  #if we made an empty group then delete that group frmo the list of groups
  if (length(Sp[[g]])==0) {
    Sp=Sp[-g]
  }
  #pick a new group to drop "a" into
  gp=sample(0:length(Sp),1)
  if (gp==0) {
    #if we choose group 0 then we form a new group (ie group "K+1")
    Sp=c(Sp,list(a))
  } else {
    #otherwise we write ZZ "a" into the chosen (existing) group
    Sp[[gp]]=c(Sp[[gp]],a)
  }

  return(Sp)
}




CooM<-function(CL,ZZ,com=NA,coc.file=NA,sam.file=NA) {
  #compute and visualise co-occurence matrix 
  #p(i,j) = prob i,j in same cluster 

  nl=length(ZZ)  
  if (is.na(com)) {
    Nsamp=dim(CL)[1] #J/SS;
    com=matrix(0,nl,nl,dimnames=list(ZZ,ZZ))
    for (i in 1:nl) {
      for (j in 1:nl) {
        for (m in 1:Nsamp) {
          com[i,j]=com[i,j]+(CL[m,i]==CL[m,j])
        }
      }
    }
    com=com/Nsamp
    com[com==1]<-1-1/Nsamp
    diag(com)<-1
  }

  #display COOM as an image
  d=as.dist(1-com)
  pth=order_tsp(d, method = "nearest", cycle=FALSE,improve=TRUE,path_dir = path_cor)

  if (!is.na(coc.file)) {pdf(coc.file)} else {windows();}
  image(com[pth,pth],xaxt='n',yaxt='n')
  axis(1, at=seq(0,1,length.out=nl),label=rownames(com)[pth],cex.axis=0.7,las=2)
  axis(2, at=seq(0,1,length.out=nl),label=colnames(com)[pth],cex.axis=0.7,las=2)
  if (!is.na(coc.file)) {dev.off()}

  #dispaly MDS projection of COOM - color/label by ZZ properties
  ZZ.dat=read.csv("ZZP.csv",header=TRUE) #load some info about ZZ's in same order as ZZ
  xy=sammon(d,y=matrix(runif(2*nl),nl,2),magic=0.51,trace=FALSE)$points

  if (!is.na(sam.file)) {pdf(sam.file)} else {windows();}
  plot(xy,col=1+(ZZ.dat$MW>145),pch=4*(2+as.numeric(ZZ.dat$Side.chain.polarity)),ann=FALSE)#,
  #     xlim=c(-0.6,0.7),ylim=c(-0.4,0.6))
  text(xy,ZZ,pos=3,col=1+(ZZ.dat$MW>145),cex=0.6)
  if (!is.na(sam.file)) {dev.off()}

  #Principal components if you like
  #windows(); plot(cmd<-cmdscale(d),col=1+(ZZ.dat$MW>145),pch=4*(2+as.numeric(ZZ.dat$Side.chain.polarity)))
  #text(cmd,ZZ,pos=3,col=1+(ZZ.dat$MW>145),cex=0.6)
  #as.character(ZZ.dat[,1]) #full ZZ names

  return(com)
}


#alternative optimisation control for lme
BOB=lmeControl(maxIter = 1000, msMaxIter = 1000, tolerance = 1e-4, niterEM = 25,
           msMaxEval = 400,msTol = 1e-4, msVerbose = FALSE,
           returnObject = TRUE, gradHess = TRUE, apVar = FALSE,
	     .relStep = .Machine$double.eps^(1/3), minAbsParApVar = 0.05,
           opt = "nlminb",optimMethod = "Nelder-Mead",sigma = NULL)

#source(system.file("utils", "allFit.R", package="lme4"))
#  fm1.all <- allFit(b.m)
#  ss <- summary(fm1.all)
#  ss$ fixef               ## extract fixed effects
#  ss$ llik                ## log-likelihoods
#  ss$ sdcor               ## SDs and correlations
#  ss$ theta               ## Cholesky factors
#  ss$ which.OK            ## which fits worked

#alternative optimisation control for lmer (nlminb optim choice helps)
DAVE=lmerControl(optimizer = "optimx",
    restart_edge = TRUE,
    boundary.tol = 1e-5,
    calc.derivs=FALSE,
    use.last.params=FALSE,
    sparseX = FALSE,
    ## input checking options
    check.nobs.vs.rankZ = "ignore",
    check.nobs.vs.nlev = "stop",
    check.nlev.gtreq.5 = "ignore",
    check.nlev.gtr.1 = "stop",
    check.nobs.vs.nRE="stop",
    check.rankX = "warn+drop.cols",
    check.scaleX = "silent.rescale",
    check.formula.LHS = "stop",
    ## convergence checking options
    check.conv.grad     = .makeCC("warning", tol = 2e-3, relTol = NULL),
    check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4),
    check.conv.hess     = .makeCC(action = "warning", tol = 1e-6),
    ## optimizer args
    optCtrl = list(method="nlminb",maxit=1000,iter.max=1000) #list(method="L-BFGS-B",maxit=1000)
  )


GenerateCandidate2 <- function(ZZ,S,i){
  part = GenerateCandidate(ZZ,S,i)
  #print(length(part))
  while (length(part) == 1) {part = GenerateCandidate(ZZ,S,i)}
  return(part)
}