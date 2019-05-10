

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



##simulating y from prior given partition (list form)

priorysimulate <- function(B, listpart){
  aalable = levels(B$Name)
  #listpart = mg2gm(part, aalable)
  B$Sf = as.factor(gm2mg(listpart, B$Name))
  b.m=lmer(formula = Y ~ Day * Stage + (-1 + Stage | Sf) + (-1 + Stage | Name),
           data = B, REML = T, control = DAVE)
  modelmtrx = cbind(getME(b.m, "X"), as.matrix(getME(b.m, "Z")))
  npart = length(listpart)
  #fixed effect: normal(0,1) for all 6 variables
  beta = rnorm(6)
  G1cov = rIW(V = diag(0.25,3), nu = 5)
  G2cov = rIW(V = diag(0.25,3), nu = 5)
  G1effect = rmvnorm(12, rep(0,3), G1cov)
  G2effect = rmvnorm(npart, rep(0,3), G2cov)
  finalbeta = c(beta, as.vector(t(G1effect)), as.vector(t(G2effect)))
  #res = 1/rgamma(dim(B)[1],0.5,0.5)
  
  
  ##############CHANGE: lower the sd
  res = rnorm(dim(B)[1], mean = 0, sd = 0.1)
  #############################
  
  priory = modelmtrx %*% finalbeta + res
  return(priory)
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




#input: 
# Part: a vector partition
#y: simulated dataset
#alpha,beta; tempreture
#dist: dist.measure between y and y_simulated
#dprior: crp log prob

aisinterpartition <- function(B, initialpart, initialy, alpha, beta, dist, yobs, partpostmass, anchor, K){
  #generating new partition
  old.y = initialy
  #number of labels
  N = length(levels(B$Name))
  #group lable, lowercase a
  aalabel = levels(B$Name)  
  #group lable, uppercase a
  capitalAlable = unlist(lapply(aalabel, function(i){return(sub(pattern = "a", replacement = "A", x = i))}))
  stringpartition = prettyPart(mg2gm(initialpart, capitalAlable), capitalAlable)
  #BIC for observed data and given partition
  logp = logpr.joint(yobs, stringpartition, B)
  priorold = logp[1]
  approxpostold = logp[2]
  distold = dist(old.y, partpostmass, B, anchor, K)
  #initial partition is a numeric vector
  Sf = mg2gm(initialpart, aalabel)
  updateindx = sample(1:N, 1)
  for (i in updateindx){
    #propose a partition
    Spf=GenerateCandidate2(aalabel,Sf,i)
    #propose a observation vector
    new.y = priorysimulate(B, Spf)  
    stringpartitionnew = prettyPart(mg2gm(gm2mg(Spf, aalabel),capitalAlable), capitalAlable)
    logpnew = logpr.joint(yobs, stringpartitionnew, B)
    priornew = logpnew[1]
    approxpostnew = logpnew[2]
    distnew = dist(new.y, partpostmass, B, anchor, K)
    #print(c((1 - alpha) * (approxpostnew - approxpostold), beta * (distold - distnew)))
    MHR = priornew - priorold + (1 - alpha) * (approxpostnew - approxpostold) + beta * (distold - distnew)
    #MHR = n.crp - o.crp #CRP prior target for debugging
    
    if (log(runif(1))<MHR) {
      Sf = Spf
      old.y = new.y
      priorold = priornew
      approxpostold = approxpostnew
      distold = distnew
    }}
  return(list(gm2mg(Sf, aalabel), old.y))}



#partmass is pre-computed value of observed data
cheap_part_ks <- function(ysim, partmass, B, anchor,K){
  aa = logpr.joint(ysim, anchor[1:K], B)
  aa = aa[,1] + aa[,2]
  aa = aa - max(aa)
  aa = cumsum(exp(aa)/sum(exp(aa)))
  ks = max(abs(partmass - aa))
  return(ks)
}


centered_y_ks <- function(ysim, yobs, dummy2, dummy3, dummy4){
  ysim = ysim - mean(ysim)
  yobs = yobs - mean(yobs)
  ksdist = ks.test(ysim, yobs)$statistic
  return(ksdist)
}



aisonepost <- function(B, initialpart, alphaseq, betaseq, dist, yobs, partpostmass, anchor, K){
  initialy = priorysimulate(B, mg2gm(initialpart, levels(B$Name))) 
  partseq = matrix(NA, ncol = length(levels(B$Name)), nrow = length(betaseq))
  yseq = matrix(NA, ncol = length(yobs), nrow = length(betaseq))
  logweightseq = rep(NA, length(betaseq))
  cumwtd = rep(NA, length(betaseq))
  part = initialpart
  y = initialy
  partseq[1,] = part
  yseq[1,] = y
  aalabel = levels(B$Name)  
  capitalAlable = unlist(lapply(aalabel, function(i){return(sub(pattern = "a", replacement = "A", x = i))}))
  #find the string partition of initial part
  strpart = prettyPart(mg2gm(part,capitalAlable), capitalAlable)
  print(strpart)
  initial.dist = dist(y, partpostmass, B, anchor, K)
  logweightseq[1] = -betaseq[1] * initial.dist - alphaseq[1] * logpr.joint(yobs, strpart, B)[2]
  cumwtd[1] = logweightseq[1]
  for (i in 2:(length(betaseq))){
    prop = aisinterpartition(B, part, y, alphaseq[i-1],  betaseq[i-1], dist, B$Y, partpostmass, anchor,K)
    partseq[i,] = prop[[1]]
    yseq[i,] = prop[[2]]
    strpart = prettyPart(mg2gm(prop[[1]],capitalAlable), capitalAlable)
    #print(strpart)
    newdistance = dist(prop[[2]], partpostmass, B, anchor, K)
    logweightseq[i] = (betaseq[i-1] - betaseq[i]) * newdistance + (alphaseq[i-1] - alphaseq[i]) * logpr.joint(yobs, strpart, B)[2]
    #print(c((betaseq[i-1] - betaseq[i]) * newdistance, (alphaseq[i-1] - alphaseq[i]) * logpr.joint(yobs, strpart, B)[2]))
    cumwtd[i] = cumwtd[i-1] + logweightseq[i]
    part = prop[[1]]
    y = prop[[2]]
  }
  #return(list(part = partseq, y= yseq, logweight = logweightseq, wtd = cumwtd,
  #            finalpart = part, finaly = y, finalwtd = cumwtd[length(betaseq)]))
  #return(c(part,cumwtd[length(betaseq)]))
  return(list(part,cumwtd,partseq,yseq))
}



AISpartpost <- function(B, partition, alphaseq, betaseq, dist, yobs, partpostmass, anchor, K){
  #thetavec: is length Nvector from prior
  #hashy: from likelihood with thetavec as parameter
  #betaseq: sequence of tempreture
  #dist: distance metric
  #dprior: prior density
  M = length(betaseq)
  N = dim(partition)[1]
  raw = foreach(i=1:N, .export = ls()) %dopar% aisonepost(B, partition[i,], alphaseq, betaseq, dist, yobs, partpostmass, anchor, K)
  return(raw)}




wtdpart <- function(rua,indxx = NULL){

  N = dim(rua)[1]
  if(is.null(indxx)) {indxx = dim(rua)[2]}
  wtd = rua[,indxx]
  wtdnorm = exp(wtd - max(wtd))/sum(exp(wtd - max(wtd)))
  ESS = 1/sum(wtdnorm^2)
  AAlable = c("AIA", "GLU", "GLY", "ILE", "LEU", "LYS", "MET", "PHE", "THR", "TRP", "TYR" ,"VAL")
  nameee = rep(NA,N)
  for (i in 1:N){
    nameee[i] = prettyPart(mg2gm(rua[i,1:12], AAlable), AAlable)
  }
  indx = rep(NA,N)
  for (i in 1:N){
    indx[i] = is.element(nameee[i], anchor)
  }
  return(c(sum(indx * wtdnorm), ESS))
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
