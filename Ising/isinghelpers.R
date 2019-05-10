###To find the neighbour of a given vertex
nbrs = function(M,N,type = "free"){
  nb = lapply(1:(M*N), function(i){i})
  if (type =="free"){
    for (i in 1:M){
      for (j in 1:N){
        yN = NULL; yS = NULL;yW = NULL;yE = NULL;
        ind = M*(j-1)+i
        #north neighbour
        if (i == 1){yN = NULL}
        else {yN = ind-1}
        #south neighbour
        if (i == M) {yS = NULL}
        else {yS = ind+1}
        #east neighbour
        if (j == 1) {yE = NULL}
        else {yE = ind-M}
        #west neighbour
        if (j == N) {yW = NULL}
        else {yW = ind+M}
        nb[[ind]] = c(yN,yS,yW,yE)
      }
    }
  }
  else{
    for (i in 1:M){
      for (j in 1:N){
        yN = NULL; yS = NULL;yW = NULL;yE = NULL;
        ind = M*(j-1)+i
        #north neighbour
        if (i == 1){yN = ind+M-1}
        else {yN = ind-1}
        #south neighbour
        if (i == M) {yS = ind-M+1}
        else {yS = ind+1}
        #east neighbour
        if (j == 1) {yE = M*N - M + i}
        else {yE = ind-M}
        #west neighbour
        if (j == N) {yW = i}
        else {yW = ind+M}
        nb[[ind]] = c(yN,yS,yW,yE)
      }
    }
  }
  return(nb)
}


numtoind <- function(M,N,num){
  return(c((num-1)%%M+1, floor((num-1)/M)+1))
}

indtonum <- function(M,N,ind){
  return(M*(ind[2]-1)+ind[1])
}

#get the matrix element whith linear index num
get01 <- function(X,M,N,num){
  indx = numtoind(M,N,num)
  return(X[indx[1],indx[2]])
}


#compute the hash value of a binary image
hashX <- function(X, nb){
  M = dim(X)[1]
  N = dim(X)[2]
  hashxx = 0
  for (i in 1:(M*N)){
    x = get01(X,M,M,i)
    xnb = X[nb[[i]]]
    hashxx = hashxx + sum(x != xnb)
  }
  return(hashxx/2)
}


##sampling from the ising likelihood
ising <- function(theta,L,SS,M,X = "rua", nb){
  ###Note that cylindrical and free boundary
  ###depends on the input nb, careful!
  N = M
  gettd = rep(NA,L/SS)
  if (all(X == "rua")) {X = matrix(round(runif(M*N)),M,N)}
  td = hashX(X, nb)
  for (i in 1:L){
    pix = ceiling(runif(1,0,M*N))
    pixnb = X[nb[[pix]]]
    agree = sum(pixnb == X[pix])
    disagree = sum(pixnb != X[pix])
    if (log(runif(1)) < theta*(disagree - agree)){
      X[pix] = 1 - X[pix]
      td = td + agree - disagree
    }
    if (i%%SS == 0){
      gettd[i/SS] = td
    }
  }
  return(list(gettd, X, nb))
}



#######################Cylindrical approximation##################

##compute the normalizing constant for toroidal likelihood
logZ <- function(m,n,K){
  if (K <= 0.001){
    return(m*n*log(2))
  }
  else{
    c = m*n
    k = 0:(n-1)
    g = c(log(exp(2*K)*tanh(K)), acosh(cosh(2*K)^2/sinh(2*K)-cos(pi*(1:(2*n))/n)))
    Y1=sum(log(2*cosh(m*g[2*k+2]/2)*sinh(2*K)^(m/2)))
    #print(Y1)
    Y2=sum(log(2*sinh(m*g[2*k+2]/2)*sinh(2*K)^(m/2)))
    #print(Y2)
    Y3=sum(log(2*cosh(m*g[2*k+1]/2)*sinh(2*K)^(m/2)))
    #print(Y3)
    ###complex logarithm, log(z) = log(r) + theta*i, theta = pi in this case
    ###only take the maginude r, equivalent to log(abs(-2))
    Y4=sum(log(complex(real = 2*sinh(m*g[2*k+1]/2)*sinh(2*K)^(m/2), imaginary = 0)))
    #print(Y4)
    lZf= Re((c/2-1)*log(2)+Y1+log(1+exp(Y2-Y1)+exp(Y3-Y1)+exp(Y4-Y1)))
    return(lZf)
  }
}



#uniform prior
logprior <- function(theta){
  return(rep(0,length(theta)))
}



#loglikelihood for toroidal ising
loglkd <- function(theta, td, m, n, nume){
  #note that our parameterization is actually theta/2, be consistent with cylindrical paper
  lt = rep(0,length(theta))
  for (i in 1:length(lt)){
    #here we use td - nume, the number of disagreeing particles under cynlidrical conditions,
    #so it is the pdf of the true torodial case
    #but can also try #y, only replace the normalising constant.
    #but we can also try #y, and theta/2 because the original paper double counted the #x
    lt[i] = -theta[i]*(td - nume/2) - logZ(m,n,theta[i]/2)
  }
  return(lt)
}



#Approximated posterior, riemann sum kind of normalization
normpost <- function(theta, td, m, n, free = F){
  if (free){
    nume = 2*m*n - m - n
  }
  else{
    nume = 2*m*n
  }
  lkd = loglkd(theta, td, m, n, nume) + logprior(theta)
  up = exp(lkd - max(lkd))
  np = up/(sum(up)*(theta[2]-theta[1]))
  return(np)  
}



#Approximately sample from approximation
rpost <- function(N, cdf, theta){
  r = runif(N)
  samp = rep(NA,N)
  for (i in 1:N){
    indx = min(which(r[i] < cdf))
    samp[i] = theta[indx]
  }
  return(samp)
}





##################################################################################
###################### Simulation and Calibration ################################
##################################################################################


getCI = function(th, post, alpha){
  upquant = 1 - alpha/2
  lwquant = alpha/2
  postecdf = cumsum(post) * (th[2] - th[1])
  lowerindx = max(which(postecdf < lwquant))+1
  upperindx = min(which(postecdf > upquant))-1
  a = th[lowerindx]
  b = th[upperindx]
  if (is.na(a)) {a = 0} 
  if (is.na(b)) {b = 2} 
  return(c(a, b))
}



isingsampling.seq = function(theta,useX = "rua"){
  theta = sort(theta)
  N = length(theta)
  hashhh = rep(NA,N)
  isingmtrxlist = lapply(1:N, function(i){1})
  if (is.character(useX)){
    useX = matrix(round(runif(200*200)),200,200)
  }
  for (i in 1:N){
    isi = ising(theta[i], 80000000, 80000000, 200, useX, nbrs(200,200))
    hashhh[i] = isi[[1]]
    isingmtrxlist[[i]] = isi[[2]]
    useX = isi[[2]]
  }
  return(list(phi = theta, s = hashhh, mtrx = isingmtrxlist))
}





#######################################################################################
####given the simulated theta and s, compute the calibration based on approximated post.
#######################################################################################


calibalter = function(phi, s){
  N = length(s)
  th = seq(0,2,length.out = 3000)
  LW = UP = cc = ksdist = rep(NA,N)
  postOBSCDF = cumsum(normpost(th, 9270, 200, 200, F) * (th[2] - th[1]))
  loglkdd = loglkd(phi,9270,200,200,80000)
  for (i in 1:N){
    if (i%%100 == 0){print(i)}
    postpdf = normpost(th, s[i], 200, 200, F)
    #postcdf = cumsum(postpdf) * (th[2] - th[1])
    CI = getCI(th, postpdf, 0.05)
    LW[i] = CI[1]
    UP[i] = CI[2]
    cc[i] = ((phi[i] >= CI[1]) & (phi[i] <= CI[2]))
    #ksdist[i] = max(abs(postOBSCDF - postcdf))
  }
  return(list(c = cc, lower = LW, upper = UP, ksdist = ksdist, loglkd = loglkdd))
}




d.prior = function(i){
  if (i <= 0 | i >= 2) {return(0)}
  else(return(0.5))
}

r.prior = function(N){return(as.matrix(runif(N,0,2)))}

model = function(theta,lastsample){
  if (theta < 0 | theta > 2) {return(-200000000)}
  return(ising(theta,500000,5000000,200,lastsample,nbrs(200,200))[[1]])}


dist_metric_ising_ais = function(a,b){
  if (a == -200000000) {return(2)}
  else{
    th = seq(0,2,length.out = 10000)
    postOBSCDF = cumsum(normpost(th, b, 200, 200, F) * (th[2] - th[1]))
    postcdf = cumsum(normpost(th, a, 200, 200, F) * (th[2] - th[1]))
    return(max(abs(postOBSCDF - postcdf)))}
}


########propose a given pair through T_n, which admits f_n as invariant distribution

aisinterpost <- function(theta, hashy, alpha, beta, dist, dprior, originalisingmtrx){
  thetaproposal = runif(1, max(theta - 0.01,0), min(theta + 0.01, 2))
  isingproposal = ising(thetaproposal, 7000000, 7000000, 200, originalisingmtrx, nbrs(200,200))
  hashyproposal = isingproposal[[1]]
  isingmtrxproposal = isingproposal[[2]]
  a = beta*(dist(hashy , 9270) - dist(hashyproposal, 9270)) + 
    (1 - alpha) * (loglkd(thetaproposal, 9270, 200, 200, 80000) - loglkd(theta, 9270, 200, 200, 80000)) + 
    log(dprior(thetaproposal)) - log(dprior(theta))
  print(a)
  if (log(runif(1)) < a) {theta = thetaproposal; hashy = hashyproposal; originalisingmtrx = isingmtrxproposal}
  return(list(theta, hashy, originalisingmtrx))
}


###propogate one pair through T_1,...,T_N-1 M-H kernels

aisonepost <- function(initialtheta,initialhashy, alphaseq, betaseq, dist, dprior, initialisingmtrx){
  thetaseq = rep(NA, length(betaseq))
  hashyseq = rep(NA, length(betaseq))
  logweightseq = rep(NA, length(betaseq))
  cumwtd = rep(NA, length(betaseq))
  theta = initialtheta
  hashy = initialhashy
  thetaseq[1] = theta
  hashyseq[1] = hashy
  logweightseq[1] = -betaseq[1] * dist(hashy, 9270) - alphaseq[1] * loglkd(theta, hashy, 200, 200, 80000)
  cumwtd[1] = logweightseq[1]
  seqmtrx = initialisingmtrx
  for (i in 2:(length(betaseq))){
    prop = aisinterpost(theta, hashy, alphaseq[i-1], betaseq[i-1], dist, dprior,seqmtrx)
    thetaseq[i] = prop[[1]]
    hashyseq[i] = prop[[2]]
    seqmtrx = prop[[3]]
    logweightseq[i] = (betaseq[i-1] - betaseq[i]) * dist(prop[[2]], 9270) + (alphaseq[i-1] - alphaseq[i]) * loglkd(prop[[1]], 9270, 200, 200, 80000)
    cumwtd[i] = cumwtd[i-1] + logweightseq[i]
    theta = prop[[1]]
    hashy = prop[[2]]
  }
  return(t(cbind(theta = thetaseq, hashy= hashyseq, logweight = logweightseq, wtd = cumwtd)))
}














#### wrapper, repeat the procedure for N points using foreach package

AISparpost <- function(thetavec, hashvec, alphaseq, betaseq, dist, dprior){
  #thetavec: is length Nvector from prior
  #hashy: from likelihood with thetavec as parameter
  #betaseq: sequence of tempreture
  #dist: distance metric
  #dprior: prior density
  M = length(betaseq)
  N = length(thetavec)
  raw = foreach(i=1:N, .combine = rbind, .export = c("loglkd", "logZ","dist_metric_ising_ais","aisonepost", "aisinterpost", "ising", "nbrs", "hashX", "dist", "dprior", "thetavec", "hashvec", "betaseq", "alphaseq")) %dopar% aisonepost(thetavec[i],hashvec[i],alphaseq, betaseq, dist, dprior)
  finpar = raw[seq(1,4*N,4),M]
  logwtd = rowSums(raw[seq(3,N*4,4),])
  logwtd = logwtd - max(logwtd)
  finwtd = exp(logwtd)/sum(exp(logwtd))
  storeinter = lapply(1:M, function(i){return(NA)})
  for (j in 1:M){
    interpar = raw[seq(1,4*N,4),j]
    lwtd = raw[seq(4,N*4,4),j]
    lwtd = lwtd - max(lwtd)
    fwtd = exp(lwtd)/sum(exp(lwtd))
    storeinter[[j]] = cbind(par = interpar, wtd = fwtd)
  }
  
  return(list(par = finpar, weight = finwtd, inter = storeinter, raw = raw))
}



#handeling raw output of AIS

AISisingprocess <- function(raw, betalength, thetalength){
  M = betalength
  N = thetalength
  finpar = raw[seq(1,4*N,4),M]
  logwtd = rowSums(raw[seq(3,N*4,4),])
  logwtd = logwtd - max(logwtd)
  finwtd = exp(logwtd)/sum(exp(logwtd))
  storeinter = lapply(1:M, function(i){return(NA)})
  for (j in 1:M){
    interpar = raw[seq(1,4*N,4),j]
    lwtd = raw[seq(4,N*4,4),j]
    lwtd = lwtd - max(lwtd)
    fwtd = exp(lwtd)/sum(exp(lwtd))
    storeinter[[j]] = cbind(par = interpar, wtd = fwtd)
  }
  
  return(list(par = finpar, weight = finwtd, inter = storeinter, raw = raw))
}