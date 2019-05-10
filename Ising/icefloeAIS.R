ttt = proc.time()
library(EBImage)
library(Matrix)
img = readImage("floe.png")
img = channel(img, "grey")
img@.Data = img@.Data[-c(1:3),]
img = resize(img,200,200)

img@.Data = round(img@.Data)
image(img)
binarymtrx = img@.Data
image(t(1-Matrix(binarymtrx))[200:1,])
source("isinghelpers.R")

nbrs200=nbrs(200,200)
hashobs = hashX(t(1-binarymtrx),nbrs200)

library(Hmisc)
library(foreach)
library(doParallel)
registerDoParallel()
getDoParWorkers()



#need: 
#post par vector thetavec
#post hashy vector hashvec
#tempreture seequence alphaseq,betaseq
#list of ising matrix correspoding to each post hashy value




####compute the observed CDF using normpost
th = seq(0,2,length.out = 1000)
postOBSCDF = cumsum(normpost(th,hashobs,200,200,F))*(th[2] - th[1])


N = 1000

rpostsample = sort(rpost(N, postOBSCDF, th))

#initialize:
initialising = ising(min(rpostsample), 20000000,10000,200,binarymtrx,nbrs(200,200))
startingmtrx = initialising[[2]]

#sample s(y_i) from samples of approx.post
ising_approxpost = isingsampling.seq(rpostsample, startingmtrx)

thetavec = ising_approxpost$phi
hashvec = ising_approxpost$s
isingmtrxlist = ising_approxpost$mtrx

approxpostpair = data.frame(phi = thetavec, hashy = hashvec)

#tempreture sequence
betaseq = 1.05^c(1:60)
alphaseq = c(seq(0,1,length.out = 45),rep(1,15))

#run AIS sampler
ais1000raw = foreach(i=1:N, .combine = rbind, .export = ls()) %dopar% aisonepost(thetavec[i],hashvec[i],alphaseq, betaseq, dist_metric_ising_ais, d.prior,isingmtrxlist[[i]])

AISisingresult = AISisingprocess(aistry1000,60,1000)


save("approxpostpair", "AISisingresult", file = "AISising.RData")

