
##inputs
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



####sampling from joint prior usung MCMC

isingsimprior <- function(N){
  
  rpriorsample = sort(2*runif(N))
  
  #initialize:
  initialising = ising(min(rpriorsample), 100000000,10000,200,"rua",nbrs(200,200))
  startingmtrx = initialising[[2]]
  
  
  isingprior = isingsampling.seq(rpriorsample, startingmtrx)
  
  thetavec = isingprior$phi
  hashvec = isingprior$s
  return(cbind(thetavec, hashvec))
}



##simulating 1000 pairs from the joint prior.
priorising = isingsimprior(1000)

colnames(priorising) = c("theta","hashy")
priorising = as.data.frame(priorising)

save("priorising", file = "priorising.RData")

