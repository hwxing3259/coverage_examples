library(EBImage)
library(Matrix)
library(ggplot2)
library(Hmisc)
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
#############################################################
#to only reproduce the plots#################################
load("isingfiguredata.RData")

#AIS
#plot the estimated c(y) for each iteration of AIS
p1= ggplot(isingaismtrx) + geom_point(aes(x = 1:60,y = wtdmean)) + geom_ribbon(aes(x = 1:60, ymin = CIlo, ymax = CIup),alpha = 0.2) + 
  theme_classic() + geom_hline(aes(yintercept = c1)) + geom_hline(aes(yintercept = c1 - 2*c1sd),linetype = "dotted")+ xlim(1.9,60) + 
  geom_hline(aes(yintercept = c1+2*c1sd),linetype = "dotted") + ylab(paste("Estimated coverage at observed data")) + xlab("Number of iteration")

print(p1)

#GAM, fitteed curve
p2 = ggplot(data = gamcurve) + geom_line(aes(x = s, y = chat)) + ylim(0.4,1) + 
  labs(x = "Sufficient Statistic S(y)", y = "Estimated coverage") + 
  theme_classic() + 
  geom_segment(aes(x = hashobs, xend = hashobs, y = twosigerrorbar[1],yend = twosigerrorbar[2]), linetype = 2) + 
  geom_point(aes(x = hashobs, y = cy))

print(p2)

#results from IS estimator in Lee et.al 2018
rho = (1:10)/10
issampler(rho)
########################################################################
########################################################################






#the CDF is computed using Reimann sum, since the density can be evaluated pointwisely.
th = seq(0,2,length.out = 5000)
postOBSPDF = normpost(th,hashobs,200,200,FALSE)
postOBSCDF = cumsum(postOBSPDF)*(th[2] - th[1])
#compute the approximated CI
approxCI = getCI(th,postOBSPDF,0.05)


###############################################################
###########approximately exact exchange algorithm##############
###############################################################

isingexchange <- function(L,LSS,nbr,X){
  thetavec = rep(NA,L/LSS)
  theta = 0.91 #initial theta
  for (i in 2:L){
    proposal = rnorm(1,mean = theta, sd = 0.007) #random walk proposal
    isingsimulateX = ising(proposal, 30000000,30000000,200,X, nbr)$X
    z = hashX(isingsimulateX, nbr)
    alpha = (z - 9270) * (proposal - theta)
    if (log(runif(1)) < alpha){
      theta = proposal
      X = isingsimulateX
    }
    if (i%%LSS == 0){
      thetavec[i/LSS] = theta
    }
  }
  return(thetavec)
}

#samples from the exchange algorithm
exchangepost = read.csv("exchangepost.csv")$x

#simulating from the exact posterior:
#exchangepost = isingexchange(7500,5,nbrs200,binarymtrx)
acf(exchangepost)
plot(exchangepost, type = "l") #looks promising
#write.csv(exchangepost, file = "exchangepost.csv")


#treat c1 as the ground truth
c1 = mean((exchangepost>approxCI[1]) * (exchangepost<approxCI[2]))
#sd if c1
c1sd = sqrt(var((exchangepost>approxCI[1]) * (exchangepost<approxCI[2]))/1500)

###########################################################
###########Annealed importance sampling####################
###########################################################
load("AISising.RData")

#get weighted mean of ci using ALgorithm2
ci = (AISisingresult$par > approxCI[1]) * (AISisingresult$par < approxCI[2])
ci = as.numeric(ci)
cy = sum(ci * AISisingresult$weight)

ESS = 1/sum(AISisingresult$weight^2)

cysd = sqrt(wtd.var(ci, AISisingresult$weight, normwt = T)/ESS)

twosigerrorbar  = c(cy - 2*cysd, cy + 2*cysd)


#wtdci computes the weighted mean and two-sigma error bar for each iteration of AIS
wtdci <- function(isingresult){
  th = seq(0,2,length.out = 5000)
  postpdf = normpost(th,9270,200,200,FALSE)
  postcdf = cumsum(postpdf)*(th[2] - th[1])
  approxCI = getCI(th,postpdf,0.05)
  opt = matrix(NA,ncol = 3,nrow = 60)
  for (i in 1:60){
    ci = (isingresult$inter[[i]][,1] > approxCI[1]) * (isingresult$inter[[i]][,1] < approxCI[2])
    ci = as.numeric(ci)
    cy = sum(ci * isingresult$inter[[i]][,2])
    ESS = 1/sum(isingresult$inter[[i]][,2]^2)
    cysd = sqrt(wtd.var(ci, isingresult$inter[[i]][,2], normwt = T)/ESS)
    appCI = c(cy - 2*cysd, cy + 2*cysd)
    opt[i,] = c(cy,appCI)
  }
  opt = data.frame(opt)
  colnames(opt) = c("wtdmean","CIlo","CIup")
  return(opt)
}


isingaismtrx = wtdci(AISisingresult)

#plot the estimated c(y) for each iteration of AIS
p1= ggplot(isingaismtrx) + geom_point(aes(x = 1:60,y = wtdmean)) + geom_ribbon(aes(x = 1:60, ymin = CIlo, ymax = CIup),alpha = 0.2) + 
  theme_classic() + geom_hline(aes(yintercept = c1)) + geom_hline(aes(yintercept = c1 - 2*c1sd),linetype = "dotted")+ xlim(1.9,60) + 
  geom_hline(aes(yintercept = c1+2*c1sd),linetype = "dotted") + ylab(paste("Estimated coverage at observed data")) + xlab("Number of iteration")

print(p1)



##################################################################
##########GAM#####################################################
##################################################################
library(mgcv) #version  mgcv 1.8-17
load("priorising.RData")


#compute c_i for each s_i, takes roughly 5 min
isingpriorci = calibalter(priorising$theta, priorising$hashy)
c01 = isingpriorci$c

#fit a GAM model
cov.gam = gam(y~s(x),family=binomial,data=as.data.frame(cbind(x = priorising$hashy, y = c01)))
lp=predict(cov.gam)
print(predict(cov.gam, data.frame(x = 9270), type = "response", se.fit = T))
plot(sort(priorising$hashy),(exp(lp)/(1+exp(lp)))[order(priorising$hashy)],xlab="Sufficient Statistic S(y';E_F)",
     ylab="Estimated Coverage", type = "l", ylim = c(0.4,1), col = 1)


gamcurve = data.frame(s = priorising$hashy, chat = (exp(lp)/(1+exp(lp))))

p2 = ggplot(data = gamcurve) + geom_line(aes(x = s, y = chat)) + ylim(0.4,1) + 
  labs(x = "Sufficient Statistic S(y)", y = "Estimated coverage") + 
  theme_classic() + 
  geom_segment(aes(x = hashobs, xend = hashobs, y = twosigerrorbar[1],yend = twosigerrorbar[2]), linetype = 2) + 
  geom_point(aes(x = hashobs, y = cy))

print(p2)

library(gridExtra)
grid.arrange(p1,p2,ncol = 2)



##################################################################
#########IS sampler from Lee et al (2018)#########################
##################################################################

#(phi,s) pairs sampled from approx.post distribution and exact likelihood
hist(approxpostpair[,1], xlab = "appeoximated posterior")

#compute the ks-distance between s_i and s_obs, takes roughly 3 mins
ksdist = rep(NA,1000)
for (i in 1:1000){
  if (i %% 100 == 0) {print(i)}
  postcdf = cumsum(normpost(th, approxpostpair$hashy[i], 200, 200, FALSE) * (th[2] - th[1]))
  ksdist[i] = max(abs(postOBSCDF - postcdf))
}

#compute tilde(p) and ci for pairs from approx.post. for the IS sampler, takes 5 mins
lkd = loglkd(approxpostpair$phi,hashobs,200,200,80000)
F = FALSE #just in case
approxpostci = calibalter(approxpostpair$phi, approxpostpair$hashy)$c


#issampler compute the IS estimate of c(y) described in Lee et al (2018)
#for given (vector of) tolerance levels rho

issampler <- function(rho){
  returnmtrx = matrix(NA,nrow = length(rho), ncol = 2)
  colnames(returnmtrx) = c("chat","ESS")
  for (i in 1:length(rho)){
    r = rho[i]
    selectindx = which(ksdist < r)
    ciselected = as.numeric(approxpostci[selectindx])
    lkdselected = lkd[selectindx]
    ptilde = exp(lkdselected - max(lkdselected))
    w = 1/ptilde
    w = w/(sum(w))
    chat = sum(ciselected * w)
    ESS = 1/sum(w^2) 
    returnmtrx[i,] = c(chat, ESS)
  }
  returnmtrx = cbind(rho,returnmtrx)
  colnames(returnmtrx)[1] = "rho"
  return(returnmtrx)
}

##output of the IS sampler for different tolerance level rho = 0.1,0.2,...,1. 
##note that ESS is low
rho = (1:10)/10
issampler(rho)

