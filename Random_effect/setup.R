
library(nlme) #lme
library(lme4) #lmer() (exploreBIC.R), simulate() (synth.R)
library(MASS) #for eg sammom MDS in the mcmc visualisation
library(lattice)
library(car)    #scatter3d 
library(arm)    #maybe the lmer() output MCMC routine
library(BayesFactor) #for its BF estimator if not using BIC
library(optimx) #for the lmer optimisation in auxiliary.R
library(LearnBayes) #I forget, maybe some of the MCMC plotting
library(coda)      #as.mcmc etc
library(PairViz)  #for the TSP ording of co-occurence matrices for viewing
library(MCMCglmm)
library(foreach)
library(mvtnorm)


