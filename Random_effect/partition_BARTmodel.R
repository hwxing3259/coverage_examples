library(bartMachine)
source("partitionhelpers.R")

load(".RData")
load(".RData")

#######################################################################
#############Regression with 17 dimensional summary statistics#########
#######################################################################

###construct (c_i, W_i) pairs
regdata17 = partitionreg17(dataregression)
colnames(regdata17)[2:18] = as.character(1:17)
regdata17 = as.data.frame(regdata17)

###observation W_obs
partobs17 = unlist(lapply(quadrandseq(B$Y,B),sort))
partobs17 =  data.frame(t(partobs17))
colnames(partobs17) = as.character(1:17)

###fit the BART model
bartpart17 = bartMachine(X = regdata17[,-1], as.factor(1 - regdata17[,1]))
bartpredpart17 = predict(bartpart17, partobs17, type = "prob")
print(bartpredpart17)
calc_credible_intervals(bartpart17, partobs17)




######same model on randomly chosen subsets of regdata17, checking stability

###partition the synthetic data
partitionindex = createFolds(1:dim(regdata17)[1], k = 4)

###fit 4 submodels
for (i in 1:4){
  regdata17_sub = regdata17[partitionindex[[i]],]
  bartpart17_sub = bartMachine(X = regdata17_sub[,-1], as.factor(1 - regdata17_sub[,1]))
  bartpredpart17_sub = predict(bartpart17_sub, partobs17, type = "prob")
  print(bartpredpart17_sub)
  calc_credible_intervals(bartpart17_sub, partobs17)
}


#######################################################################
#############Regression with 72 dimensional summary statistics#########
#######################################################################

###construct (c_i,T_i) pairs
regdata72 = partitionreg72(dataregression,B)
colnames(regdata72)[2:73] = as.character(1:72)
regdata72 = as.data.frame(regdata72)

###The observation T_obs
partobs72 = aggregate(B$Y, list(B$SSSS, B$DDDD,B$NNNN),mean)[,4]
partobs72 =  data.frame(t(partobs72))
colnames(partobs72) = as.character(1:72)

###fit the bart model
bartpart72 = bartMachine(X = regdata72[,-1], y = as.factor(1 - regdata72[,1]))
predict(bartpart72, new_data = partobs72)
bartpart72$pred_type
calc_credible_intervals(bartpart72, new_data = partobs72)



