library(bartMachine)
source("Aminoacidhelpers.R")

load("aminoregression.RData")
load("PEF_subset.RData")

#######################################################################
#############Regression with 17 dimensional summary statistics#########
#######################################################################
regdata17 = partitionreg17(aminoregression)
colnames(regdata17)[2:18] = as.character(1:17)
regdata17 = as.data.frame(regdata17)

partobs17 = unlist(lapply(quadrandseq(B$Y,B),sort))
partobs17 =  data.frame(t(partobs17))
colnames(partobs17) = as.character(1:17)


bartpart17 = bartMachine(X = regdata17[,-1], as.factor(1 - regdata17[,1]))
bartpredpart17 = predict(bartpart17, partobs17, type = "prob")
print(bartpredpart17)
calc_credible_intervals(bartpart17, partobs17)


#######################################################################
#############Regression with 72 dimensional summary statistics#########
#######################################################################


regdata72 = partitionreg72(aminoregression,B)
colnames(regdata72)[2:73] = as.character(1:72)
regdata72 = as.data.frame(regdata72)


partobs72 = aggregate(B$Y, list(B$Stage, B$Day,B$Name),mean)[,4]
partobs72 =  data.frame(t(partobs72))
colnames(partobs72) = as.character(1:72)


bartpart72 = bartMachine(X = regdata72[,-1], y = as.factor(1 - regdata72[,1]))
predict(bartpart72, new_data = partobs72)
bartpart72$pred_type
calc_credible_intervals(bartpart72, new_data = partobs72)

