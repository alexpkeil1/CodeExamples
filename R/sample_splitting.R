library(qgcomp)
# fake data
alldata = cbind(case=rbinom(100, 1, 0.5), X=matrix(runif(100*10), nrow=100, ncol=10))
head(alldata)

set.seed(100)
datlist = qgcomp:::.split.iid.data(alldata, prop.train=0.6)

trainingdata = datlist$traindata
testdata = datlist$validdata
dim(trainingdata)
dim(testdata)

head(trainingdata)
