## This is a program that I wrote in the Laboratory of Jan Krumsiek for predicting binary phenotype(diabetes) from metabolomics data.
## You can download these 2 files from https://figshare.com/articles/dataset/Qatar_Metabolomics_Study_on_Diabetes/5904022?file=10531345

library(readxl)
library(tidyverse)

#Data 
metabolome=read_excel("/home/spirpinias/Desktop/QMDiab_metabolomics_Preprocessed.xlsx")
phenotypes=read_excel("/home/spirpinias/Desktop/QMDIAB/Urine/QMDiab_phenotypes.xlsx")

#Cleaning
phenotypes=phenotypes %>% remove_rownames() %>% column_to_rownames(var="QMDiab-ID") %>% as.data.frame()
metabolome=metabolome %>% remove_rownames() %>% column_to_rownames(var="QMDiab-ID") %>% as.data.frame() %>% select(-c(AGE,GENDER,BMI,ETH,T2D))

##Implement Linear Model

straightLines=function(predictors,targets,numCores,numFoldsOuter,numFoldsInner,alpha){
  
#Predictors = Independent Variables
#Targets = Dependent Variable 
#numCores = Threads Available
#NumFoldsOuter = How many folds to split the data into.
#NumFoldsInner = How many folds to conduct CV for choice of Lambda.
#Alpha = Sparsity [0,1]
#Family = Binomial (Designed to stay Binomial)
  
  
library(glmnet)
library(parallel)
library(mltools)
library(caret)
  
testingSeq=seq(1,length(targets))
testingFolds=createFolds(y = targets,k = as.numeric(numFoldsOuter),list = TRUE,returnTrain = FALSE)

followMe=1:length(testingFolds)
boot_fx=function(j){
  
    #Predict at LambdaMin
    result1=glmnet::cv.glmnet(x = as.matrix(predictors[testingSeq[-testingFolds[[j]]],]),y = targets[testingSeq[-testingFolds[[j]]]],alpha=alpha,family="binomial",nfolds = as.numeric(numFoldsInner))
    yHat=predict(result1, newx=as.matrix(predictors[testingSeq[testingFolds[[j]]],]),s=result1$lambda.min,type="class")
    myCoefs=coef(result1, s=result1$lambda.min)
    res=list(myFuture=data.frame(yHat=yHat,yTrue=targets[testingSeq[testingFolds[[j]]]],Index=testingSeq[testingFolds[[j]]]),myCoefs=myCoefs)
    return(res)
}

#Models
results=mclapply(followMe,boot_fx,mc.cores=numCores)

#Calculate AUC for Every Fold
aucScore=lapply(results,function(x) mltools::auc_roc(pred=x$myFuture$X1,actual=x$myFuture$yTrue))

#Average the 10 Models into One
AverageModel=Reduce('+',lapply(results,function(x) x$myCoefs))/as.numeric(numFoldsOuter)

#Grab the Indices of Sampling
Indexes=lapply(results, function(x) x$myFuture)

#Send it All Out 
sendOut=list(aucScore,AverageModel,Indexes)

return(sendOut)
}

#At the time it was doing Lasso, but you could change it. (Alpha=1 is Lasso in glmNet)
checkThis=straightLines(predictors=metabolome,targets=phenotypes$T2D,numCores=4,numFoldsOuter=10,numFoldsInner=5,alpha=1)
