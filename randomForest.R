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

#Implement Random Forest

rfk2Fold=function(predictors,targets,numCores,numTree,numFolds){

#Returns Average Model of the 10 Folds.   
  
library(randomForest)
library(parallel)
library(mltools)
library(caret)
  
testingSeq=seq(1,length(targets))
testingFolds=createFolds(y = targets,k = as.numeric(numFolds),list = TRUE,returnTrain = FALSE)


followMe=1:length(testingFolds)
boot_fx=function(j){
  
    modelFit=randomForest(x = as.matrix(predictors[testingSeq[-testingFolds[[j]]],]),y = as.factor(targets[testingSeq[-testingFolds[[j]]]]),ntree=numTree,importance=TRUE)
    yHat=predict(modelFit,as.matrix(predictors[testingSeq[testingFolds[[j]]],]),type="prob")
    myCoefs=importance(modelFit,scale=TRUE)
    res=list(Future=data.frame(yHat=max.col(m = yHat,ties.method = "random")-1,yTrue=targets[testingSeq[testingFolds[[j]]]],Index=testingSeq[testingFolds[[j]]]),myCoefs=myCoefs)
    
    return(res)
    }

results=mclapply(followMe,boot_fx,mc.cores=numCores)

AverageModel=Reduce('+',lapply(results,function(x) x$myCoefs))/as.numeric(numFolds)
Indexes=lapply(results, function(x) x$Future)

sendOut=list(AverageModel,Indexes)

return(sendOut)
}

CheckThis=rfk2Fold(predictors=metabolome,targets=phenotypes$T2D,numCores=4,numTree=10,numFolds=10)
