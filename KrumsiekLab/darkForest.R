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

#Implement Gradient Boosted Decision Trees

dark2Foldest=function(predictors,targets,numDepth,numEta,numThread,numRounds,numFolds){

## IMPLEMENTS A kFold GRADIENT BOOSTED DECISION TREES FOR FEATURE EXTRACTION
  
## NOTES FOR STEPHEN  

#General Parameters
# Booster is which one to use - gbTree is default
# Verbosity is chatter default=1 for Errors, 0 for No Talk
# nThread is numCores

  
# Parameters for Tree Booster
# ETA is a learning rate [0,1]
# Gamma Conservative Leaf Expansion [0,infinity] -- Higher for More Conservative Rules to Expand Leaves
# Max Depth is the max depth of the Trees [0,infinity]
# Min Child Weight 
  
#Many More Options this function is simple

library(xgboost)
library(mltools)
library(ppcor)
library(caret)  
library(purrr)  

  
testingSeq=seq(1,length(targets))
testingFolds=createFolds(y = targets,k = as.numeric(numFolds),list = TRUE,returnTrain = FALSE)
  
  
theList=lapply(1:length(testingFolds),function(x) {
  
  param=list(max_depth=as.numeric(numDepth),eta=as.numeric(numEta),nthread=as.numeric(numThread),objective="binary:logistic",gamma=0,min_child_weight=1,max_delta_step=0,lambda=5,tree_method="auto")
  dtrain=xgb.DMatrix(as.matrix(predictors[testingSeq[-testingFolds[[x]]],]),label=as.numeric(targets[testingSeq[-testingFolds[[x]]]]))
  dtest=xgb.DMatrix(as.matrix(predictors[testingSeq[testingFolds[[x]]],]),label=as.numeric(targets[testingSeq[testingFolds[[x]]]]))
  watchList=list(train=dtrain,eval=dtest)
  modelFit=xgb.train(params = param,dtrain,nrounds=as.numeric(numRounds),watchList)
  yHat=predict(modelFit,dtest)
  myFeatures=xgb.importance(model=modelFit)[1:10] #Extract top 10 Features.
  res=list(yHat=data.frame(yHat=as.numeric(yHat>0.50),yTrue=as.numeric(targets[testingSeq[testingFolds[[x]]]]),Index=testingSeq[testingFolds[[x]]]),myFeatures=data.frame(features = xgb.importance(model=modelFit)[1:10][,1], gain = xgb.importance(model=modelFit)[1:10][,2],cover=xgb.importance(model=modelFit)[1:10][,3],frequency=xgb.importance(model=modelFit)[1:10][,4]))

}

)


return(theList)
}

secondModel=dark2Foldest(predictors = metabolome,targets = phenotypes$T2D,numDepth = 10,numEta = 0.1,numThread = 4,numRounds = 2,numFolds = 10)

#Check to see if it predicts well using AUC.
mean(unlist(lapply(secondModel,function(x) mltools::auc_roc(preds = x$yHat$yHat,actuals = x$yHat$yTrue))))

#I was investigating which predictions failed and return indices to check the phenotypes for an explanation..
pmap(.l = list(lapply(secondModel,function(x) x$yHat$yHat),lapply(secondModel,function(x) x$yHat$yTrue),lapply(secondModel,function(x) x$yHat$Index)),function(x,y,z) z[which(x!=y)])

#From the phenotypes the samples incorrectly predicted had low BMI. 
View(phenotypes[sort(mySeq[unlist(pmap(.l = list(lapply(secondModel,function(x) x$yHat$yHat),lapply(secondModel,function(x) x$yHat$yTrue),lapply(secondModel,function(x) x$yHat$Index)),function(x,y,z) z[which(x!=y)]))],decreasing = TRUE),])

#People with low BMI were harder to predict having diabetes. Which made sense... 
                                                                                                                                               
