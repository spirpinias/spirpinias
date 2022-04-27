#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  Packages
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require(data.table)
require(tibble)
require(purrr)
require(randomForest)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  Set Directory
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dir="/home/spirpinias/Desktop/BinZhang/pn/"
setwd(dir)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                           Import Expression Dataset
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``

exp=fread(file = "msbb.BM_36.CDR_adjusted.PMI_AOD_race_sex_RIN_exonicRate_rRnaRate_batch_adj.tsv") %>% as.data.frame()
exp=column_to_rownames(.data = exp,var = "V1")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  Partition the Data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

kPartition=function(numFolds,numSamples) {
  
  ## Break it up into 4 folds 
  clusters=as.factor(sample(x = 1:as.numeric(numFolds),size = as.numeric(numSamples),replace = TRUE))
  
  ## Grouping
  myGroups=data.table(Index=1:as.numeric(numSamples),Group=clusters)
  
  ## Put Lists into folds.
  myFolds=map(1:as.numeric(numFolds),function(x) 
    myGroups[Group==as.numeric(x)]$Index
    )
    
  ## Sanity Check.
  if (length(reduce(myFolds,intersect))== 0) {
    return(myFolds)
  }
  else{
    print('Folds are Not Unique')
  }

}

## Create 4 equal folds of the dataset.
stratify=kPartition(numFolds = 4,numSamples = 215)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  Fit the Models
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fittingModel=function(stratify,expression){
  
  res=map(stratify,function(a) {
    
    ## 80/20 Splitting
    indices=sample(x = 1:length(a),size = length(a)*0.80)
    
    train=expression[,a[indices]]
    test=expression[,a[-indices]]
    
    modelFit=randomForest(x = as.matrix(t(train[-1,])),y = as.matrix(train[1,]),ntree=100,importance=TRUE)
    
    yHat=predict(modelFit,as.matrix(t(test[-1,])),type="response")
    
    myCoefs=importance(modelFit,scale=TRUE)
    myCoefs=as.data.frame(myCoefs[-which(myCoefs[,1]==0),])
    myCoefs=myCoefs[order(-myCoefs$IncNodePurity),]
    
    list(Gene=row.names(exp)[1],MSE=mean((as.numeric(yHat)-as.numeric(test[1,])^2)),myCoefs=myCoefs,maxImpurity=myCoefs[1:5,])
  })  
    
  return(res)
  
}  
 

testingModel=fittingModel(stratify = stratify,expression = exp)
