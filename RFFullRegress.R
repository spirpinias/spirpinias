## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                Code for making life easy
##                                      Regress every variable using everything else
##                                   Just one of many ways you may want to use a dataset
##                                  Simple implemenation that can be changed in any way shape or form
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

rfkFull=function(features,numTrees,split){
  
  require(randomForest)

  results=lapply(1:dim(features)[1],function (j){
    
    splitData=dim(exp)[2]
    
    samples=sample(1:splitData,size=splitData*split)
    
    features=scale(x = features,center = TRUE,scale = TRUE)
    
    train=features[,samples]
    test=features[,-samples]
    
    modelFit=randomForest(x = as.matrix(t(train[-j,])),y = as.matrix(train[j,]),ntree=numTrees,importance=TRUE)
    
    yHat=predict(modelFit,as.matrix(t(test[-j,])),type="response")
    
    myCoefs=importance(modelFit,scale=TRUE)
    myCoefs=as.data.frame(myCoefs[-which(myCoefs[,1]==0),])
    myCoefs=myCoefs[order(-myCoefs$IncNodePurity),]
    
    res=list(MSE=mean((yHat-test[j,])^2),myCoefs=myCoefs)
  })
  return(results)
}
