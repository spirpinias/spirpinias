## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                Code for making life easy
##                                      Regress every variable using everything else
##                                   Just one of many ways you may want to use a dataset
##                                  Simple implemenation that can be changed in any way shape or form
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 


gbdtFull=function(features,numDepth,numEta,numThread,numRounds,split){
  
  require(xgboost)
  
  results=lapply(1:dim(features)[1],function (j){
    
    splitData=dim(exp)[2]
    
    samples=sample(1:splitData,size=splitData*split)
    
    features=scale(x = features,center = TRUE,scale = TRUE)
    
    train=features[,samples]
    test=features[,-samples]
    
    param=list(max_depth=as.numeric(numDepth),eta=as.numeric(numEta),nthread=as.numeric(numThread),objective="reg:squarederror",gamma=0,min_child_weight=1,max_delta_step=0,lambda=5,tree_method="auto")
    
    dtrain=xgb.DMatrix(as.matrix(t(train[-j,])),label=as.numeric(train[j,]))
    
    dtest=xgb.DMatrix(as.matrix(t(test[-j,])),label=as.numeric(test[j,]))
    
    
    watchList=list(train=dtrain,eval=dtest)
    modelFit=xgb.train(params = param,dtrain,nrounds=as.numeric(numRounds),watchList)
    yHat=predict(modelFit,dtest)
    myFeatures=xgb.importance(model=modelFit)[1:10]
    res=list(MSE=mean((as.numeric(yHat)-as.numeric(test[j,]))^2),myFeatures=data.frame(features = xgb.importance(model=modelFit)[1:10][,1], gain = xgb.importance(model=modelFit)[1:10][,2],cover=xgb.importance(model=modelFit)[1:10][,3],frequency=xgb.importance(model=modelFit)[1:10][,4]))
  })
  return(results)
}
