## This is an algorithm that I wrote in the Laboratory of Bin Zhang. 
## It operates on the ranked gene list of MEGENA prior to embedding the relationships onto a planar graph. 
## Very simple concept. 
## Ranked MEGENA EdgeList -> Directed Graph -> Iterate every node and build local models using neighborhoods.
## This implementatio uses linear models, but I have 2 other variants that are structurally identical except using randomForest and gradient boosted decision trees.

#M_network is the ranked MEGENA edgelist. I have withheld the actual code for simplicity. You get the idea.

network=as.data.frame(M_network)

#Some genes were removed by MEGENA, we no longer have to visit them. 
goTo=sort(unique(union(x=network[,1],y=network[,2])))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(igraph)

##Convert Pairs of Genes to Directed Graph.
graph=graph_from_edgelist(as.matrix(network[,1:2],ncol=2), directed =T)

#Create the degree vector
de=degree(graph = graph,v = goTo)

#If the target less than k degrees for k in R, do not build models. In this case k<3.
goTo=goTo[-which(de<3)]
de=de[-which(de<3)]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
## Objective : In a network of variables, use the neighborhood to predict the Home.
## We use no sparsity and maintain simplicity with a Ridge/EN/Lasso Regression for starters.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## FIRST IMPLENTATION MAY 21st 2021 by Stephen Pirpinias

library(glmnet)
library(caret)
library(purrr)
library(furrr)
library(tictoc)
library(future)

boot_fx=function(j){
  
  ##GET INDICES OF THE NEIGHBORHOOD AROUND YOUR NODE J.
  indices=unlist(ego(graph = graph,order = 3,nodes = goTo[j],mode = c("in"),mindist = 0)[1])[-1]
   
  ##WE HAVE TO SUBSET THE DATASET BEFORE FITTING THE MODEL..
  predictors=t(exp)[,indices]
  target=t(exp)[,goTo[j]]
  
  ##70/30 RANDOM SPLIT
  choice=caret::createDataPartition(y = 1:214,times = 1,p = 0.70)[[1]]
  
  ##PARTITIONING
  training=predictors[-choice,]
  testing=predictors[choice,]
  
  #FIT A RIDGE/EN/LASSO REGRESSION WITH CROSS VALIDATION (10 folds).
  result1=glmnet::cv.glmnet(x = as.matrix(training),y = as.matrix(target[-choice]),alpha=0.5,family="gaussian")
  
  #MAKE PREDICTIONS USING LAMBDA MIN.
  yHat=predict(result1, newx=as.matrix(testing),s=result1$lambda.min,type="response")
  
  #COEFFICIENTS AT LAMBDA MIN.
  myCoefs=coef(result1, s=result1$lambda.min)
  
  #DEPENDS ON HOW YOU WANT TO SEE THE PROBLEM.
  #res=list(myFuture=data.frame(yPred=yHat,yTrue=after[choice]),myCoefs=myCoefs)
  res=list(Gene=row.names(exp)[goTo[j]],MSE=mean(yHat-target[choice])^2,myCoefs=data.frame(name = myCoefs@Dimnames[[1]][which(myCoefs!=0)], coefficient = myCoefs[which(myCoefs!=0)]))
  
  return(res)}  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ANALYSIS  
  
## TESTING USING 4 CORES.
plan(multisession,workers=4)

##USING PURRR POSSIBLY FOR EXCEPTION HANDELING.
##USING FUTURE_MAP to parallelize MAP.  
tic()
possiblySomeFunction=possibly(boot_fx,otherwise = "Less than or Equal to 2 Variables in Neighborhood")
possiblySomeInfo=future_map(1:length(goTo),possiblySomeFunction)
toc()

## Which Indices are not models.
## If Character, the element has No Model.
first=which(lapply(possiblySomeInfo,function(x) is.character(x))==TRUE)

## Create a new list of only the models.
cleanList=possiblySomeInfo[-first]
