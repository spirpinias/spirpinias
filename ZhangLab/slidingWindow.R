require(purrr)
require(glmnet)
require(igraph)
require(dequer)
require(reshape2)
require(ggplot2)
require(dplyr)

## Set Directory.
dir="/home/spirpinias/Desktop/BinZhang/pn/"
setwd(dir)

## Import expression data.
exp=read.table("msbb.BM_36.CDR_adjusted.PMI_AOD_race_sex_RIN_exonicRate_rRnaRate_batch_adj.tsv",sep = '\t', header = TRUE)
exp=exp[,-1]
gene_names=row.names(exp)

## Import Meta Data
ATP6V1Ameta=read.table(file = '/home/spirpinias/Desktop/BinZhang/ANN/Rcode/BM36/msbb.meta.BM_36.tsv', sep = '\t', header = TRUE)
ATP6V1Ameta=ATP6V1Ameta[-1,]

## Recoding Diagnosis as Numerical
ATP6V1Ameta$Dx.by.braak.cerad[which(ATP6V1Ameta$Dx.by.braak.cerad=="Normal")]=0
ATP6V1Ameta$Dx.by.braak.cerad[which(ATP6V1Ameta$Dx.by.braak.cerad=="ADpp")]=1
labels=as.numeric(ATP6V1Ameta$Dx.by.braak.cerad)

## Perturbation - Validation
ATP6V1A=read.table(file = '/home/spirpinias/Downloads/ATP6V1A.KD.tsv', sep = '\t', header = TRUE)
ATP6V1A=ATP6V1A %>% filter(Contrast=="KD.vs.WT.in.unchallenged")


## Removing Labels which Diagnosis is NA.
features=t(exp)[-which(is.na(labels)==TRUE),]
labels=labels[-which(is.na(labels)==TRUE)]

## Import Keydrivers.
MEGENA_keydrivers=read.delim("all_BN_keydrivers_CTD_subtypes.txt",header = T)
MEGENA_keydrivers=MEGENA_keydrivers$keydrivers
MEGENA_keydrivers=as.numeric(na.omit(match(x = MEGENA_keydrivers,table = row.names(exp))))
MEGENA_keydrivers=unique(MEGENA_keydrivers)

## Import Network
MEGENA_network=read.delim("BM36_CDR_adj_himem.onlyclust.megena_out.tsv",header = T)
M_network=cbind(match(MEGENA_network$row,gene_names),match(MEGENA_network$col,gene_names),MEGENA_network$weight)

## Remove NA from MEGENA network and M_network
MEGENA_network=MEGENA_network[-which(is.na(M_network),arr.ind = T)[,1],]
M_network=M_network[-which(is.na(M_network),arr.ind = T)[,1],]

## Assign Indices of the Edge Information to a Dataframe.
network=as.data.frame(M_network)

## Some genes were removed by MEGENA, we no longer have to visit them. 
## We have to iterate this graph from a "start" to "finish".
goTo=as.numeric(sort(union(x=network[,1],y=network[,2])))

## Remove the Keydrivers from the set of Genes to Iterate.
goTo=setdiff(x = goTo,y = MEGENA_keydrivers)

## Remove the original expression and meta no longer needed. Refer to Features and Labels.
rm(list = c("exp"))

## Network to iGraph Object.
graph=graph_from_edgelist(as.matrix(network[,1:2],ncol=2), directed=T)
E(graph)$weight=network[,3]

## Prep the queue for appending.
keydrivers=as.queue(as.list(MEGENA_keydrivers))

## Estimate node degree average to create a decision point on the next step of what to keep.
checkAVG=floor(mean(as.numeric(degree(graph = graph))))

## Iterate the Graph and Find the Neighbors with any connectivity.
#  check=lapply(X = 1:length(goTo),FUN = function(j) {
#   Neighbors=unlist(ego(graph = graph,order = 1,nodes = goTo[j],mode = c("in"),mindist = 0)[1])[-1]
#   return(list(Gene=goTo[j],Neighbors=Neighbors))
# })

## Iterate the Graph and Find the Neighbors with above average connectivity.
check=lapply(X = 1:length(goTo),FUN = function(j) {
  Neighbors=unlist(ego(graph = graph,order = 1,nodes = goTo[j],mode = c("in"),mindist = 0)[1])[-1]
  Neighbors=Neighbors[which(Neighbors>checkAVG)]
  return(list(Gene=goTo[j],Neighbors=Neighbors))
  })

##Remove elements with ZERO neighbors.
check=compact(.x = check,.p = function(x) x$Neighbors)

## To Save Time, Go Look for the Intersections. Specific to the Perturb Gene.
FindIt=lapply(X = check,function(p){
  XN_k=intersect(p$Neighbors,as.numeric(unlist(as.list(keydrivers))))
  if (length(XN_k)==0) {
    #print("Zero Intersection, Moving On")
  }
  else if (length(XN_k) == 1 & 12864 %in% XN_k) {
    pushback(x = keydrivers,data = p$Gene)
    #which(XN_k==12864)
    return(list(Node=p$Gene,Predictors=XN_k))
  }
  else {
    pushback(x = keydrivers,data = p$Gene)
    if (12864 %in% XN_k) {
      return(list(Node=p$Gene,Predictors=XN_k))  
    }
    else {
      return("Nothing")
    }
  }
})

## Remove the Emptiness.
FindItClean=compact(discard(FindIt, .p = is.character))

## Now we can iterate through and fit the models.
## But we are going to remove out any double single predictors. This is a waste of time.
FindItClean=keep(FindItClean,function(x) length(x$Predictors)>1)

## Remove unnecessary memory.
rm(list = c("FindIt","check","M_network","MEGENA_network","network","checkAVG"))

## Collect the layers.
layerList=lapply(1:1,function(h) unlist(ego(graph = graph,order = h,nodes = 12864,mode = c("all"),mindist = 0)[1])[-1])
counter=2
while (has_element(.x = layerList,.y = integer(0))!=TRUE) {
  layerList=lapply(1:counter,function(h) setdiff(unlist(ego(graph = graph,order = h,nodes = 12864,mode = c("all"),mindist = 0)[1])[-1],unlist(ego(graph = graph,order = h-1,nodes = 12864,mode = c("all"),mindist = 0)[1])[-1]))
  counter=counter+1
}

## Remove the trailing zero that halted the While Loop.
layerList=purrr::compact(.x = layerList)

## Find the paths from drivers to the Layers.
## From the documentation, beware, "Note that potentially there are exponentially many paths between two vertices of a graph, and you may run out of memory when using this function, if your graph is lattice-like."
## we have to prevent  ^^ that from happening, if it does.
findPaths=lapply(FindItClean, function(a){
  findPaths2=lapply(layerList, function(b){
    paths=igraph::all_simple_paths(graph = graph,from = a$Predictors,to = b)
    if (length(paths)==0) {
      
    }
    else{
      return(as.numeric(paths[[1]]))
    }
  })
  return(findPaths2)
})

## Remove elements without any paths. These are disconnected nodes.
findPaths=rlist::list.clean(findPaths,recursive=TRUE)

## We have to remove the paths WITHOUT perturb Gene / ATP6V1A and only 2 member pathways.
## This whole thing is built on you having a perturb gene to begin with.
for (p in 1:length(findPaths)) {
  findPaths[[p]]=keep(findPaths[[p]],function(x) 12864 %in% x)
  findPaths[[p]]=keep(findPaths[[p]],function(x) length(x)>2)
}

## Rolling Regression from drivers via viable pathways to distant gene layers.
withWindow=slider::slide(.x = findPaths[[1]][[1]],.after = 2,.before = Inf,.complete = TRUE,.f = function(c) {
  
  findSize=length(c)
  findGene=as.numeric(c[findSize])
  predictors=as.numeric(c[1:(findSize-1)])
  
  depend=features[,findGene]
  indep=features[,predictors]
  
  ## For indexing
  theSeq=seq(1:length(labels))
  
  ## For Case
  theSamples=which(labels==1)
  toChoice=sample(x = theSamples,size = floor(length(theSamples)*.90))
  ## For Case -- Hold Out
  toPredict=theSamples[-toChoice]
  
  ## For Control
  theSamples2=theSeq[-theSamples]
  toChoice2=sample(x = theSamples2, size = floor(length(theSamples2)*.90))
  ## For Control -- Hold Out
  toPredict2=theSamples2[-toChoice2]
  
  ## Now merge it back together and shuffle it.
  toChoiceFinal=features[c(toChoice,toChoice2),]
  toChoiceFinal=toChoiceFinal[sample(nrow(toChoiceFinal)),]
  
  ## Partition. 
  choice=caret::createDataPartition(y = 1:dim(toChoiceFinal)[1],times = 1,p = 0.75)[[1]]
  
  ## Partitioning into Train vs Test.
  training=indep[-choice,]
  testing=indep[choice,]
  
  ## FIT A RIDGE/EN/LASSO REGRESSION WITH CROSS VALIDATION (k folds).
  modelFit=glmnet::cv.glmnet(x = as.matrix(training),y = as.numeric(depend[-choice]),alpha=0.0,family="gaussian",nfolds = 5)
  
  ## Make predictions on test.
  yHat=predict(modelFit, newx=as.matrix(testing),s=modelFit$lambda.min,type="response")
  
  #Collect the information.
  res=list(Gene=findGene,model=modelFit,CasePredict=toPredict,ControlPredict=toPredict2,predictors=predictors)
  #res=list(Gene=findGene,MSE=mean(yHat-depend[choice])^2)
  
  return(res)
})

## Drop if any Null.
withWindow=compact(withWindow)
