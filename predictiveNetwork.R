## This is currently a project of mine. The idea is simple.
## A linear model makes a prediction of a target using variables you think are important. 
## What if I used that prediction to parameterize a normal distribution?
## From that distribution, draw the number of samples needed to update the target vector, and continue on modeling.
## I am introducing a ripple into the pond. 
## I have stuck with elastic net because of the unique solution. 
## "Moreover, for any α < 1 and λ > 0, the elastic-net problem (4.2) is strictly convex: a unique
  ## solution exists irrespective of the correlations or duplications in the Xj ." - Statistical Learning with Sparsity Hastie, Tibshirani, and Wainwright.



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                      Predictive Networks 
##                  Developed by Stephen Pirpinias
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
require(igraph)
require(purrr)
require(dplyr)
require(tibble)
require(ggplot2)
require(tidyr)
require(limma)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                          Data Preparation 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
## Set directory for files.
dir="/home/spirpinias/Desktop/BinZhang/pn/"
setwd(dir)

## Import Expression Dataset.
exp=read.table("msbb.BM_36.CDR_adjusted.PMI_AOD_race_sex_RIN_exonicRate_rRnaRate_batch_adj.tsv",sep = '\t', header = TRUE)
exp=exp[,-1]

## MEGENA network 
MEGENA_network=read.delim("BM36_CDR_adj_himem.onlyclust.megena_out.tsv",header = T)

## Use only genes in MEGENA and subset expression.
genesInvolved=union(x = MEGENA_network[,1],y = MEGENA_network[,2])
exp=exp[match(genesInvolved,row.names(exp)),]
exp=exp[grepl("^NA", rownames(exp))==F,]

## Encode the network based on indices.
MEGENA_network=sapply(MEGENA_network[,-3],function(x) match(x,row.names(exp)))

## Remove NA from MEGENA 
MEGENA_network=tidyr::drop_na(data = as.data.frame(MEGENA_network))

## Ready
network=MEGENA_network

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                           Meta Data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Importing Expression Metadata to extract AD/Normal Labels.
ATP6V1Ameta=read.table(file = '/home/spirpinias/Desktop/BinZhang/ANN/Rcode/BM36/msbb.meta.BM_36.tsv', sep = '\t', header = TRUE)
ATP6V1Ameta=ATP6V1Ameta[-1,]

## Recoding Diagnosis as Factoral.
ATP6V1Ameta$Dx.by.braak.cerad=dplyr::recode(ATP6V1Ameta$Dx.by.braak.cerad, "ADpp"=1, "Normal"=0)
labels=as.numeric(ATP6V1Ameta$Dx.by.braak.cerad)

## Removing Labels which Diagnosis is NA for the samples.
exp=exp[,-which(is.na(labels)==TRUE)]
labels=labels[-which(is.na(labels)==TRUE)]

## Index to build the Layer List
aroundMe=which(row.names(exp)=="TYROBP")

## Subet for Control Only. Testing phase.
exp0=exp[,which(labels==0)]
exp1=exp[,which(labels==1)]

## Simulation.
splitHere=which(as.numeric(exp0[aroundMe,]<mean(as.numeric(exp0[aroundMe,])))==1)
exp0Low=exp0[,splitHere]
exp0High=exp0[,-splitHere]

splitHere2=which(as.numeric(exp1[aroundMe,]<mean(as.numeric(exp1[aroundMe,])))==1)
exp1Low=exp1[,splitHere2]
exp1High=exp1[,-splitHere2]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                             Create the Graph
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Read in the Directed Graph.
graph=graph_from_edgelist(as.matrix(network[,1:2],ncol=2), directed =T)

## Adjacency Matrix.
ad.matrix=as_adjacency_matrix(graph, type = c("both"))

## Degree Matrix.
de=degree(graph)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                               Building the Layers
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Collect the layers outside a gene. 
layerList=lapply(1:1,function(h) c(unlist(ego(graph = graph,order = h,nodes = aroundMe,mode = c("all"),mindist = 0)[1])[-1],aroundMe))
counter=2
while (has_element(.x = layerList,.y = integer(0))!=TRUE) {
  layerList=lapply(1:counter,function(h) setdiff(unlist(ego(graph = graph,order = h,nodes = aroundMe,mode = c("all"),mindist = 0)[1])[-1],unlist(ego(graph = graph,order = h-1,nodes = aroundMe,mode = c("all"),mindist = 0)[1])[-1]))
  counter=counter+1
}

## Remove the trailing zero that halted the While Loop.
layerList=purrr::compact(.x = layerList)

## Containers. 
collectTables=list()
collectUpdates=list()
SimControl=list()
euclidList=list()

## Vectors for KD or OE.
vectorKnockDown=c(1.00,1/2,1/4,1/8,1/16,1/32)
OverExpress=c(1.00,1.25,1.50,1.75,2)
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  Modeling
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#length(vectorKnockDown)
for (e in 1:1) {
  
  ## Make a copy.
  updateMe=exp0
  
  ## Perturb it.
  updateMe[aroundMe,]=updateMe[aroundMe,]*vectorKnockDown[e]
  
  testingModel=map(layerList,function (d){
    map(d, function(j){
      
      splitData=dim(updateMe)[2]
      
      samples=sample(1:splitData,size=splitData*0.80)
      
      train=updateMe[,samples]
      test=updateMe[,-samples]
      
      modelFit=glmnet::glmnet(x = as.matrix(t(train[-j,])),y = as.numeric(train[j,]),alpha = 0.5,weights=labels[samples],nlambda = 100)
      
      myCoefs=coef(modelFit,s=0.05,exact=FALSE)  
      myCoefs=myCoefs[myCoefs[,1]!=0,]
      
      yHat=predict(modelFit,as.matrix(t(test[-j,])),s=0.05,type="response")
      MSE=mean((as.numeric(test[j,])-as.numeric(yHat))^2)

      top2=rnorm(n = dim(updateMe)[2],mean = mean(yHat),sd = sd(yHat))

      updateMe[j,]<<-top2
      
      print(paste("Gene",j,"Iteration",e))
      res=list(MSE=MSE,myCoefs=myCoefs)
    })
  })
  
  ## Calculate the Distances
  euclidMap=sapply(1:dim(updateMe)[1],function(a) {
    dist(rbind(updateMe[a,], exp0[a,]))
  })
  names(euclidMap)=row.names(exp0)
  
  ## Remove the Zero Distances.
  euclidMap=euclidMap[-which(euclidMap==0)]
  
  ## We have to re-attach the case to the modeled.
  design=model.matrix(~0+c(rep(1,dim(updateMe)[2]),rep(0,dim(exp0)[2]))+c(rep(0,dim(updateMe)[2]),rep(1,dim(exp0)[2])))
  colnames(design)=c("Sim","Org")
  SimControl[[e]]=list(cbind(updateMe,exp0),design)
  
  ## Do differential expression.
  contr.matrix=makeContrasts(SimvsOrg=Sim-Org,levels=colnames(design))
  vfit=lmFit(SimControl[[e]][[1]],design)
  vfit=contrasts.fit(vfit,contrasts=contr.matrix)
  tfit=treat(vfit,lfc=1)
  theTable=topTreat(tfit,coef=1,n=Inf)
  
  collectTables[[e]]=theTable
  collectUpdates[[e]]=updateMe
  euclidList[[e]]=euclidMap
}
