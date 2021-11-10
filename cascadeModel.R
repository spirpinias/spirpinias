 ## Stephen And Qi Predictive/Cascade Network
## Adaptive Lasso/Elastic Net with graphical (MEGENA) prior.
## Debatable eQTL usage.

require(igraph)
require(purrr)
require(magrittr)
require(dplyr)
require(tibble)
require(ggplot2)

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
exp=exp[match(intersect(genesInvolved,row.names(exp)),row.names(exp)),]

## Encode the network based on indices.
MEGENA_network=sapply(MEGENA_network[,-3],function(x) match(x,row.names(exp)))

## Remove NA from MEGENA 
MEGENA_network=MEGENA_network[-which(is.na(MEGENA_network),arr.ind = T)[,1],]

## Ready
network=MEGENA_network

#Importing Expression Metadata to extract AD/Normal Labels.
ATP6V1Ameta=read.table(file = '/home/spirpinias/Desktop/BinZhang/ANN/Rcode/BM36/msbb.meta.BM_36.tsv', sep = '\t', header = TRUE)
ATP6V1Ameta=ATP6V1Ameta[-1,]

## Recoding Diagnosis as Numerical
ATP6V1Ameta$Dx.by.braak.cerad[which(ATP6V1Ameta$Dx.by.braak.cerad=="Normal")]=0
ATP6V1Ameta$Dx.by.braak.cerad[which(ATP6V1Ameta$Dx.by.braak.cerad=="ADpp")]=1
labels=as.numeric(ATP6V1Ameta$Dx.by.braak.cerad)

## Removing Labels which Diagnosis is NA.
exp=exp[,-which(is.na(labels)==TRUE)]
labels=labels[-which(is.na(labels)==TRUE)]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Read in the Directed Graph.
graph=graph_from_edgelist(as.matrix(network[,1:2],ncol=2), directed =T)

## Adjacency Matrix.
ad.matrix=as_adjacency_matrix(graph, type = c("both"))

## Degree Matrix.
de=degree(graph)

## Collect the layers outside ATP6V1A or 12864. Element of List is a Layer.
layerList=lapply(1:1,function(h) c(unlist(ego(graph = graph,order = h,nodes = 12864,mode = c("all"),mindist = 0)[1])[-1],12864))
counter=2
while (has_element(.x = layerList,.y = integer(0))!=TRUE) {
  layerList=lapply(1:counter,function(h) setdiff(unlist(ego(graph = graph,order = h,nodes = 12864,mode = c("all"),mindist = 0)[1])[-1],unlist(ego(graph = graph,order = h-1,nodes = 12864,mode = c("all"),mindist = 0)[1])[-1]))
  counter=counter+1
}

## Remove the trailing zero that halted the While Loop.
layerList=purrr::compact(.x = layerList)

## Perturb the expression of your gene of interest/Reduce by half.
exp[12864,]=exp[12864,]*1/2

##Get dimension.
samples=1:dim(exp)[2]

## Create a Distance Metric 
euclidean=function(a, b) sqrt(sum((a - b)^2))

for (z in 1:1) {
  
  ## October 27 11:10:32
  ## This model predicts Gene K in Layer A using -K Genes.
  
  ##Recopy the expression matrix to start again.
  updateMe=exp
  
  ## Whatever this is but will soon be gone. -Qi
  network=pred.res=pred.res0=NULL
  
  ## List of Lists - Layer x Genes.
  tryMe=purrr::map(.x = layerList,.f = function(a) {
    
    tryMeInner=purrr::map(.x = a, .f = function(x) {
      
    ## Random Sampling.
    intrain<-sample(1:length(samples),size=0.8*length(samples))
    
    ## Training Test Split.
    training<-updateMe[,intrain]
    testing<-updateMe[,-intrain]
    
    ## Gene to predict
    id=x
    
    ## Retrieving indices
    ydata0=training[id,]
    
    xdata0=training[-id,]
    
    ydata1=testing[id,]
    
    xdata1=testing[-id,]
    
    ## Feature Selection
    ## Target : id
    
    ## ACCORDING TO THE GRAPH, SELECT 3 NODES IN NEIGHBORHOOD OF VERTEX ID=i
    neighbor.id3=ego(graph, 3, nodes = id, mode = c("in"),
                     mindist = 0)[[1]][-1]
    
    ## ACCORDING TO THE GRAPH, SELECT 2 NODES IN NEIGHBORHOOD OF VERTEX ID=i
    neighbor.id2=ego(graph, 2, nodes = id, mode = c("in"),
                     mindist = 0)[[1]][-1]
    
    ## ACCORDING TO THE GRAPH, SELECT 1 NODES IN NEIGHBORHOOD OF VERTEX ID=i
    neighbor.id1=ego(graph, 1, nodes = id, mode = c("in"),
                     mindist = 0)[[1]][-1]
    
    
    ## GET VERTICES IN NEIGHBORHOODS THAT DO NOT SHARE SIMILARITIES. 3 and 2.
    neighbor.id3=setdiff(neighbor.id3,neighbor.id2)
    
    ## GET VERTICES IN NEIGHBORHOODS THAT DO NOT SHARE SIMILARITIES. 2 and 1.
    neighbor.id2=setdiff(neighbor.id2,neighbor.id1)
    
    ## FINAL VERTICES CONVERT THESE NUMBERS TO NUMERIC -- PROBABLY FOR LATER INDEXING
    neighbor.id1=as.numeric(neighbor.id1)
    
    ## IF THE DEGREE AT VERTEX i is < 2. There is only 1 EDGE.
    if(de[id]<2){
      
      selected.id=neighbor.id1[de[neighbor.id1]>de[id]]
      
      if(length(selected.id)>0){
        
        network=rbind(network,cbind(selected.id,id))
        
      }
      #if the degree is >2 do this.
    }else{
      ## THIS ALGORITHM IS DEFAULTING DEPENDING ON WEIGHT NUMERICAL VALUES//NON SPECIFIC
      ## IF NO WEIGHTS BFS used, otherwise if ALL POSITIVE, Dijkstra's.
      ## Mode is in. Has others.
      p=distances(graph, v = id, to = V(graph), mode = "in")
      
      ## Use 10 if path is infinity.
      p[p==Inf]=10
      
      ## Has something to do with weight..
      W=log(p+1)+2
      
      ## Re-Label the Degree Matrix.
      d=de
      
      ## It's doing something with the degrees.
      x3=d[neighbor.id3][d[neighbor.id3]>d[id]]
      
      x2=d[neighbor.id2][d[neighbor.id2]>d[id]]
      
      ## If the neighbor has higher degree than the target.
      x1=d[neighbor.id1][d[neighbor.id1]>d[id]]
      
      if(length(c(1/x1,2/x2,4/x3))>0){
        ## If there are neighbors who have higher degrees than the targets.
        W[d>d[id]]=max(1/x1,2/x2,4/x3)*5+p[d>d[id]]/d[d>d[id]]
        
        W[neighbor.id3[d[neighbor.id3]>d[id]]]=4/x3
        
        W[neighbor.id2[d[neighbor.id2]>d[id]]]=2/x2
        
        W[neighbor.id1[d[neighbor.id1]>d[id]]]=1/x1
        
      }else{
        ## If the neighbors have less degree than the target.
        W[d>d[id]]=p[d>d[id]]/d[d>d[id]]
        
      }
      
      W=W[-id]

      ## Fit the First Model.
      modelFit=glmnet::glmnet(x = t(xdata0),y = as.numeric(ydata0),family = "gaussian",alpha = 0.5,penalty.factor = W,weights = labels[intrain],type.measure = "mse")
      
      ## Extract Coefficients
      features=coef(modelFit,s=0.05,exact=FALSE)  
      
      ## Extract Nonzero Coefficients.
      features=features[features[,1]!=0,]
      
      ## Get Coefficients without Bias.
      geneSelected=features[-1]
      
      ## If the predictors are greater than 10, refit. We need to slim it down.
      if (length(geneSelected)>10) {
        
        ## K Top Predictors and refit the model.
        posFeature=names(head(sort(geneSelected,decreasing = TRUE)))
        negFeature=names(tail(sort(geneSelected,decreasing = TRUE)))
        newPredictors=c(posFeature,negFeature)
        newPredictors=match(newPredictors,row.names(exp))
        newPredictors=unique(posFeature,negFeature)
        
        ## Subset data for K Top Predictors.
        ydata0=training[id,]
        
        xdata0=training[newPredictors,]
        
        ydata1=testing[id,]
        
        xdata1=testing[newPredictors,]
        
        ## Re-Fit the model.
        modelFit2=glmnet::glmnet(x = t(xdata0),y = as.numeric(ydata0),family = "gaussian",alpha = 0.5,weights = labels[intrain],type.measure = "mse")
        
        ## Arbitrary choice of lambda. It should be studied at min before anything.
        yHat2=predict(modelFit2,t(xdata1),s=0.05,type="response")
        
        ## Extract Coefficients
        features2=coef(modelFit2,s=0.05,exact=FALSE)
        
        ## Extract Nonzero Coefficients.
        features2=features2[features2[,1]!=0,]
        
        ## Get Coefficients without Bias.
        geneSelected2=features2
        
        ## Update the Expression (UpdateMe) Matrix by parameterizing a normal distribution.
        top2=rnorm(n = dim(updateMe)[2],mean = mean(yHat2),sd = sd(yHat2))
        updateMe[id,]<<-top2
        
        ## Collect the information and keep track of progress.
        print(paste("Iteration",z,"Gene",x,"in Layer",length(a)))
        
        return(list(yHat=as.numeric(yHat2),actual=as.numeric(ydata1),gene=id,coefficients=geneSelected2))
      }
      else{
        
        ## If there are less than 10 predictors, no refit necessary.
        ## Proceed with predictions, updates, and saving information.
        
        ## Predict the test using lambda.min
        yHat=predict(modelFit,t(xdata1),s=0.05,type="response")
        
        ## Update the Expression Matrix by parameterizing a normal distribution.
        top2=rnorm(n = dim(updateMe)[2],mean = mean(yHat),sd = sd(yHat))
        updateMe[id,]<<-top2
        
        ## Collect the information and keep track of progress.
        print(paste("Iteration",z,"Gene",x,"in Layer",length(a)))
        
        return(list(yHat=as.numeric(yHat),actual=as.numeric(ydata1),gene=id,coefficients=geneSelected))
      }
    }
    })
    return(tryMeInner)
  })
}
