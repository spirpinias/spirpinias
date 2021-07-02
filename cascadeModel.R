##This is an iterative model that I created in the Laboratory of Bin Zhang.
##We use the edgelist from MEGENA to read into iGraph. 
##We start with the first node of the graph, call neighbors intersect with keydrivers, and if intersection is not null; build a model.
##The keydrivers are genes of importance from the MEGENA clustering package.


##Output of calculate.PFN from MEGENA package.
#Keydrivers are a set of genes MEGENA finds informative.
MEGENA_keydrivers=read.delim("/home/spirpinias/Desktop/BinZhang/pn/all_BN_keydrivers_CTD_subtypes.txt",header = T)
MEGENA_keydrivers=MEGENA_keydrivers$keydrivers
MEGENA_keydrivers=as.numeric(na.omit(match(x = MEGENA_keydrivers,table = row.names(exp))))

##Get the Network.
MEGENA_network=read.delim("BM36_CDR_adj_himem.onlyclust.megena_out.tsv",header = T)
M_network=cbind(match(MEGENA_network$row,gene_names),match(MEGENA_network$col,gene_names),MEGENA_network$weight)

#Remove NA from MEGENA network and M_network
M_network=M_network[-which(is.na(M_network),arr.ind = T)[,1],]

##Assign Indices of the Edge Information to a Dataframe.
network=as.data.frame(M_network)

#Some genes were removed by MEGENA, we no longer have to visit them.
goTo=as.numeric(union(x=network[,1],y=network[,2]))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(igraph)
library(dequer)
library(purrr)
library(tictoc)
library(ggplot2)
library(randomForest)

##Convert Pairs of Genes with Weights to Directed Graph.
graph=graph_from_edgelist(as.matrix(network[,1:2],ncol=2), directed=T)

##Remove the Keydrivers from the set of Genes to Predict.
goTo=setdiff(x = goTo,y = MEGENA_keydrivers)

##Keydrivers in queue for the ease of appending whilst maintaing memory efficiency and speed.
##The append function takes copies of the entire array and I did not want that behavior. The Queue was much better option.
keydrivers=as.queue(as.list(MEGENA_keydrivers))

#Neighborhood search is done only once to save time.
check=lapply(X = 1:length(goTo),FUN = function(j) list(Gene=goTo[j],Neighbors=unlist(ego(graph = graph,order = 1,nodes = goTo[j],mode = c("in"),mindist = 0)[1])[-1]))

##Remove elements with ZERO neighbors.
check=compact(.x = check,.p = function(x) x$Neighbors)

##cascadeModel
fitIt=lapply(X = check,function(p){
  XN_k=intersect(p$Neighbors,as.numeric(unlist(as.list(keydrivers))))
  if (length(XN_k)==0) {
    pushback(x = keydrivers,data = p$Gene)
    print("Zero Intersection, Appending, Moving on.")
  }
  else if (length(XN_k) == 1) {
    pushback(x = keydrivers,data = p$Gene)
    print("One Intersection, Appending, No Model")
  }
  else{
    print("Model Detected")
    first=p$Gene
    second=XN_k
    
    depend=t(exp)[,first] 
    indep=t(exp)[,second]
    
    #Partition.  
    choice=caret::createDataPartition(y = 1:214,times = 1,p = 0.70)[[1]]
    
    ##Partition it, again.
    training=indep[-choice,]
    testing=indep[choice,]
    
    #Fit a Random Forest with arbitrary amount of trees.
    modelFit=randomForest(x = as.matrix(training),y = as.numeric(depend[-choice]),ntree=1000,importance=TRUE)
    
    #Make Predictions
    yHat=predict(modelFit,as.matrix(testing),type="response")
    
    ##Append Gene to the Queue
    pushback(x = keydrivers,data = first)
    
    #Importance Matrix
    myCoefs=importance(modelFit,scale=TRUE)
    
    
    #Information
    res=list(Gene=first,MSE=mean(yHat-depend[choice])^2,myCoefs=myCoefs)
    return(res)
  }
})

## Clean it up. If you see characters, they are not with model. Delete them.
first=which(lapply(fitIt,function(x) is.character(x))==TRUE)

## Create a new list of only the models.
cleanList=fitIt[-first]

## Extract MSE and make a pretty plot.
cleanMSE=unlist(lapply(cleanList,function(x) x$MSE))

##Plotting
x=seq(1,length(cleanList))
y=cleanMSE
data=data.frame(x=x,y=y)
ggplot(data, aes(x=x, y=y)) +
  geom_point() + 
  geom_segment( aes(x=x, xend=x, y=0, yend=y))+xlab("Genes")+ylab("MSE")+ggtitle("Cascade Model")
