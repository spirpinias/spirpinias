## This is currently a project of mine. The idea is simple.
## A linear model makes a prediction of a target using variables you think are important. 
## What if I used that prediction to parameterize a normal distribution?
## From that distribution, draw the number of samples needed to update the target vector, and continue on modeling.
## I am introducing a ripple into the pond. 
## I have stuck with elastic net because of the unique solution. 
## "Moreover, for any α < 1 and λ > 0, the elastic-net problem (4.2) is strictly convex: a unique
  ## solution exists irrespective of the correlations or duplications in the Xj ." - Statistical Learning with Sparsity Hastie, Tibshirani, and Wainwright.



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                      Predictive Networks with Weights 
##                       Developed by Stephen Pirpinias
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require(igraph)
require(purrr)
require(dplyr)
require(tibble)
require(ggplot2)
require(tidyr)
require(limma)
require(data.table)
require(magrittr)
require(gt)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                          Functions 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

findGene=function(geneName){
  
  ## Example testing=findGene("ATP6V1A") 
  firstX=which(row.names(lookUpTable)==as.character(geneName))
  secondX=lookUpTable[firstX,]
  thirdX=collectUpdates[[1]][[as.numeric(secondX$Inner)]][[as.numeric(secondX$Outer)]]
  thirdX$myCoefs=names(thirdX$myCoefs)[-1]
  return(thirdX)
}

buildSubNetwork=function(target,predictors){
  require(visNetwork)
  
  ## Example testing2=buildSubNetwork(testing$Gene,testing$myCoefs)
  center=which(row.names(exp)==as.character(target))
  neigh=match(x = predictors,table = row.names(exp))
  subGraph=induced.subgraph(graph = graph,v = as.numeric(c(center,neigh)),impl = "auto")
  V(subGraph)[1]$color="red"
  outGraph=plot(x = subGraph,vertex.size=28,vertex.shape="circle",vertex.label=row.names(exp)[c(center,neigh)],vertex.label.cex=1)
  
  return(outGraph)
}

buildScatterPlot=function(Gene){
  
  ## UpdatedMe/Original
  scat=which(row.names(exp)==as.character(Gene))
  Skim=data.frame(Expression=c(as.numeric(SimControl[[1]][[1]][scat,1:27]),as.numeric(SimControl[[1]][[1]][scat,28:54])),
                  Labels=SimControl[[1]][[2]][,1])
  Skim$Labels=dplyr::recode(Skim$Labels, "1"="Simulated", "0"="Original")
  p=ggplot(data = Skim,mapping = aes(x=Expression,color=Labels))+geom_density()+labs(title=paste("Distribution of",as.character(Gene)))
  return(p)  
}


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                          Data Preparation 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Set directory for files.
dir="/home/spirpinias/Desktop/BinZhang/pn/"
setwd(dir)

## Import Expression Dataset.
exp=fread(file = "msbb.BM_36.CDR_adjusted.PMI_AOD_race_sex_RIN_exonicRate_rRnaRate_batch_adj.tsv") %>% as.data.frame()
exp=column_to_rownames(.data = exp,var = "V1")

## MEGENA network 
MEGENA_network=fread("BM36_CDR_adj_himem.onlyclust.megena_out.tsv",header = T,drop=1) %>% as.data.frame
colnames(MEGENA_network)=c("Start","Finish","Weight")

## Use only genes in MEGENA and subset expression. Some are NA, remove them.
genesInvolved=union(x = MEGENA_network[,1],y = MEGENA_network[,2])
exp=exp[match(genesInvolved,row.names(exp)),]
exp=exp[grepl("^NA", rownames(exp))==F,]

## MEGENA keydrivers
MEGENA_keydrivers=read.delim("all_BN_keydrivers_CTD_subtypes.txt",header = T)
MEGENA_keydrivers=data.frame(Keydriver=intersect(MEGENA_keydrivers$keydrivers,row.names(exp)),Index=match(intersect(MEGENA_keydrivers$keydrivers,row.names(exp)),row.names(exp)))


## Encode the network based on indices.
MEGENA_network=data.frame(Start=match(MEGENA_network$Start,row.names(exp)),Finish=match(MEGENA_network$Finish,row.names(exp)),Weight=MEGENA_network$Weight) %>% na.omit() %>% as.data.frame()


## Ready
network=MEGENA_network

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                             Create the Graph
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Read in the Directed Graph.
graph=graph_from_edgelist(as.matrix(network[,1:2],ncol=2), directed =T)
#graph=graph_from_data_frame(d = network,directed = FALSE)

## Adjacency Matrix. Maybe used this to weight the variables.
#ad.matrix=as_adjacency_matrix(graph, type = c("both"))

## Weight Adjanceny Matrix.
##as_adjacency_matrix(graph = graph,attr = "weight")

## Degree Matrix. Maybe use this to weight the variables.
de=degree(graph)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                           Meta Data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Importing Expression Metadata to extract AD/Normal Labels.
ATP6V1Ameta=read.table(file = '/home/spirpinias/Desktop/BinZhang/ANN/Rcode/BM36/msbb.meta.BM_36.tsv', sep = '\t', header = TRUE)
ATP6V1Ameta=ATP6V1Ameta[,-1]

## Recoding Diagnosis as Factor.
ATP6V1Ameta$Dx.by.braak.cerad=dplyr::recode(ATP6V1Ameta$Dx.by.braak.cerad, "ADpp"=1, "Normal"=0)
labels=as.numeric(ATP6V1Ameta$Dx.by.braak.cerad)

## Removing Labels which Diagnosis is NA for the samples.
exp=exp[-which(de==0),-which(is.na(labels)==TRUE)]
ATP6V1Ameta=ATP6V1Ameta[-which(is.na(labels)==TRUE),]
labels=labels[-which(is.na(labels)==TRUE)]


## Remove the genes with zero degree from connectivity and graph.
graph=delete_vertices(graph = graph,v = which(de==0))
de=de[-which(de==0)]

## Index to build the Layer List
aroundMe=which(row.names(exp)=="ATP6V1A")

## Subet for Control Only. Testing phase.
exp0=exp[,which(labels==0)]
exp1=exp[,which(labels==1)]

# ## Simulation.
splitHere=which(as.numeric(exp0[aroundMe,]<mean(as.numeric(exp0[aroundMe,])))==1)
exp0Low=exp0[,splitHere]
exp0High=exp0[,-splitHere]

splitHere2=which(as.numeric(exp1[aroundMe,]<mean(as.numeric(exp1[aroundMe,])))==1)
exp1Low=exp1[,splitHere2]
exp1High=exp1[,-splitHere2]


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
collectTablesBoot=list()
lookUpTable=data.frame()

## Vectors for KD or OE.
vectorKnockDown=c(1/2)
OverExpress=c(1.50)
numNeighbor=1

## Only Keydrivers in 3 layers. TEMPORARY.
MEGENA_keydrivers=MEGENA_keydrivers[match(intersect(unlist(layerList[1:3]),MEGENA_keydrivers$Index),MEGENA_keydrivers$Index),]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  Modeling
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Boot Strapping for 10 iterations.
for (b in 1:2) {
  
  for (e in 1:length(vectorKnockDown)) {
    
    ## Make a copy.
    updateMe=exp0High

    testingModel=imap(head(layerList,2),function (d,d2){
      
        
      imap(d, function(j,j2){
        
        if (b==1) {
          lookUpTable<<-rbind(lookUpTable,c(row.names(exp)[j],c(d2,j2)))
        }
        
        ## Recopy TYROBP to correct after perturbing.
        updateMe[aroundMe,]=exp0High[aroundMe,]
        
        ## Split the Data.
        splitData=dim(updateMe)[2]
        
        ## Partition
        samples=sample(1:splitData,size=splitData*0.80)
        
        ## Train/Test
        train=updateMe[,samples]
        test=updateMe[,-samples]
        
        ## Intersect with keydrivers to find the variables.
        indep=intersect(unlist(ego(graph = graph,order = numNeighbor,nodes = j,mode = c("all")))[-1],MEGENA_keydrivers$Index)
        
        if (length(indep)>=2) {
          
          # Modeling
          modelFit=glmnet::glmnet(x = as.matrix(t(train[indep,])),y = as.numeric(train[j,]),alpha = 0.5,nlambda = 100)
          
          ## Coefficients
          myCoefs=coef(modelFit,s=0.05,exact=FALSE)  
          myCoefs=myCoefs[myCoefs[,1]!=0,]
          
          ## I have to perturb the gene right here and capture the effect.
          test[aroundMe,]=test[aroundMe,]*vectorKnockDown[e]
          
          ## Make the prediction.
          yHat=predict(modelFit,as.matrix(t(test[indep,])),s=0.05,type="response")
          MSE=mean((as.numeric(test[j,])-as.numeric(yHat))^2)
          
          ## Updating scheme.
          top2=rnorm(n = dim(updateMe)[2],mean = mean(yHat),sd = sd(yHat))
          
          ## Save it.
          updateMe[j,]<<-top2
          
          print(paste("Gene",j,"Iteration",b))
          
          res=list(Gene=row.names(exp)[j],MSE=MSE,myCoefs=myCoefs,Driver=TRUE)
        }
        else{
          
          ## Key Drivers were null. Use connectivity. 
          connectMe=data.frame(Index=unlist(ego(graph = graph,order = numNeighbor,nodes = j,mode = c("all")))[-1],Connect=de[unlist(ego(graph = graph,order = numNeighbor,nodes = j,mode = c("all")))[-1]])
          
          ## Sort the dataframe in order of highest connectivity.
          connectMe=connectMe[order(-connectMe$Connect),]
          indepConnect=head(connectMe$Index,n=5)
          
          if(length(indepConnect)>= 2) {
            ## Modeling
            modelFit=glmnet::glmnet(x = as.matrix(t(train[indepConnect,])),y = as.numeric(train[j,]),alpha = 0.5,nlambda = 100)
            
            ## Coefficients
            myCoefs=coef(modelFit,s=0.05,exact=FALSE)  
            myCoefs=myCoefs[myCoefs[,1]!=0,]
            
            ## I have to perturb the gene right here and capture the effect.
            test[aroundMe,]=test[aroundMe,]*vectorKnockDown[e]
            
            ## Make the prediction.
            yHat=predict(modelFit,as.matrix(t(test[indepConnect,])),s=0.05,type="response")
            MSE=mean((as.numeric(test[j,])-as.numeric(yHat))^2)
            
            ## Updating scheme.
            top2=rnorm(n = dim(updateMe)[2],mean = mean(yHat),sd = sd(yHat))
            
            ## Save it.
            updateMe[j,]<<-top2
            
            print(paste("Gene",j,"Iteration",b))
            
            res=list(Gene=row.names(exp)[j],MSE=MSE,myCoefs=myCoefs,Driver=FALSE)
            
          }
          else{
            print('Not enough Genes in Neighborhood')
          }
        }
        
      })
    })
  
    ## We have to re-attach the case to the modeled.
    design=model.matrix(~0+c(rep(1,dim(updateMe)[2]),rep(0,dim(exp0High)[2]))+c(rep(0,dim(updateMe)[2]),rep(1,dim(exp0High)[2])))
    colnames(design)=c("Sim","Org")
    SimControl[[e]]=list(cbind(updateMe,exp0High),design)
    
    ## Do differential expression.
    contr.matrix=makeContrasts(SimvsOrg=Sim-Org,levels=colnames(design))
    vfit=lmFit(SimControl[[e]][[1]],design)
    vfit=contrasts.fit(vfit,contrasts=contr.matrix)
    tfit=treat(vfit,lfc=1)
    theTable=topTreat(tfit,coef=1,n=Inf)
    theTable=tibble::rownames_to_column(.data = theTable,var = "Geneid")
    theTable=data.table(theTable)
    #theTable=theTable[adj.P.Val<0.05] 
    
    ## Load the containers. 
    collectTables[[e]]=theTable
    collectUpdates[[e]]=testingModel
  }
  
  collectTablesBoot[[b]]=collectTables 
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                              Cell Line Validation - TYROBP
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## To collect the bootstraps and fix the look up Table.
colnames(lookUpTable)=c("Gene","Inner","Outer")
lookUpTable=tibble::column_to_rownames(.data = lookUpTable,var = "Gene")

collectTablesBoot=pmap(collectTablesBoot,rbind)
collectTablesBoot=map(collectTablesBoot, function(b) unique(b,by="Geneid"))


## TYROBP Cell Line ... do unique here.
TYROBP=fread("/home/spirpinias/Downloads/TYROBP.KO.tsv",header = TRUE,sep='\t')

## Subset TYROBP Contrasts into respective dataframes and clean out any repeats.
TYROKOWT=TYROBP[Contrast=="KO.vs.WT"]
TYROKOWT=unique(TYROKOWT,by="Geneid")
TYROAPPWT=TYROBP[Contrast=="AppPs1KO.vs.WT"]
TYROAPPWT=unique(TYROAPPWT,by="Geneid")
TYROPs1KO1WT=TYROBP[Contrast=="AppPs1KO.vs.AppPs1WT"]
TYROPs1KO1WT=unique(TYROPs1KO1WT,by="Geneid")

## Uppercase the gene names.
TYROKOWT$Geneid=toupper(TYROKOWT$Geneid)
TYROAPPWT$Geneid=toupper(TYROAPPWT$Geneid)
TYROPs1KO1WT$Geneid=toupper(TYROPs1KO1WT$Geneid)

## As to ensure the tables are in identical order.
first=TYROKOWT[TYROAPPWT,on=.(Geneid),.(Geneid,logFC,i.logFC)]
second=first[TYROPs1KO1WT,on=.(Geneid),]
TYROBP=second[,.(Geneid,KOvsWT=logFC,APPvsWT=i.logFC,s1K01vsWT=i.logFC.1)]

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                              Cell Line Validation - ATP6V1A 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## To collect the bootstraps and fix the lookup table.
colnames(lookUpTable)=c("Gene","Inner","Outer")
lookUpTable=tibble::column_to_rownames(.data = lookUpTable,var = "Gene")

collectTablesBoot=pmap(collectTablesBoot,rbind)
collectTablesBoot=unique(collectTablesBoot,by="Geneid")

## Temporary -- March 29 -- For Subset Layers
collectTablesBoot[[1]]=collectTablesBoot[[1]][match(row.names(exp)[unlist(layerList[1:2])],collectTables[[1]]$Geneid),]


## ATP Neuronal Line
ATP6V1A=fread("/home/spirpinias/Downloads/ATP6V1A.KD.tsv",header=TRUE,sep='\t')
colnames(ATP6V1A)[2]="Geneid"

## Subset ATP6V1A 
ATPkdVSwt_in_ab=ATP6V1A[Contrast=="KD.vs.WT.in.Ab"][,-1]
ATPkdVSwt_in_ab=unique(ATPkdVSwt_in_ab,by="Geneid")
ATPab_kdVSWT_unchallen=ATP6V1A[Contrast=="Ab_KD.vs.WT_unchallenged"][,-1]
ATPab_kdVSWT_unchallen=unique(ATPab_kdVSWT_unchallen,by="Geneid")
ATPkd_vs_wd_in_unchallenged=ATP6V1A[Contrast=="KD.vs.WT.in.unchallenged"][,-1]
ATPkd_vs_wd_in_unchallenged=unique(ATPkd_vs_wd_in_unchallenged,by="Geneid")

## Uppertcase the Gene Names
ATPkdVSwt_in_ab$Geneid=toupper(ATPkdVSwt_in_ab$Geneid)
ATPab_kdVSWT_unchallen$Geneid=toupper(ATPab_kdVSWT_unchallen$Geneid)
ATPkd_vs_wd_in_unchallenged$Geneid=toupper(ATPkd_vs_wd_in_unchallenged$Geneid)

## ATP6V1A
first=ATPkdVSwt_in_ab[ATPab_kdVSWT_unchallen,on=.(Geneid),.(Geneid,logFC,i.logFC)]
second=first[ATPkd_vs_wd_in_unchallenged,on=.(Geneid),]
ATP6V1A=second[,.(Geneid,ATPkdVsWT_inAB=logFC,ATPab_kdVsWT_unchallen=i.logFC,ATPkdVsWd_in_Unchallenged=i.logFC.1)]


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                              Results & Visualizations
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## collectTablesBoot

## Results. TYROBP.
ResultsFinal=map(.x = collectTablesBoot,.f = function(b) {
  finalFrame=inner_join(x = TYROBP,y = b,by="Geneid")
  finalFrame=finalFrame[,.(Geneid,KOvsWT,APPvsWT,s1K01vsWT,Simulated=logFC)]
  return(finalFrame)  
})

## Results. ATP6V1A
ResultsFinal=map(.x = collectTablesBoot,.f = function(b) {
  finalFrame=inner_join(x = ATP6V1A,y = b,by="Geneid")
  finalFrame=finalFrame[,.(Geneid,ATPkdVsWT_inAB,ATPab_kdVsWT_unchallen,ATPkdVsWd_in_Unchallenged,Simulated=logFC)]
  return(finalFrame)  
})


## Plot the Results.
colors=c('Simulated'='black','KOvsWT'='red','APPvsWT'='green','s1K01vsWT'='yellow')

ResultsPlots=map(.x = ResultsFinal,.f = function(c) {
  require(ggplot2)
  ggplot(c,aes(x=1:dim(c)[1]))+
    geom_line(aes(y=KOvsWT, color="KOvsWT"),size=0.4)+
    geom_line(aes(y=APPvsWT, color="APPvsWT"),size=0.4)+
    geom_line(aes(y=s1K01vsWT, color="s1K01vsWT"),size=0.4)+
    geom_line(aes(y=Simulated, color="Simulated"),size=0.4)+
    labs(x="Genes Predicted",y="Log Fold Change",title="Cascade Effect for TYROBP Knockout",subtitle="Elastic Net",color="Legend")+
    scale_color_manual(values=colors)
  #ggsave("/home/spirpinias/Desktop/ATPknockDown",device = "png",width=6,height=6)
})


colors=c('Simulated'='black','ATPkdVsWT_inAB'='red','ATPab_kdVsWT_unchallen'='green','ATPkdVsWd_in_Unchallenged'='yellow')

ResultsPlots=map(.x = ResultsFinal,.f = function(c) {
  require(ggplot2)
  ggplot(c,aes(x=1:dim(c)[1]))+
    geom_line(aes(y=ATPkdVsWT_inAB, color="ATPkdVsWT_inAB"),size=0.4)+
    geom_line(aes(y=ATPab_kdVsWT_unchallen, color="ATPab_kdVsWT_unchallen"),size=0.4)+
    geom_line(aes(y=ATPkdVsWd_in_Unchallenged, color="ATPkdVsWd_in_Unchallenged"),size=0.4)+
    geom_line(aes(y=Simulated, color="Simulated"),linetype="dashed",size=0.4)+
    labs(x="Genes Predicted",y="Log Fold Change",title="Cascade Effect for ATP6V1A Knockdown",subtitle="Elastic Net",color="Legend")+
    scale_color_manual(values=colors)
  #ggsave("/home/spirpinias/Desktop/ATPknockDown",device = "png",width=6,height=6)
})
