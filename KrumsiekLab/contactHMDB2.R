## This is a program that I wrote in the Laboratory of Jan Krumsiek for using contactHMDB.py ID's.
## The File returned by contactHMDB.py is fed into this script for further processing.

## InputFile = Output of HMDB.py ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

generateIDs=function(inputFile){
#inputFile='/home/spirpinias/myTest.csv'
  
#Importing the ID's for Analysis
library(hmdbQuery)
metaboliteIDs=read.csv(inputFile)
metaboliteIDs=metaboliteIDs[,-1]
metaboliteIDs=metaboliteIDs[!(is.na(metaboliteIDs$Match) | metaboliteIDs$Match==""),]
return(metaboliteIDs)
}

#Contact the Database.~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

metaboliteIDs=generateIDs(inputFile = '/home/spirpinias/myTest.csv')
generateDiseases=function(j){
  testingVariable=hmdbQuery::HmdbEntry(prefix = "http://www.hmdb.ca/metabolites/",id=metaboliteIDs[j,3],keepFull = TRUE)
  hmdbQuery::diseases(testingVariable)
}

#Error Exception Handling ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
library(purrr)

possiblySomeFunction=possibly(generateDiseases,otherwise = "No Information")
possiblySomeInfo=map(1:8,possiblySomeFunction)

#Convert to Disease by Citation for Plotting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
testFrame=do.call(rbind.data.frame,map2(.x = possiblySomeInfo[[5]][1,],.y = possiblySomeInfo[[5]][3,],~cbind(print(.x),length(.y))))
colnames(testFrame)=c("Disease","Cited")
testFrame$Cited=as.numeric(testFrame$Cited)

#Generate the Circular Plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Entirely off example given on ggplot2 circular part tutorial. I wanted the identical plot on the page. It fit the graphic perfectly.
data = data.frame(
  id=seq(1,6),
  individual=testingFrame$Disease,
  value=testingFrame$Cited
)
label_data=data
number_of_bar=nrow(label_data)
angle=90-360*(label_data$id-0.5)/number_of_bar     
label_data$hjust=ifelse( angle < -90, 1, 0)
label_data$angle=ifelse(angle < -90, angle+180, angle)
  
# Start the plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p <- ggplot(data, aes(x=as.factor(id), y=value)) + geom_bar(stat="identity", fill=alpha("skyblue", 0.7))+ylim(-40,120) + theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm")      # Adjust the margin to make in sort labels are not truncated!
  ) + coord_polar(start = 0) +
  
  geom_text(data=label_data, aes(x=id, y=value+10, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) 
 
  # Print the Plot
  p
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
