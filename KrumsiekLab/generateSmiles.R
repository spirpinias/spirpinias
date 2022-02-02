## This is a program that I wrote for converting metabolites into SMILES. 
## I was trying to find another way to study the problem. The first hurdle was the convert the Metabolites to Smiles.
## You can download these 2 files from https://figshare.com/articles/dataset/Qatar_Metabolomics_Study_on_Diabetes/5904022?file=10531345

library(purrr)
library(readxl)
library(tidyverse)

#Data 
metabolome=read_excel("/home/spirpinias/Desktop/QMDiab_metabolomics_Preprocessed.xlsx")
phenotypes=read_excel("/home/spirpinias/Desktop/QMDIAB/Urine/QMDiab_phenotypes.xlsx")

#Cleaning
phenotypes=phenotypes %>% remove_rownames() %>% column_to_rownames(var="QMDiab-ID") %>% as.data.frame()
metabolome=metabolome %>% remove_rownames() %>% column_to_rownames(var="QMDiab-ID") %>% as.data.frame() %>% select(-c(AGE,GENDER,BMI,ETH,T2D))

#Get Metabolite Name Vector
first=colnames(metabolome)

#Remove Unknowns
cleanMetabolites=first[grep("[a-z]",first)]

#Remove String in between Parenthesis
cleanMetabolites=gsub("[(][^)]*[)]","",cleanMetabolites)

#Remove String after last comma
cleanMetabolites=sub('^([^,]+,[^,]+).*', '\\1', cleanMetabolites)

#Remove String in Brackets
cleanMetabolites=sub("\\[.*]", "", cleanMetabolites)

#Delete the whitespace
cleanMetabolites=gsub(" ", "", cleanMetabolites, fixed = TRUE)

#Remove trailing Asterisk
cleanMetabolites=gsub("[*].*$","",cleanMetabolites)

#Contact OPSIN -- This is a website service that generates SMILES from metabolites/IUPAC chemical names.
#So I figured out a way to talk to it.
generateSMILES=function(j){
  url=download.file(url = paste0('http://opsin.ch.cam.ac.uk/opsin/',cleanMetabolites[j],".smi"),destfile = paste0(cleanMetabolites[j],'.smi'))
}

#Error Exception Handling to continue if OPSIN does not have information on your metabolite.
possiblySomeSMILES=possibly(generateSMILES,otherwise = "No Information")
possiblySomeInfo=map(1:length(cleanMetabolites),possiblySomeSMILES)
