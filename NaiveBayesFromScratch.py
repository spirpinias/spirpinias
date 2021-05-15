## This is a homework assignment I did during my course of Applied Machine Learning at Cornell Tech in the Fall of 2020.
## This course was taught by Volodymyr Kuleshov and it was fantastic. See his Youtube series for the exact course.
## We had to create a Naive Bayes Classifier from Scratch for competing in Kaggle Competition.
## This can likely be done a lot better but at the time I had about 6 months of real Python experience and possibly 3 months doing OOP, but I never gave up.

class NaiveBayes2(object):
    
    #Preprocessing
    def __init__(self,link,maxFeatures,minDF,maxDF):
        import numpy as np
        import csv
        from nltk.corpus import stopwords 
        from nltk.tokenize import word_tokenize
        import string
        import re
        import random
        from string import digits
        from sklearn.feature_extraction.text import CountVectorizer   
        
        self.link=link
        text,labels=[],[]
        stop_words=set(stopwords.words('english'))
        
        with open(self.link) as csv_file:
            csv_reader=csv.reader(csv_file,delimiter=',')
            next(csv_reader,None) #Skip the Header
            translator=str.maketrans('','',string.punctuation) #remove punctuation
            remove_digits=str.maketrans('','',digits) #remove digits
            for row in csv_reader:
                first=word_tokenize(re.sub(r"http\S+","",(row[3].lower()).translate(translator)))
                second=[]
                for w in first:
                    if w not in stop_words: 
                        second.append(re.sub(r'[^\x00-\x7f]',r'',(w.translate(remove_digits)))) #remove numerical digits/funkystuff
                text.append(' '.join(second))
                labels.append(row[4])
        
        self.training=text
        self.labels=labels
        
        #Vectorizing and Stuff
        self.vectorizer= CountVectorizer(binary=True, max_features=maxFeatures,min_df=minDF,max_df=maxDF)         
        
        #Fit to Training Data
        vector=self.vectorizer.fit_transform(self.training).toarray()    
        
        #Split into Classes~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        n=vector.shape[0]
        d=vector.shape[1]
        self.K=len(np.unique(self.labels))
        
        #These are the shapes of the parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        self.psis=np.zeros([self.K,d])
        self.phis=np.zeros([self.K])
        
        ## Dimensions are good. Proceed~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
        lblwords=np.array(self.labels,dtype=np.int)
        
        for k in range(self.K):
            
            X_k=vector[lblwords==k]
            self.psis[k]=np.mean(X_k,axis=0)
            self.phis[k]=X_k.shape[0]/float(n)
                        
#Make Predictions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def Predict(self,newX):
        import numpy as np
        
        test=newX
        future=[]
        
        for j in test:
            
            newData=self.vectorizer.transform([j]).toarray()
            n,d=newData.shape
            x=np.reshape(newData,(1,n,d))
            psis=np.reshape(self.psis,(self.K,1,d))
            psis=psis.clip(1e-14,1-1e-14)
            logpy=np.log(self.phis).reshape([self.K,1])
            logpxy=x*np.log(psis)+(1-x)*np.log(1-psis)
            logpyx=logpxy.sum(axis=2)+logpy
            futures,_= logpyx.argmax(axis=0).flatten(),logpyx.reshape([self.K,n])
            future.append(futures[0])
        
        return(np.array(future))     
        
    #Check Accuracy~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def accuracy(self,predictions,truth):
        
        correct=0
        incorrect=0
        #F1=np.mean(truth==predictions)
        
        data=zip(predictions,truth)

        for a in data:
            if a[0]==a[1]:
                correct+=1
            
        if correct == 0:
            return('Accuracy is 0%')

        else:
            return correct/float(len(truth))
          
          
          
   #Running it~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #Sigare is Greek for "Go Slow, man"
  sigare=NaiveBayes2('/home/spirpinias/Desktop/AMLHW2/train.csv',5000,0.01,0.75)
  sigare.Predict(targetX)
