# My friend Andrew wrote a book called Chronicle or Kaplan1.pdf.
# I wanted to practice NLP on this book so I naturally exercised simple tools to play around.
# This analysis is under construction as I am going to add more methods to look for information in the book.


import PyPDF2
from gensim.utils import simple_preprocess
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.decomposition import TruncatedSVD
import matplotlib.pyplot as plt


## Preprocessing
pdfileObj=open('/home/spirpinias/Downloads/Kaplan1.pdf','rb')
pdfReader=PyPDF2.PdfFileReader(pdfileObj)
theBook=[pdfReader.getPage(x).extractText() for x in range(8,pdfReader.numPages)]
theBook=[simple_preprocess(theBook[x]) for x in range(0,len(theBook))]
theFlatBook=[sentence for sub in theBook for sentence in sub]

## Create the Dictionary for BoW
corpus=[' '.join(ele) for ele in theBook]

## Modeling
vectorizer=CountVectorizer()
X=vectorizer.fit_transform(corpus)

clf = TruncatedSVD(100)
Xpca = clf.fit_transform(X)
clf.explained_variance_ratio_

## Nothing much Here. 
