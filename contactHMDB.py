## This is a program that I wrote in the Laboratory of Jan Krumsiek for converting Metabolites into Human Metabolome Disease ID's.
## I used a web scraper to submit a form to the website metaboanylst.

from selenium import webdriver
import pandas as pd

#Initiate the Driver
browser=webdriver.Chrome()
browser.get('https://www.metaboanalyst.ca/faces/upload/ConvertView.xhtml')

#Send the Driver to Submit the Form
chooseFile=browser.find_element_by_id('form1:j_idt48').send_keys("/home/spirpinias/Desktop/myTest.dat")
submitFile=browser.find_element_by_id('form1:j_idt61').click()
elems=browser.find_element_by_xpath("//*[@id='form']/table[2]/tbody/tr[2]/td/a").get_attribute("href")
    
#Download the CSV and Create Table
listofIDS=pd.read_csv(elems)
