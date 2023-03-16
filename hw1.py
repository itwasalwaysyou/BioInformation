import csv
import numpy as np
from collections import Counter
import pandas
import math
import matplotlib.pyplot as plt
from sklearn.preprocessing import FunctionTransformer
from collections import Counter
import pylab as pl
from scipy import stats

from scipy.stats import wilcoxon
import scipy.stats as stats
from matplotlib.ticker import StrMethodFormatter
import operator
from scipy.stats import ranksums
from collections import OrderedDict



def question2():
    with open("HW1-GSE62944-count.csv",'r') as file: 



        reader = csv.reader(file)

        countPatient=len(next(reader))-1
        file.seek(0)  # point to the beginning of file
        
        next(reader)  # skip the first row
        countGen=0
        
        for row in reader:
            countGen+=1

        print("-------------------------question(2) part------------------------------")


        print(str(countGen)+" genes are included in the dataset, and "+str(countPatient)+" patient samples are in the dataset") 

    with open("HW1-GSE62944-clinical.csv",'r') as file2: 
            reader = csv.DictReader(file2)
            countShort= 0

            for row in reader:
                if row['Group']=="short":
                    countShort+=1
                
            countLong= countPatient-countShort

            

            print(str(countShort)+" patients survived short-(<=1 year) and "+str(countLong)+ " patients survived long-term(>1 year)")


def question3():
        with open("HW1-GSE62944-count.csv",'r') as file: 
            reader = csv.DictReader(file)
                        
            dic={}

            
            # first row data --- g(ij) : 
            for row in reader:
                dic=row
                #pop(0,0)
                dic.pop('')
                break
            
            
            
            # first raw (key---sample, value---#gene)
            dic= dict((k,int(v)) for k, v in dic.items())
            
        


            # total number of reads mapping to all genes in sample i : 
            for row in reader:
                # pop the gene name part
                row.pop('')
                
                # char --> interger for each value
                row= dict((k,int(v)) for k, v in row.items())
    
                dic=Counter(dic)+Counter(row)
                



            med=43114673.5
            
            
            #read again at the beging of file
        with open("HW1-GSE62944-count.csv",'r') as file: 
            reader_again = csv.DictReader(file)
        
            
            dic_copy={}

            
            # first row data --- g(ij) : 
            for row in reader_again:

                dic_copy=row
                #pop (0,0)
                dic_copy.pop('')
                break
            
            # first raw (key---sample, value---#gene)

            dic_copy= dict((k,int(v)*(med/dic[k])) for k, v in dic_copy.items())
            
            
            # sample_names= [i for i in dic_copy.value]
            
        

            
            
            # total number of reads mapping to all genes in sample i : 
            for row in reader_again:
                # pop the gene name part
                row.pop('')
                
                # char --> interger for each normalization value g'(ij)
                row= dict((k,int(v)*(med/dic[k])) for k, v in row.items())
    
                
                dic_copy=Counter(dic_copy)+Counter(row)
            
            print(dic)



sample_names=[]
genes=[]


def getSampleList():
    sample_names=[]
    with open("HW1-GSE62944-count.csv",'r') as file: 
        reader = csv.DictReader(file)
    
        
        dic_copy={}

        
        # first row data --- g(ij) : 
        for row in reader:

            dic_copy=row
            #pop (0,0)
            dic_copy.pop('')
            break
        
        # first raw (key---sample, value---#gene)
        sample_names= [i for i in dic_copy.keys()]
        return sample_names    



def log():
    sample_names=getSampleList()  
    

    
    
    
    df= pandas.read_csv("HW1-GSE62944-count.csv")

    
    med=43114673.5
    for curr in sample_names:
        # -- normalized count for gene j in sample i
        df[curr] = df[curr]* (med/df[curr].sum())
        
        # -- log-transform for gene j in sample i
        df[curr] = df[curr].apply(np.log1p).round(decimals=4)
            
    return df
    
    




def quantile():  
    
    df= log()
    sample_names=getSampleList()

    # sort 
    df_copy = df.apply(lambda x: x.sort_values(ascending=False).values)
    
    # add mean column in df 
    df['mean']= df_copy.mean(axis=1)
    # .apply("{:.4f}".format) 
    df_copy['mean']= df_copy.mean(axis=1)
    
    # .apply("{:.4f}".format) 
    
    
    # 3(e) - for first five samples 
    FIVE_SAMPLES=["Unnamed: 0","TCGA-02-2483-01A-01R-1849-01","TCGA-02-2485-01A-01R-1849-01","TCGA-02-2486-01A-01R-1849-01","TCGA-06-0129-01A-01R-1849-01","TCGA-06-0178-01A-01R-1849-01"]
 


    
    FIVE_SAMPLES2=["TCGA-02-2483-01A-01R-1849-01","TCGA-02-2485-01A-01R-1849-01","TCGA-02-2486-01A-01R-1849-01","TCGA-06-0129-01A-01R-1849-01","TCGA-06-0178-01A-01R-1849-01"]


    # dic_curr={}
    count=1;
    

    print(df)

    for curr in sample_names:
        
        dic_curr= {row[0]:row[count] for row in df.values}
        
        dic_curr=dict(sorted(dic_curr.items(), key=lambda item: item[1], reverse=True))  

        
        df[curr]= pandas.Series(dic_curr.keys())
                
        dic_curr= {row[count]:row[155] for row in df.values}       

        
        df[curr]=df['Unnamed: 0']   
        df[curr]=df[curr].map(dic_curr)
        count+=1
    


#-----------------------------------------------------
    

    print(df)
    
    short_df=df.loc[:,'Unnamed: 0':"TCGA-76-4932-01A-01R-1850-01"]
    long_df=df.loc[:,"TCGA-02-0047-01A-01R-1849-01":"TCGA-76-4931-01A-01R-1850-01"]
    
    p = []

    genes=df["Unnamed: 0"].tolist()


    for i in range(23368): # iterate gene index
        
        data = ranksums(short_df.loc[i, "TCGA-02-2483-01A-01R-1849-01":].values.astype(float).tolist(), long_df.loc[i, "TCGA-02-0047-01A-01R-1849-01":].values.astype(float).tolist())
        p.append({"gene" :genes[i], "p" : data.pvalue.astype(float), "stat" : data.statistic })
            # df.iloc[[i],[0]]})
        
    wilc_df = pandas.DataFrame(p)
    wilc_df = wilc_df.set_index("gene")


    print(wilc_df)
    
    top10=wilc_df.nsmallest(10, "p")
    print(top10)
    

    
    DESeq_df= pandas.read_csv("HW1-DESeq2.csv")

    
    
    
    DESeq_df = DESeq_df.set_index("Unnamed: 0")
    
    list1=(wilc_df[wilc_df["p"] < 0.05].index).tolist()

    print(len(list1))
    
    list2=(DESeq_df[DESeq_df["pvalue"] <0.05].index).tolist()
    
    print(len(list2))

    print(len(set(list1)&set(list2)))



if __name__ == '__main__':

    quantile()