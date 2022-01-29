###this is a python script to calculate the cell-cycle score, the input data is a expression matrix with each row for a gene and each column for a cell 
import sys
import os
import pandas as pd
import numpy as np
import scipy.stats as stats
import math
from sklearn import preprocessing

reference_folder = sys.argv[1]

df1=pd.read_csv(os.path.join(reference_folder, 'G1S.csv'),index_col=1,header=None)
df2=pd.read_csv(os.path.join(reference_folder, 'S.csv'),index_col=1,header=None)
df3=pd.read_csv(os.path.join(reference_folder, 'G2M.csv'),index_col=1,header=None)
df4=pd.read_csv(os.path.join(reference_folder, 'M.csv'),index_col=1,header=None)
df5=pd.read_csv(os.path.join(reference_folder, 'MG1.csv'),index_col=1,header=None)
df=pd.read_csv(sys.argv[2],index_col=0)
g1s=pd.concat([df1,df],axis=1,join='inner')
s=pd.concat([df2,df],axis=1,join='inner')
g2m=pd.concat([df3,df],axis=1,join='inner')
m=pd.concat([df4,df],axis=1,join='inner')
mg1=pd.concat([df5,df],axis=1,join='inner')

def phasescore(pds):
    a=[]
    pds2=pds.iloc[:,3:]
    mean=pds2.mean(axis=0)
    for i in range(pds2.shape[0]):
        a.append(stats.pearsonr(pds2.iloc[i,:],mean)[0])
        
    pds2['mean']=a
    pds3=pds2[pds2['mean']>0.3]
    pds4=pds3.iloc[:,:-1]+1.0
    pds4=np.log2(pds4)
    score=pds4.mean(axis=0)
    score=preprocessing.scale(score)
    return score

g1sscore=phasescore(g1s)
sscore=phasescore(s)
g2mscore=phasescore(g2m)
mscore=phasescore(m)
mg1score=phasescore(mg1)

toarrray=np.array([g1sscore,sscore,g2mscore,mscore,mg1score])

cellscore=pd.DataFrame(toarrray,columns=df.columns,index=['G1/S','S','G2/M','M','M/G1'])
cellscore.to_csv(sys.argv[3]+'/cellscore.csv')