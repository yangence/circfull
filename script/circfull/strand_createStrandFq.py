import numpy as np,pandas as pd,sys
def createStrandFq(RG_dir,outPrefix,fastqFile):
    predictStrandDf=pd.read_csv(outPrefix+'strandProbability.txt',sep='\t',names=['ID','p'])
    pSdict=dict(zip(predictStrandDf.iloc[:,0],predictStrandDf.iloc[:,1]))
    oldFq=open(fastqFile)
    newFq=open(outPrefix+'strandedFq.fastq','w')

    while True:
        fq1=oldFq.readline()
        if fq1:
            fq2=oldFq.readline().strip()
            fq3=oldFq.readline()
            fq4=oldFq.readline().strip()
            ID=fq1.split()[0].split('_')[0][1:]
            if pSdict.__contains__(ID):
                pro=pSdict[ID]
                if pro<0.5:
                    fq2=fq2[::-1].replace('A','t').replace('T','a').replace('C','g').replace('G','c').upper()
                    fq4=fq4[::-1]
                newFq.write(fq1+fq2+'\n'+fq3+fq4+'\n')   
        else:
            break

    oldFq.close()
    newFq.close()
