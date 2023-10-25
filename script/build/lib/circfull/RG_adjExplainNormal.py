import pandas as pd, numpy as np,pyfasta,sys,os
import intervals as I # Note: use "pip install python-intervals" to install not intervals or pyinterval
from multiprocessing import Pool
hangLen=300

def createIntervals(s,e):
    x=I.empty()
    s=[int(i) for i in s.split(',')]
    e=[int(i) for i in e.split(',')]
    for i in range(len(s)):
        x=x|I.closed(s[i],e[i])
    return(x)

def outExon(exon):
    exon=list(exon)
    exonS=[]
    exonE=[]
    for i in exon:
        exonS.append(i.lower)
        exonE.append(i.upper)
        exonLen=sum([exonE[i]-exonS[i] for i in range(len(exonS))])
    return([exonS[0],exonE[-1],','.join([str(i) for i in exonS]),','.join([str(i) for i in exonE]),exonLen])

def getNewFL(i):
    idx=targetID_dict[i]
    tmp=FLdf_2.iloc[idx,:].copy()
    each1=tmp.iloc[0,:].copy()
    each2=tmp.iloc[1,:].copy()
    if abs(each1['reference_start']-each2['reference_start'])>hangLen:
        exon1=createIntervals(each1['exon_start'],each1['exon_end'])
        exon2=createIntervals(each2['exon_start'],each2['exon_end'])
        newExon=exon1 | exon2
        each1['reference_start'],each1['reference_end'],each1['exon_start'],each1['exon_end'],each1['exon_length']=outExon(newExon)
        each1['query_start']=min(each1['query_start'],each2['query_start'])
        each1['query_end']=max(each1['query_end'],each2['query_end'])
        each1['leftSeq']=genome.sequence({'chr': each1['chr'], 'start':each1['reference_start']-1, 'stop':each1['reference_start']}).upper()
        each1['rightSeq']=genome.sequence({'chr': each1['chr'], 'start':each1['reference_end']+1, 'stop':each1['reference_end']+2}).upper()
        return(each1.to_list())
    return(list(np.array([each1.to_list(),each2.to_list()]).ravel()))

def adjExplainNormal(genomeFile,outPrefix,thread,isSecond=False):
    global genome,FLdf_2, targetID_dict
    genome = pyfasta.Fasta(genomeFile)
    FLdf=pd.read_csv(outPrefix+"explainFL_Normal.txt",sep='\t')
    FLdf=FLdf.sort_values(by=["ID","query_start"],ascending=True)
    FLdf_counts=FLdf['ID'].value_counts()
    FLdf.index=FLdf.ID
    FLdf_2=FLdf.loc[FLdf_counts.index[FLdf_counts.values==2],:].copy()
    FLdf_no2=FLdf.loc[FLdf_counts.index[FLdf_counts.values!=2],:].copy()
    targetID=list(set(FLdf_2.index))
    targetID_dict={}

    for i in range(FLdf_2.shape[0]):
        if targetID_dict.__contains__(FLdf_2.index[i]):
            targetID_dict[FLdf_2.index[i]].append(i)
        else:
            targetID_dict[FLdf_2.index[i]]=[i]
    
    pool=Pool(processes=thread)
    result=pool.map(getNewFL,targetID)
    pool.close()
    pool.join()
    newFL_2=[]
    for i in result:
        for j in i:
            newFL_2.append(j)
    newFL_2=pd.DataFrame(np.array(newFL_2).reshape(-1,12))
    newFL_2.columns=FLdf_no2.columns
    result=pd.concat([newFL_2,FLdf_no2])
    result.to_csv(outPrefix+"explainFL_Normal_adj.txt",sep="\t",header=True,index=False)
    if isSecond:
        circOrigin=result.loc[:,['ID','strand']].copy()
        strandScore=[]
        for i in circOrigin.strand:
            if i=='+':
                strandScore.append(1)
            elif i=='-':
                strandScore.append(-1)
            else:
                strandScore.append(0)
        circOrigin.strand=strandScore
        circOrigin.columns=['ID','score']
        circOrigin=circOrigin.drop_duplicates('ID')
        return(circOrigin)

