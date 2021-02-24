import pandas as pd, numpy as np,pyfasta,sys
from .RG_fusionGeneric import *
from multiprocessing import Pool
errorLen=40

def getLeftSeq(x):
    if len(x['exon_start'])==0:
        return('')
    chr=x['chr']
    x=[int(i) for i in x['exon_start'].split(',')]
    y=[]
    for i in x:
        y.append(genome.sequence({'chr': chr, 'start':i-2, 'stop':i-1}).upper())
    return(','.join(y))

def getRightSeq(x):
    if len(x['exon_end'])==0:
        return('')
    chr=x['chr']
    x=[int(i) for i in x['exon_end'].split(',')]
    y=[]
    for i in x:
        y.append(genome.sequence({'chr': chr, 'start':i+1, 'stop':i+2}).upper())
    return(','.join(y))
def getInterExon(i):
    id=seqID[i]
    tmpBS=BS_Normal_adj.loc[id]
    tmpFL=FL_Normal.loc[id]
    #tmpFL=tmpFL.iloc[range(1,tmpFL.shape[0]-1)]
    maxError=100
    tmpFL=tmpFL.iloc[list(set(np.where(abs(tmpFL['reference_start'] - tmpBS['start']) <maxError)[0]) & set(np.where(abs(tmpFL['reference_end'] - tmpBS['end']) <maxError)[0]))]
    tmpChr=tmpBS['chr']
    if tmpFL.shape[0]<1:
        return([id,tmpChr,'','',False])
    tmpExonNum=tmpFL['exonNum']
    tmpExonNumVC=tmpExonNum.value_counts()
    tmpExonNumVC_value=tmpExonNumVC.values
    if len(tmpExonNumVC_value)>1:
        if tmpExonNumVC_value[0]==tmpExonNumVC_value[1]:
            return([id,tmpChr,'','',False])
    tmpFL=tmpFL[tmpExonNum==tmpExonNumVC.index[0]]
    tmpFL_maxLen=tmpFL[tmpFL['exon_length']==max(tmpFL['exon_length'])]
    tmpFL_maxLen_exonNum=list(set(tmpFL_maxLen['exonNum']))
    if 1 in tmpFL_maxLen_exonNum:
            tmpStart=[tmpBS['start']]
            tmpEnd=[tmpBS['end']]
    elif 2 in tmpFL_maxLen_exonNum:
        tmpStart=[tmpBS['start']]
        tmpEnd=getCommonPosRight(tmpFL['exon_end'])
        tmpStart.extend(getCommonPosLeft(tmpFL['exon_start']))
        tmpEnd.extend([tmpBS['end']])
    else:
        # transform to 0 based
        tmpStart=[tmpBS['start']]
        #tmpEnd=getCommonPosRight(tmpFL['exon_end'])
        #tmpStart.extend(getCommonPosLeft(tmpFL['exon_start']))
        commonPos=getCommonPosBoth(tmpFL['exon_start'],tmpFL['exon_end'])
        tmpStart.extend(commonPos[0])
        tmpEnd=commonPos[1]
        tmpEnd.extend([tmpBS['end']])
    if len(np.where((np.array(tmpEnd)-np.array(tmpStart))<0)[0])>0:
        return([id,tmpChr,'','',False])
    return([id,tmpChr,','.join([str(i) for i in tmpStart]),','.join([str(i) for i in tmpEnd]),True])



def getComBSAdj(i):
    tmp=BS_Normal_adj_no.iloc[i].copy()
    tmpChr=tmp['chr']
    tmpStart=tmp['start']
    tmpEnd=tmp['end']
    startKey=tmpChr+'_'+str(tmpStart)
    endKey=tmpChr+'_'+str(tmpEnd)
    if ExonSdict.__contains__(startKey) and ExonEdict.__contains__(endKey):
        cValue=list(set(ExonSdict[startKey]) & set(ExonEdict[endKey]))
        cValueLen=len(cValue)
        if cValueLen==0:
            return(tmp)
        elif cValueLen==1:
            keyValue=cValue[0].split('|')
            tmp['start']=int(keyValue[1])
            tmp['end']=int(keyValue[2])
            tmp['leftSeq']=genome.sequence({'chr': tmpChr, 'start':tmp['start']-2, 'stop':tmp['start']-1}).upper()
            tmp['rightSeq']=genome.sequence({'chr': tmpChr, 'start':tmp['end']+1, 'stop':tmp['end']+2}).upper()
        else:
            keyValue=[i.split('|')  for i in cValue]
            minValue=float('Inf')
            minIndex=0
            for j in range(len(keyValue)):
                eachKeyValue=keyValue[j]
                eachDiff=abs(tmpStart-int(eachKeyValue[1]))+abs(tmpEnd-int(eachKeyValue[2]))
                if eachDiff < minValue:
                    minValue=eachDiff
                    minIndex=j
            tmp['start']=int(keyValue[minIndex][1])
            tmp['end']=int(keyValue[minIndex][2])
            tmp['leftSeq']=genome.sequence({'chr': tmpChr, 'start':tmp['start']-2, 'stop':tmp['start']-1}).upper()
            tmp['rightSeq']=genome.sequence({'chr': tmpChr, 'start':tmp['end']+1, 'stop':tmp['end']+2}).upper()
    return(tmp)

def fusion1ConstructFL(genomeFile,outPrefix,adjFusionNum,thread):
    global genome,BS_Normal_adj,FL_Normal,ExonSdict,ExonEdict,BS_Normal_adj_no,seqID
    genome = pyfasta.Fasta(genomeFile)
    # start site 0 based
    FLFIle=outPrefix+"explainFL_Fusion1.txt"
    FL_Normal=pd.read_csv(FLFIle,sep='\t')
    FL_Normal=FL_Normal.sort_values(['ID','query_start']).set_index('ID')
    # start site 1 based
    BSadjFile=outPrefix+"BS_Fusion1_adj_"+adjFusionNum+".txt"
    BS_Normal_adj=pd.read_csv(BSadjFile,sep="\t",names=['ID','chr','start','end','leftSeq','rightSeq'])
    BS_Normal_adj['circID']=BS_Normal_adj['chr']+'|'+BS_Normal_adj['start'].map(str)+'|'+BS_Normal_adj['end'].map(str)
    BS_Normal_adj['motif']=BS_Normal_adj['leftSeq']+BS_Normal_adj['rightSeq']
    BS_Normal_adj.index=BS_Normal_adj['circID']


    BS_Normal_adj_can=BS_Normal_adj.loc[[i in ['ACCT','AGGT'] for i in BS_Normal_adj['motif'] ]].copy()

    BS_Normal_adj_can_counts=BS_Normal_adj_can['circID'].value_counts()
    BS_Normal_adj_can_index=BS_Normal_adj_can_counts.loc[BS_Normal_adj_can_counts>1].index
    BS_Normal_adj_temp=BS_Normal_adj.loc[BS_Normal_adj_can_index]
    BS_Normal_adj_temp_nodup=BS_Normal_adj_temp.drop_duplicates('circID')
    BS_Normal_adj_no=BS_Normal_adj.loc[list(set(BS_Normal_adj.index) - set(BS_Normal_adj_can_index))]


    ExonSdict={}
    ExonEdict={}
    for i in range(BS_Normal_adj_temp_nodup.shape[0]):
        tmp=BS_Normal_adj_temp_nodup.iloc[i]
        tmpChr=tmp['chr']
        tmpStart=tmp['start']
        tmpEnd=tmp['end']
        tmpValue=tmp.name
        for j in range(errorLen*2+1):
            startKey=tmpChr+'_'+str(tmpStart-errorLen+j)
            endKey=tmpChr+'_'+str(tmpEnd-errorLen+j)
            if ExonSdict.__contains__(startKey):
                ExonSdict[startKey].append(tmpValue)
            else:
                ExonSdict[startKey]=[tmpValue]
            if ExonEdict.__contains__(endKey):
                ExonEdict[endKey].append(tmpValue)
            else:
                ExonEdict[endKey]=[tmpValue]
                
    BS_Normal_adj_no=pd.DataFrame([getComBSAdj(i) for i in range(BS_Normal_adj_no.shape[0])])
    BS_Normal_adj=pd.concat([BS_Normal_adj_temp,BS_Normal_adj_no])
    BS_Normal_adj['motif']=BS_Normal_adj['leftSeq']+BS_Normal_adj['rightSeq']

    BS_Normal_adj_Can=BS_Normal_adj.loc[[i in ['ACCT','AGGT'] for i in BS_Normal_adj['motif'] ]].copy()
    BS_Normal_adj_nonCan=BS_Normal_adj.loc[[i not in ['ACCT','AGGT'] for i in BS_Normal_adj['motif'] ]].copy()
    chrNum_nonCan=list(set(BS_Normal_adj_nonCan['chr']))
    if len(chrNum_nonCan)>1:
        pool=Pool(processes=thread)
        BS_Normal_adj_nonCan=pd.concat(pool.map(selfCluster,[BS_Normal_adj_nonCan[BS_Normal_adj_nonCan['chr']==i].copy() for i in chrNum_nonCan]))
        pool.close()
        pool.join()
    elif len(chrNum_nonCan)==1:
        BS_Normal_adj_nonCan=selfCluster(BS_Normal_adj_nonCan.copy())
        
    chrNum_Can=list(set(BS_Normal_adj_Can['chr']))
    if len(chrNum_Can)>1:
        pool=Pool(processes=thread)
        BS_Normal_adj_Can=pd.concat(pool.map(selfCluster,[BS_Normal_adj_Can[BS_Normal_adj_Can['chr']==i].copy() for i in chrNum_Can]))
        pool.close()
        pool.join()
    elif len(chrNum_Can)==1:
        BS_Normal_adj_Can=selfCluster(BS_Normal_adj_Can.copy())

    BS_Normal_adj=pd.concat([BS_Normal_adj_nonCan,BS_Normal_adj_Can])
    chrNum_all=list(set(BS_Normal_adj['chr']))
    if len(chrNum_all)>1:
        pool=Pool(processes=thread)
        BS_Normal_adj=pd.concat(pool.map(selfCluster,[BS_Normal_adj[BS_Normal_adj['chr']==i].copy() for i in chrNum_all]))
        pool.close()
        pool.join()
    else:
        BS_Normal_adj=selfCluster(BS_Normal_adj.copy())
    BS_Normal_adj=BS_Normal_adj.set_index('ID')
    seqID=list(BS_Normal_adj.index)
    FL_Normal['exonNum']=FL_Normal['exon_start'].map(mysplit)
    pool=Pool(processes=thread)
    FL_con=pool.map(getInterExon,[i for i in range(BS_Normal_adj.shape[0])])
    pool.close()
    pool.join()
    FL_con_df=pd.DataFrame(FL_con,columns=['ID','chr','exon_start','exon_end','pass'])
    if FL_con_df.shape[0]==1 and (not FL_con_df.iloc[0,4]):
        return(False)
    FL_con_df=FL_con_df.set_index('ID')
    BS_Normal_adj_con=BS_Normal_adj.loc[FL_con_df.index][['start','end','leftSeq','rightSeq']].copy()
    FL_result=pd.concat([FL_con_df,BS_Normal_adj_con],axis=1)
    FL_result['exon_leftSeq']=FL_result.apply(getLeftSeq,axis=1)
    FL_result['exon_rightSeq']=FL_result.apply(getRightSeq,axis=1)
    FL_result['circID']=FL_result['chr']+'|'+FL_result['start'].map(str)+'|'+FL_result['end'].map(str)
    FL_result=FL_result[FL_result['start']!=0]
    FL_result=FL_result[FL_result['end']!=0]
    FL_result=FL_result[FL_result['exon_start']!='']
    FL_result['len']=FL_result.apply(getExonLen,axis=1)
    FL_result['motif']=FL_result.apply(lambda x:x['leftSeq']+x['rightSeq'],axis=1)
    FL_result['exonNum']=FL_result.apply(lambda x:len(x['exon_start'].split(',')),axis=1)
    FL_result=FL_result[['circID','chr','start','end','len','exonNum','exon_start','exon_end','motif','leftSeq','rightSeq','exon_leftSeq','exon_rightSeq']]
    FL_result=FL_result.sort_values(['circID','len'])
    FL_result.to_csv(outPrefix+"constructFL_Fusion1_"+adjFusionNum+".txt",sep="\t")
    return(True)

