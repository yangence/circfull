import pandas as pd, numpy as np,pyfasta,sys,os
from multiprocessing import Pool
errorLen=40
hangLen=300 # consistent with detectBS.py

def mysplit(x):
    return(len(x.split(',')))

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

def getCommonPosLeft(pos):    
    pos=pos.map(lambda x: [int(i) for i in x.split(',')])
    exonNum=len(pos.iloc[0])
    commonPos=[]
    for i in range(1,exonNum):
        commonPos.append(pos.map(lambda x:x[i]).value_counts().index[0]+1)
    return(commonPos)

def getCommonPosRight(pos):    
    pos=pos.map(lambda x: [int(i) for i in x.split(',')])
    exonNum=len(pos.iloc[0])
    commonPos=[]
    for i in range(exonNum-1):
        commonPos.append(pos.map(lambda x:x[i]).value_counts().index[0])
    return(commonPos)

def getCommonPosBoth(left,right):
    start=[]
    startFinal=int(left.map(lambda x: [i for i in x.split(',')][-1]).value_counts().index[0])
    end=[int(right.map(lambda x: [i for i in x.split(',')][0]).value_counts().index[0])]
    left=left.map(lambda x: [i for i in x.split(',')][1:-1])
    right=right.map(lambda x: [i for i in x.split(',')][1:-1])
    exonNum=len(left.iloc[0])
    both=[]
    for i in range(left.shape[0]):
        tmpLeft=left.iloc[i]
        tmpRight=right.iloc[i]
        tmpBoth=[]
        for j in range(exonNum):
            tmpBoth.append(tmpLeft[j]+'_'+tmpRight[j])
        both.append(tmpBoth)
    for i in range(exonNum):
        commonPos=pd.Series(both).map(lambda x:x[i]).value_counts().index[0].split('_')
        start.append(int(commonPos[0])+1)
        end.append(int(commonPos[1]))
    start.append(startFinal+1)
    return([start,end])

def getInterExon(i):
    id=seqID[i]
    tmpBS=BS_Normal_adj.loc[id]
    tmpFL=FL_Normal.loc[[id]].copy()
    #tmpFL=tmpFL.iloc[range(1,tmpFL.shape[0]-1)]
    tmpID=list(set(np.where(abs(tmpFL['reference_start'] - tmpBS['start']) <hangLen)[0]) & set(np.where(abs(tmpFL['reference_end'] - tmpBS['end']) <hangLen)[0]))
    tmpFL=tmpFL.iloc[tmpID,:]
    tmpChr=tmpBS['chr']
    if tmpFL.shape[0]<1:
        return([id,tmpChr,'','',False])
    tmpExonNum=tmpFL['exonNum']
    tmpExonNumVC=tmpExonNum.value_counts()
    tmpExonNumVC_value=tmpExonNumVC.values
    #if len(tmpExonNumVC_value)>1:  # This option may lost some true position, mostly, max exonNum is good
    #    if tmpExonNumVC_value[0]==tmpExonNumVC_value[1]:
    #        return([id,tmpChr,'','',False])
    tmpFL=tmpFL[tmpExonNum==max(tmpExonNumVC.index)]
    
    #tmpFL=tmpFL[tmpExonNum==tmpExonNumVC.index[0]]
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

def getExonLen(x):
    end=np.array([int(i) for i in x['exon_end'].split(',')])
    start=np.array([int(i) for i in x['exon_start'].split(',')])
    return(sum(end-start)+len(end))

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
    
def selfCluster(x):
    x.index=range(x.shape[0])
    nonCan_count=x['circID'].value_counts()
    nonCan_countDf=pd.DataFrame({'circID':list(nonCan_count.index),'num':list(nonCan_count.values)})
    nonCan_countDf['chr']=pd.DataFrame(nonCan_countDf['circID'].map(lambda x: x.split('|')[0]))
    nonCan_countDf['start']=pd.DataFrame(nonCan_countDf['circID'].map(lambda x: int(x.split('|')[1])))
    nonCan_countDf['end']=pd.DataFrame(nonCan_countDf['circID'].map(lambda x: int(x.split('|')[2])))
    nonCan_countDf=nonCan_countDf.sort_values(['chr','start','end'])
    nonCan_numDict={}
    for i in range(nonCan_countDf.shape[0]):
        tmp=nonCan_countDf.iloc[i,].copy()
        circID=tmp['circID']
        num=tmp['num']
        nonCan_numDict[circID]=num
    nonCan_dict={}
    for i in range(nonCan_countDf.shape[0]):
        tmp=nonCan_countDf.iloc[i,]
        chr1=tmp['chr']
        start1=tmp['start']
        end1=tmp['end']
        circID1=tmp['circID']
        num1=nonCan_numDict[circID1]
        num0=nonCan_numDict[circID1]
        targetID=circID1
        targetJ=i
        adj=False
        for j in range(max(0,i-20),min(i+20,nonCan_countDf.shape[0])):
            tmp=nonCan_countDf.iloc[j,]
            chr2=tmp['chr']
            start2=tmp['start']
            end2=tmp['end']
            circID2=tmp['circID']
            num2=nonCan_numDict[circID2]
            if chr1==chr2 and abs(start1-start2)<errorLen and abs(end1-end2)<errorLen:
                adj=True
                if num1<=num2:
                    targetID=circID2
                    num1=num2
                    targetJ=j
            else:
                if adj:
                    break
        if targetID != circID1:
            nonCan_numDict[targetID]+=num0
            nonCan_numDict[circID1]=0
        nonCan_dict[circID1]=targetID
    x.index=x['circID']
    for i in nonCan_dict.keys():
        value=nonCan_dict[i]
        if i !=value:
            target=x.loc[value].copy()
            if len(target.shape)>1:
                target=target.iloc[0]
            x.loc[i,'start']=target['start']
            x.loc[i,'end']=target['end']
            x.loc[i,'leftSeq']=target['leftSeq']
            x.loc[i,'rightSeq']=target['rightSeq']
            x.loc[i,'motif']=target['motif']
            
    x.index=x['ID']
    x['circID']=x['chr']+'|'+x['start'].map(str)+'|'+x['end'].map(str)
    return(x)
########## End of functions ##########

def constructFL(options):
    global isSecond,genome,genome,BS_Normal_adj,BS_Normal_adj_no,ExonSdict,ExonEdict,seqID,FL_Normal
    isSecond=False
    genomeFile=options[0]
    outPrefix=options[1]
    thread=int(options[2])
    if len(options)>3:
        strandFile=options[3]
        strandFile.index=strandFile['ID']
    genome = pyfasta.Fasta(genomeFile)
    # start site 0 based
    FLFIle=outPrefix+"explainFL_Normal_adj.txt"
    FL_Normal=pd.read_csv(FLFIle,sep='\t')
    #FLdf_NC=pd.read_csv(outPrefix+"explainFL_NC.txt",sep='\t')
    #FL_Normal=pd.concat([FL_Normal,FLdf_NC]).copy() # include candidate to increase sensitive
    FL_Normal=FL_Normal.sort_values(['ID','query_start']).set_index('ID')

    # start site 1 based
    BSadjFile=outPrefix+"BS_Normal_adj.txt"
    BS_Normal_adj=pd.read_csv(BSadjFile,sep="\t",names=['ID','chr','start','end','leftSeq','rightSeq','leftAdj','rightAdj'])

    BS_Normal_adj['circID']=BS_Normal_adj['chr']+'|'+BS_Normal_adj['start'].map(str)+'|'+BS_Normal_adj['end'].map(str)
    BS_Normal_adj['motif']=BS_Normal_adj['leftSeq']+BS_Normal_adj['rightSeq']
    BS_Normal_adj.index=BS_Normal_adj['circID']


    BS_Normal_adj_can=BS_Normal_adj.loc[[i in ['ACCT','AGGT'] for i in BS_Normal_adj['motif'] ]].copy()
    #print(BS_Normal_adj_can.shape)
    if isSecond:
        BS_Normal_adj_can['score']=strandFile.loc[BS_Normal_adj_can['ID'],'score'].values
        passIdx=[]
        for i in range(BS_Normal_adj_can.shape[0]):
            tmp=BS_Normal_adj_can.iloc[i]
            tmpScore=tmp['score']
            tmpMotif=tmp['motif']
            if (tmpScore>0 and tmpMotif == 'AGGT') or (tmpScore<0 and tmpMotif == 'ACCT'):
                passIdx.append(i)
        BS_Normal_adj_can=BS_Normal_adj_can.iloc[passIdx].copy()
    #print(BS_Normal_adj_can.shape)
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


        
    if BS_Normal_adj_nonCan.shape[0]>0:
        pool=Pool(processes=thread)
        BS_Normal_adj_nonCan=pd.concat(pool.map(selfCluster,[BS_Normal_adj_nonCan[BS_Normal_adj_nonCan['chr']==i].copy() for i in list(set(BS_Normal_adj_nonCan['chr']))]))
        pool.close()
        pool.join()
    if BS_Normal_adj_Can.shape[0]>0:
        pool=Pool(processes=thread)
        BS_Normal_adj_Can=pd.concat(pool.map(selfCluster,[BS_Normal_adj_Can[BS_Normal_adj_Can['chr']==i].copy() for i in list(set(BS_Normal_adj_Can['chr']))]))
        pool.close()
        pool.join()
    BS_Normal_adj=pd.concat([BS_Normal_adj_nonCan,BS_Normal_adj_Can])
    pool=Pool(processes=thread)
    BS_Normal_adj=pd.concat(pool.map(selfCluster,[BS_Normal_adj[BS_Normal_adj['chr']==i].copy() for i in list(set(BS_Normal_adj['chr']))]))
    pool.close()
    pool.join()
    BS_Normal_adj=BS_Normal_adj.set_index('ID')
    seqID=list(BS_Normal_adj.index)
    FL_Normal['exonNum']=FL_Normal['exon_start'].map(mysplit)

    pool=Pool(processes=thread)
    FL_con=pool.map(getInterExon,[i for i in range(BS_Normal_adj.shape[0])])
    pool.close()
    pool.join()
    FL_con_df=pd.DataFrame(FL_con,columns=['ID','chr','exon_start','exon_end','pass'])
    FL_con_df=FL_con_df.set_index('ID')
    BS_Normal_adj_con=BS_Normal_adj.loc[FL_con_df.index][['start','end','leftSeq','rightSeq','leftAdj','rightAdj']].copy()
    FL_result=pd.concat([FL_con_df,BS_Normal_adj_con],axis=1)
    FL_result['exon_leftSeq']=FL_result.apply(getLeftSeq,axis=1)
    FL_result['exon_rightSeq']=FL_result.apply(getRightSeq,axis=1)
    FL_result['circID']=FL_result['chr']+'|'+FL_result['start'].map(str)+'|'+FL_result['end'].map(str)
    FL_result=FL_result[FL_result['start']!=0]
    FL_result=FL_result[FL_result['end']!=0]
    FL_result=FL_result[FL_result['exon_start']!='']
    FL_result['len']=FL_result.apply(getExonLen,axis=1)
    FL_result['motif']=FL_result.apply(lambda x:str(x['leftSeq'])+str(x['rightSeq']),axis=1)
    FL_result['exonNum']=FL_result.apply(lambda x:len(x['exon_start'].split(',')),axis=1)
    FL_result=FL_result[['circID','chr','start','end','len','exonNum','exon_start','exon_end','motif','leftSeq','rightSeq','exon_leftSeq','exon_rightSeq','leftAdj','rightAdj']]
    FL_result=FL_result.sort_values(['circID','len'])
    FL_result.to_csv(outPrefix+"constructFL_Normal.txt",sep="\t")
