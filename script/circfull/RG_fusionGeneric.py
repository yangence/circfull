import pandas as pd, numpy as np,pyfasta,sys
errorLen=40
def mysplit(x):
    return(len(str(x).split(',')))

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
    
def getExonLen(x):
    end=np.array([int(i) for i in x['exon_end'].split(',')])
    start=np.array([int(i) for i in x['exon_start'].split(',')])
    return(sum(end-start)+len(end))

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