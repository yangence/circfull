import pandas as pd, numpy as np,sys,os
import intervals as I
from math import log

def getNonStrandNum(tmpCluster):
    n=len(tmpCluster)
    for j in  range(len(tmpCluster)):
        eachID=tmpCluster[j]
        if eachReadStrand.__contains__(eachID):
            n+=-1
    return(n)

def getClusterStrand(tmpCluster):
    n=getNonStrandNum(tmpCluster)
    if n==0:
        return()
    for m in range(n):
        for j in  range(len(tmpCluster)):
            eachID=tmpCluster[j]
            global eachReadStrand
            if eachReadStrand.__contains__(eachID):
                continue
            canID=cluster2ID[ID2cluster[eachID]]
            for q in range(len(canID)):
                eachCanID=canID[q]
                if eachReadStrand.__contains__(eachCanID):
                    idxStrand=eachReadStrand[eachCanID]
                    tmpKey=eachID+'_'+eachCanID
                    tmpOK=False
                    if strandDict.__contains__(tmpKey):
                        tmpStrand=strandDictTrans[strandDict[tmpKey]]
                        tmpOK=True
                    elif strandDict.__contains__(eachCanID+'_'+eachID):
                        tmpStrand=strandDictTrans[strandDict[eachCanID+'_'+eachID]]
                        tmpOK=True
                    else:
                        continue
                    if tmpOK:
                        if idxStrand  == 0:
                            if tmpStrand == 0:
                                eachReadStrand[eachID]=1
                            else:
                                eachReadStrand[eachID]=0               
                        else:
                            if tmpStrand == 0:
                                eachReadStrand[eachID]=0
                            else:
                                eachReadStrand[eachID]=1
                    break
                else:
                    continue
def novoCluster(outPrefix):
    global eachReadStrand,cluster2ID,ID2cluster,strandDict,strandDictTrans
    thFile_Pass_HQ=pd.read_csv(outPrefix+'TideHunter_Pass.tab',sep='\t')
    targetID=thFile_Pass_HQ['ID']
    thFile_Pass_HQ=thFile_Pass_HQ.set_index('ID')
    raw2rawFile=outPrefix+'raw2raw.sort.paf'
    if os.path.exists(outPrefix+'strandProbability.txt'): # just for easy test of simulation data without strandProbability.txt file
        strandProbability=pd.read_csv(outPrefix+'strandProbability.txt',sep='\t',names=['ID','pro'])
        strandProbability['strand']=[round(i) for i in strandProbability['pro']]
        strandProbability=strandProbability.set_index('ID').sort_index()
        strandProbability_ID=list(strandProbability.index)

    clusterNumMin=2
    rawCoverage=0.5
    idDict={}
    for i in targetID:
        idDict[i]=[]
    rawDiff1=I.empty()
    rawDiff2=I.empty()
    raw2raw=open(raw2rawFile)
    line=raw2raw.readline().split('\t')
    while True:
        line[0]=line[0]
        line[1]=int(line[1])
        line[2]=int(line[2])
        line[3]=int(line[3])
        line[6]=int(line[6])
        line[7]=int(line[7])
        line[8]=int(line[8])
        rawDiff1=rawDiff1|I.closedopen(line[2],line[3])
        rawDiff2=rawDiff2|I.closedopen(line[7],line[8])
        newline=raw2raw.readline()
        if newline:
            newline=newline.split('\t')
            if newline[5] == line[5] and newline[0] == line[0]:
                line=newline
            else:
                rawDiff1=list(rawDiff1)
                rawDiff2=list(rawDiff2)
                diffLen1=0
                diffLen2=0
                for i in rawDiff1:
                    diffLen1+=(i.upper-i.lower)
                for i in rawDiff2:
                    diffLen2+=(i.upper-i.lower)
                if diffLen1/line[1]>rawCoverage and diffLen2/line[6]>rawCoverage :
                    idDict[line[0]].append(line[5])
                rawDiff1=I.empty()
                rawDiff2=I.empty()
                line=newline
        else:
            rawDiff1=list(rawDiff1)
            rawDiff2=list(rawDiff2)
            diffLen1=0
            diffLen2=0
            for i in rawDiff1:
                diffLen1+=(i.upper-i.lower)
            for i in rawDiff2:
                diffLen2+=(i.upper-i.lower)
            if diffLen1/line[1]>rawCoverage and diffLen2/line[6]>rawCoverage :
                idDict[line[0]].append(line[5])
            break
    raw2raw.close()
    idDictLen=[len(idDict[i]) for i in idDict.keys()]
    ID2cluster={}
    clusterList={}
    n=0
    for i in idDict:
        if len(idDict[i])>0:
            tmpCluster=list(set([i]+idDict[i]))
            previousCluster=[]
            for j in tmpCluster:
                if ID2cluster.__contains__(j):
                    previousCluster.append(ID2cluster[j])
            previousCluster=list(set(previousCluster))
            if len(previousCluster) == 0:
                clusterName='C'+str(n)
                clusterList[clusterName]=tmpCluster
                for j in tmpCluster:
                    ID2cluster[j]=clusterName
            elif len(previousCluster) == 1:
                clusterName=previousCluster[0]
                clusterList[clusterName]=list(set(clusterList[clusterName]+tmpCluster))
                for j in tmpCluster:
                    ID2cluster[j]=clusterName
            else:
                clusterName=previousCluster[0]
                clusterAll=[]
                for j in previousCluster:
                    clusterAll+=clusterList[j]
                for j in range(1,len(previousCluster)):
                    clusterList[previousCluster[j]]=[]
                clusterList[clusterName]=list(set(clusterAll+tmpCluster))
                for j in clusterList[clusterName]:
                    ID2cluster[j]=clusterName
        n+=1
        
    clusterListLen=[len(clusterList[i]) for i in clusterList.keys()]
    cluster_pass=[clusterList[i] for i in clusterList.keys() if len(clusterList[i])>=clusterNumMin]
    cluster_pass_len=[len(i) for i in cluster_pass]
    cluster_pass_idx=[i[0] for i in cluster_pass]
    raw2raw2=open(raw2rawFile)
    line=raw2raw2.readline().split('\t')
    strandDict={}
    while True:
        strandDict[line[0]+'_'+line[5]]=line[4]
        line=raw2raw2.readline()
        if line:
            line=line.split('\t')
            continue
        else:
            break
    raw2raw2.close()
    strandDictTrans={'+':1,'-':0}
    cluster_strand_count=[]
    for i in range(len(cluster_pass)):
        tmpCluster=cluster_pass[i]
        tmpIdx=tmpCluster[0]
        tmpPredict=[]
        for j in range(1,len(tmpCluster)):
            tmpID=tmpCluster[j]
            tmpKey=tmpIdx+'_'+tmpID
            tmpOK=False
            if strandDict.__contains__(tmpKey):
                tmpStrand1=strandDictTrans[strandDict[tmpKey]]
                tmpOK=True
            elif strandDict.__contains__(tmpID+'_'+tmpIdx):
                tmpStrand1=strandDictTrans[strandDict[tmpID+'_'+tmpIdx]]
                tmpOK=True
            else:
                continue
            if tmpOK:
                if os.path.exists(outPrefix+'strandProbability.txt'): # just for easy test of simulation data without strandProbability.txt file
                    if tmpID[0:36] in strandProbability.index:
                        tmpStrand2=strandProbability.loc[tmpID[0:36],'strand']
                        if tmpStrand1==tmpStrand2:
                            tmpPredict.append(1)
                        else:
                            tmpPredict.append(0)
                else:# just for easy test of simulation data without strandProbability.txt file
                    tmpStrand2=1
                    if tmpStrand1==tmpStrand2:
                        tmpPredict.append(1)
                    else:
                        tmpPredict.append(0)
        tmpPredict=np.array(tmpPredict)
        cluster_strand_count.append([sum(tmpPredict==0),sum(tmpPredict==1)])
    #accuraryLH=0.95/0.05
    #diffThreshold=log(len(cluster_pass),accuraryLH)
    cluster2ID={}
    for i in ID2cluster.keys():
        val=ID2cluster[i]
        if cluster2ID.__contains__(val):
            cluster2ID[val].append(i)
        else:
            cluster2ID[val]=[i]
    diffThreshold=0
    cluster_pass_idx_strand=[]
    for i in range(len(cluster_strand_count)):
        if cluster_strand_count[i][0]- cluster_strand_count[i][1]>diffThreshold:
            cluster_pass_idx_strand.append(0)
        elif cluster_strand_count[i][1]- cluster_strand_count[i][0]>diffThreshold:
            cluster_pass_idx_strand.append(1)
        else:
            cluster_pass_idx_strand.append(-1)
    eachReadStrand={}
    for i in range(len(cluster_pass)):
        idxStrand=cluster_pass_idx_strand[i]
        #if idxStrand==-1:
        #    continue
        tmpCluster=cluster_pass[i]
        tmpIdx=tmpCluster[0]
        eachReadStrand[tmpIdx]=idxStrand
        for j in range(1,len(tmpCluster)):
            #print(j)
            tmpID=tmpCluster[j]
            tmpKey=tmpIdx+'_'+tmpID
            tmpOK=False
            if strandDict.__contains__(tmpKey):
                tmpStrand=strandDictTrans[strandDict[tmpKey]]
                tmpOK=True
            elif strandDict.__contains__(tmpID+'_'+tmpIdx):
                tmpStrand=strandDictTrans[strandDict[tmpID+'_'+tmpIdx]]
                tmpOK=True
            else:
                continue
            if tmpOK:
                if idxStrand  == 0:
                    if tmpStrand == 0:
                        eachReadStrand[tmpID]=1
                    else:
                        eachReadStrand[tmpID]=0               
                else:
                    if tmpStrand == 0:
                        eachReadStrand[tmpID]=0
                    else:
                        eachReadStrand[tmpID]=1    
    for i in range(len(cluster_pass)):
        tmpCluster=cluster_pass[i]
        getClusterStrand(tmpCluster)
    clusterDf=pd.DataFrame(eachReadStrand.values(),index=eachReadStrand.keys(),columns=['strand'])
    clusterDf['cluster']=clusterDf.index.map(ID2cluster)
    clusterDf['ID']=clusterDf.index
    #clusterDf[['ID','cluster','strand']].to_csv(outPrefix+'novoCluster.txt',sep='\t',index=None)
    thFile_cluster=thFile_Pass_HQ.loc[clusterDf.index]
    clusterDf=pd.concat([clusterDf,thFile_cluster],axis=1)
    clusterDf=clusterDf[['ID','cluster','strand','readLen','start','end','consLen','copyNum','usage','consensus']]
    clusterDf.to_csv(outPrefix+'novoCluster.txt',sep='\t',index=None)