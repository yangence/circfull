import os, sys, pyfasta, pysam, pandas as pd, numpy as np
from multiprocessing import Pool
errorLen=14 # FSJ is better than BSJ
errorLen2=80

# filter microExon
def keepIndex(x,index,pos):
    new=[]
    x=x.split(',')
    for i in index:
        new.append(x[i])
    return(','.join(new),new[pos])
    
def adjustMicroExon(i):
    x=consFLraw.iloc[i].copy()
    start=[int(i) for i in x['exon_start'].split(',')]
    end=[int(i) for i in x['exon_end'].split(',')]
    keep_index=[]
    for i in range(len(start)):
        if end[i]-start[i]>10:
            keep_index.append(i)
    if len(keep_index)==len(start):
        return(x)
    elif len(keep_index)==0:
        return()
    x['exon_start'],x['start']=keepIndex(x['exon_start'],keep_index,0)
    x['exon_end'],x['end']=keepIndex(x['exon_end'],keep_index,-1)
    x['exon_leftSeq'],x['leftSeq']=keepIndex(x['exon_leftSeq'],keep_index,0)
    x['exon_rightSeq'],x['rightSeq']=keepIndex(x['exon_rightSeq'],keep_index,-1)
    x['circID']=x['chr']+'|'+str(x['start'])+'|'+str(x['end'])
    x['motif']=x['leftSeq']+x['rightSeq']
    x['len']=sum([end[i]-start[i]+1 for i in keep_index])
    x['exonNum']=len(keep_index)
    x['start']=int(x['start'])
    x['end']=int(x['end'])
    return(x)
    

def mapExon(i):
    x=consFL.iloc[i].copy()
    chr=x['chr']
    start=x['exon_start'].split(',')
    end=x['exon_end'].split(',')
    startKeyList=[chr+'_'+i for i in start]
    endKeyList=[chr+'_'+i for i in end]
    num=len(start)
    if num==0:
        return()
    adjS=[]
    adjE=[]
    seqS=[]
    seqE=[]
    for i in range(num):
        startValue=[]
        endValue=[]
        startKey=startKeyList[i]
        endKey=endKeyList[i]
        if ExonSdict.__contains__(startKey):
            startValue=ExonSdict[startKey]
        if ExonEdict.__contains__(endKey):
            endValue=ExonEdict[endKey]
        cValueCan=list(set(startValue) & set(endValue))
        if isSecond and len(cValueCan)>0:
            strandScore=strandFile.loc[x['ID'],'score']
            cValue=[]
            for j in cValueCan:
                if (strandDict[j]=='+' and strandScore>=0) or (strandDict[j]=='-' and strandScore<=0):
                    cValue.append(j)
        else:
            cValue=cValueCan
        if len(cValue)==0:
            adjS.append(int(start[i]))
            adjE.append(int(end[i]))
        elif len(cValue)==1:
            tmp_gtf=gtf_exon.loc[cValue[0]]
            if len(tmp_gtf.shape)==1:
                adjS.append(tmp_gtf['start'])
                adjE.append(tmp_gtf['end'])
            else:
                tmp_gtf_start=np.array(tmp_gtf['start'])
                tmp_gtf_start_diff=abs(tmp_gtf_start-int(start[i]))
                adjS.append(tmp_gtf_start[np.where(tmp_gtf_start_diff==min(tmp_gtf_start_diff))[0][0]])
                tmp_gtf_end=np.array(tmp_gtf['end'])
                tmp_gtf_end_diff=abs(tmp_gtf_end-int(end[i]))
                adjE.append(tmp_gtf_end[np.where(tmp_gtf_end_diff==min(tmp_gtf_end_diff))[0][0]])
        else:
            mSlist=[]
            mElist=[]
            mDifflist=[]
            for j in range(len(cValue)):
                tmp_gtf=gtf_exon.loc[cValue[j]]
                if len(tmp_gtf.shape)==1:
                    mSlist.append(tmp_gtf['start'])
                    mElist.append(tmp_gtf['end'])
                    mDifflist.append(abs(tmp_gtf['start']-int(start[i]))+abs(tmp_gtf['end']-int(end[i])))
                else:
                    tmp_gtf_start=np.array(tmp_gtf['start'])
                    tmp_gtf_start_diff=abs(tmp_gtf_start-int(start[i]))
                    mSlist.append(tmp_gtf_start[np.where(tmp_gtf_start_diff==min(tmp_gtf_start_diff))[0][0]])
                    tmp_gtf_end=np.array(tmp_gtf['end'])
                    tmp_gtf_end_diff=abs(tmp_gtf_end-int(end[i]))
                    mElist.append(tmp_gtf_end[np.where(tmp_gtf_end_diff==min(tmp_gtf_end_diff))[0][0]])
                    mDifflist.append(min(tmp_gtf_start_diff)+min(tmp_gtf_end_diff))
            mDifflist=np.array(mDifflist)
            mIdx=np.where(mDifflist==min(mDifflist))[0][0]
            adjS.append(mSlist[mIdx])
            adjE.append(mElist[mIdx])
    adjS[0]=x['start']
    adjE[-1]=x['end']
    x['len']=sum([adjE[i]-adjS[i]+1 for i in range(num)])
    for i in range(num):
        seqS.append(genome.sequence({'chr': chr, 'start':adjS[i]-2, 'stop':adjS[i]-1}).upper())
        seqE.append(genome.sequence({'chr': chr, 'start':adjE[i]+1, 'stop':adjE[i]+2}).upper())
        adjS[i]=str(adjS[i])
        adjE[i]=str(adjE[i])
    x['exon_start']=','.join(adjS)
    x['exon_end']=','.join(adjE)
    x['exon_leftSeq']=','.join(seqS)
    x['exon_rightSeq']=','.join(seqE)
    return(x)




def adjExon(i):
    x=consFL_fail.iloc[i].copy()
    id=x['circID']
    exonNum=x['exonNum']
    chr=x['chr']
    name=x.name
    adj=False
    if consFL_pass_idxDict.__contains__(name):
        adj=True
        tmp=consFL_pass.loc[name].copy()
        xExonStart=[int(i) for i in x['exon_start'].split(',')]
        xExonEnd=[int(i) for i in x['exon_end'].split(',')]
        if len(tmp.shape)==1:
            tmpExonStart=[int(i) for i in tmp['exon_start'].split(',')]
            tmpExonEnd=[int(i) for i in tmp['exon_end'].split(',')]
            for i in range(exonNum):
                diffStart=abs(xExonStart[i]-tmpExonStart[i])
                diffEnd=abs(xExonEnd[i]-tmpExonEnd[i])
                #if tmpExonStart[i]>=tmpExonEnd[i]:
                #    continue
                if diffStart < errorLen2:
                    xExonStart[i]=tmpExonStart[i]
                if diffEnd < errorLen2:
                    xExonEnd[i]=tmpExonEnd[i]     
        else:
            tmp=tmp.drop_duplicates('key2')
            tmpExonStart=tmp['exon_start'].map(lambda x: [int(i)for i in x.split(',')])
            tmpExonEnd=tmp['exon_end'].map(lambda x: [int(i) for i in x.split(',')])
            for j in range(exonNum):
                minIndex=0
                minLenStart=errorLen2
                minLenEnd=errorLen2
                for i in range(tmpExonStart.shape[0]):  
                    diffStart=abs(xExonStart[j]-tmpExonStart.iloc[i][j])
                    diffEnd=abs(xExonEnd[j]-tmpExonEnd.iloc[i][j])
                    if diffStart < minLenStart:
                        xExonStart[j]=tmpExonStart.iloc[i][j]
                        minLenStart=diffStart
                    if diffEnd < minLenEnd:
                        xExonEnd[j]=tmpExonEnd.iloc[i][j]
                        minLenEnd=diffEnd
        xSeqS=[]
        xSeqE=[]
        for i in range(exonNum):
            xSeqS.append(genome.sequence({'chr': chr, 'start':xExonStart[i]-2, 'stop':xExonStart[i]-1}).upper())
            xSeqE.append(genome.sequence({'chr': chr, 'start':xExonEnd[i]+1, 'stop':xExonEnd[i]+2}).upper())
        x['len']=sum([xExonEnd[i]-xExonStart[i]+1 for i in range(len(xExonStart))])
        xExonStart=[str(i) for i in xExonStart]
        xExonEnd=[str(i) for i in xExonEnd]
        x['exon_start']=','.join(xExonStart)
        x['exon_end']=','.join(xExonEnd)
        x['exon_leftSeq']=','.join(xSeqS)
        x['exon_rightSeq']=','.join(xSeqE)
    return(x)


def adjFL(options):
    global genome,ExonSdict,ExonEdict,strandDict,consFL,consFLraw,isSecond,gtf_exon,consFL_fail,consFL_pass,consFL_pass_idxDict,strandFile
    genomeFile=options[0]
    gtfFile=options[1]
    outPrefix=options[2]
    thread=int(options[3])
    genome = pyfasta.Fasta(genomeFile)
    consFLraw=pd.read_csv(outPrefix+'constructFL_Normal.txt',sep='\t',dtype={'exon_start':str,'exon_end':str})
    isSecond=False
    if len(options)>4:
        strandFile=options[4]
        isSecond=True
        strandFile.index=strandFile['ID']
    pool=Pool(processes=thread)
    consFL_list=pool.map(adjustMicroExon,[i for i in range(consFLraw.shape[0])])
    pool.close()
    pool.join()
    consFL_list=[i for i in consFL_list if len(i)>0]
    consFL=pd.DataFrame(consFL_list)

    # gtf2dict
    gtf = pd.read_csv(gtfFile,sep="\t",names=['chr','source','type','start','end','score','strand','phase','attributes'],comment='#')
    gtf_exon=gtf[gtf['type']=='exon'].sort_values(by=['chr','start'])
    gtf_exon['key']=gtf_exon['chr'].values+'|'+gtf_exon['start'].map(str).values+gtf_exon['strand'].values+gtf_exon['end'].map(str).values
    gtf_exon['gene']=gtf_exon.apply(lambda x: x['attributes'].split(';')[0].split(' ')[1][1:-1],axis=1)
    gtf_exon['geneKey']=gtf_exon['gene']+'|'+gtf_exon['key']
    gtf_exon=gtf_exon.set_index('gene')
    gtf_exon=gtf_exon.sort_index()
    if isSecond:
        path=outPrefix+'../'
    else:
        path=outPrefix
    if os.path.exists(os.path.dirname(gtfFile)+'/ExonSdict_fsj.npy'):
        ExonSdict=np.load(os.path.dirname(gtfFile)+'/ExonSdict_fsj.npy',allow_pickle=True).item()
        ExonEdict=np.load(os.path.dirname(gtfFile)+'/ExonEdict_fsj.npy',allow_pickle=True).item()
        strandDict=np.load(os.path.dirname(gtfFile)+'/strandDict_fsj.npy',allow_pickle=True).item()
    elif os.path.exists(path+'ExonSdict_fsj.npy'):
        ExonSdict=np.load(path+'ExonSdict_fsj.npy',allow_pickle=True).item()
        ExonEdict=np.load(path+'ExonEdict_fsj.npy',allow_pickle=True).item()
        strandDict=np.load(path+'strandDict_fsj.npy',allow_pickle=True).item()
    else:
        ExonSdict={}
        ExonEdict={}
        strandDict={}
        for i in range(gtf_exon.shape[0]):
            tmp=gtf_exon.iloc[i]
            tmpChr=tmp['chr']
            tmpStart=tmp['start']
            tmpEnd=tmp['end']
            tmpValue=tmp.name
            strandDict[tmpValue]=tmp['strand']
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
        np.save(outPrefix+'ExonSdict_fsj.npy',ExonSdict)
        np.save(outPrefix+'ExonEdict_fsj.npy',ExonEdict)
        np.save(outPrefix+'strandDict_fsj.npy',strandDict)
        
    pool=Pool(processes=thread)
    consFL_list=pool.map(mapExon,[i for i in range(consFL.shape[0])])
    pool.close()
    pool.join()
    consFL_list=[i for i in consFL_list if len(i)>0]
    consFL=pd.DataFrame(consFL_list)
    # common adjust rare internal structure
    consFL['key1']=consFL['circID']+'|'+consFL['exonNum'].map(str)+'|'+consFL['len'].map(str)
    consFL['key2']=consFL['circID']+'|'+consFL['exonNum'].map(str)
    consFL.index=consFL['key2']
    consFL_key2_uniq=list(set(consFL['key2']))
    consFL=consFL.sort_index()
    consFL_id_Pass=[]
    for i in range(len(consFL_key2_uniq)):
        tmpFL=consFL.loc[consFL_key2_uniq[i]]
        if len(tmpFL.shape)>1:
            tmpFL_count=tmpFL['key1'].value_counts()
            if len(tmpFL_count)>1:
                consFL_id_Pass.append(tmpFL_count.index[0])
    consFL_pass=consFL.loc[consFL['key1'].isin(consFL_id_Pass)].copy()
    consFL_fail=consFL.loc[consFL['key1'].isin(set(consFL['key1'])-set(consFL_id_Pass))].copy()
    consFL_pass_idxDict={}
    for i in consFL_pass.index:
        consFL_pass_idxDict[i]=1 

    
    pool=Pool(processes=thread)
    resultDf=pool.map(adjExon,[i for i in range(consFL_fail.shape[0])])
    pool.close()
    pool.join()
    resultDf=pd.DataFrame(resultDf)
    newDf=pd.concat([resultDf,consFL_pass],axis=0)
    newDf=newDf[['ID','circID','chr','start','end','len','exonNum','exon_start','exon_end','motif','leftSeq','rightSeq','exon_leftSeq','exon_rightSeq']]
    newDf=newDf.sort_values(['circID','len'])
    newDf.to_csv(outPrefix+'constructFL_Normal_adj.txt',sep='\t',index=None)
    