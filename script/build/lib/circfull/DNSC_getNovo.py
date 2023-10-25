import pandas as pd, numpy as np,sys,os
from multiprocessing import Pool
from progressbar import *

def getTrf(i):
    tmp=thFile_Pass_HQ.iloc[i].copy()
    isTrf=False
    currentID=tmp['ID']
    r2=tmp['consensus']
    consLen=tmp['consLen']
    readOut=open(outPrefix+'tmp/'+currentID+'.fa','w')
    readOut.write('>'+currentID+'\n'+r2)
    readOut.close()
    cmd='trf '+outPrefix+'tmp/'+currentID+'.fa  2 5 7 80 5 5 2000  -h -ngs >'+outPrefix+'tmp/'+currentID+'.trf'
    os.system(cmd)
    thIn=open(outPrefix+'tmp/'+currentID+'.trf')
    eachTrf=thIn.readline()
    if eachTrf:
        eachTrf=thIn.readline().split(' ')
        if ((int(eachTrf[1])-int(eachTrf[0]))>0.5*consLen) or (int(eachTrf[7])>40):
            isTrf=True
    thIn.close()
    os.remove(outPrefix+'tmp/'+currentID+'.fa')
    os.remove(outPrefix+'tmp/'+currentID+'.trf')
    return([currentID,isTrf])

def getNovo(outDir,fastqFile,thread,isCluster):
    global outPrefix
    global thFile_Pass_HQ
    outPrefix=outDir
    thFile=pd.read_csv(outPrefix+'TideHunter.tab',sep='\t',names=['readName','consN','readLen','start','end','consLen','copyNum','fullLen','consensus'])
    if thFile['readName'][0][0]=='@' or thFile['readName'][0][0]=='>':
        thFile['ID']=[i[1:] for i in thFile['readName']]
    else:
        thFile['ID']=[i for i in thFile['readName']]
    thFile=thFile.set_index('readName').sort_index()
    thFile['usage']=(thFile['end']-thFile['start']+1)/thFile['readLen']
    thFile_Pass_HQ=thFile.loc[(thFile.usage>0.8)&(thFile.copyNum>=2) &(thFile.consLen>=20)]

    pool=Pool(processes=thread)
    widgets = [ 'getTrf: ' , Percentage () , ' ' , Bar ( marker = RotatingMarker ( ) ) ,' ' , ETA ( ) , ' '  ]
    bar=progressbar.ProgressBar(widgets=widgets,maxval=thFile_Pass_HQ.shape[0]).start()
    eachBin=20*thread
    allResult_list=[]
    for i in range(0,thFile_Pass_HQ.shape[0],eachBin):
        allResult_list.extend(pool.map(getTrf,range(i,min(i+eachBin,thFile_Pass_HQ.shape[0]))))
        bar.update(i)
    result=allResult_list
    pool.close()
    pool.join()

    trf_result=pd.DataFrame(result,columns=['ID','trf'])
    trf_result=trf_result.loc[[not i for i in trf_result.trf]]
    thFile_Pass_HQ=thFile_Pass_HQ.set_index('ID')
    thFile_new=thFile_Pass_HQ.loc[trf_result.ID].copy()
    thFile_new.to_csv(outPrefix+'TideHunter_Pass.tab',sep='\t')

    if isCluster:
        readOut=open(outPrefix+'rawseq.fa','w')
        targetID=list(thFile_new.index)
        targetKey=dict(zip(targetID,[1 for i in range(len(targetID))]))
        fastqIn=open(fastqFile)
        rawList=[]
        while True:
            r1=fastqIn.readline()
            if r1:  
                r2=fastqIn.readline()
                r3=fastqIn.readline()
                r4=fastqIn.readline()
                currentID=r1.split(' ')[0][1:]
                if targetKey.__contains__(currentID):
                    readOut.write('>'+r1[1:]+r2)
            else:
                break
        fastqIn.close()
        readOut.close()
        seqFout=open(outPrefix+'novoseq.fa','w')
        for i in range(thFile_new.shape[0]):
            seqFout.write('>'+thFile_new.iloc[i].name+'\n'+thFile_new.iloc[i]['consensus']+'\n')