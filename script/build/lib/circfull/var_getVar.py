import pandas as pd,pyfasta,os,pysam,time,numpy as np,sys
from multiprocessing import Pool
def getVCFLine(vcf):
    fin=open(vcf)
    arrAD={}
    while True:
        eachLine=fin.readline()
        if not eachLine:
            return(arrAD)
        if eachLine[0]=='#':
            continue
        each=eachLine.strip().split('\t')
        refBase=each[3]
        altBase=each[4].split(',')[0]
        pos=int(each[1])
        x=each[9]
        eachAD=[int(i) for i in x.split(':')[1].split(',')]
        if len(eachAD)>3:
            currentBase=refBase
        elif altBase not in ['A','C','G','T','<*>']:
            currentBase=refBase
        elif sum(eachAD)==0:
            currentBase='N'
        elif eachAD[1]>eachAD[0]:
            currentBase=altBase
        else:
            currentBase=refBase
        arrAD[pos]=currentBase
    return(arrAD)

def getEachVCF(sam,ref):
    fin=open(sam)
    headLine=[]
    mapLine={}
    while True:
        eachLine=fin.readline()
        if not eachLine:
            break
        if eachLine[0]=='@':
            headLine.append(eachLine)
        else:
            readID=eachLine.split('\t')[0]
            if mapLine.__contains__(readID):
                mapLine[readID].append(eachLine)
            else:
                mapLine[readID]=[eachLine]
    fin.close()
    headLine=''.join(headLine)
    n=0
    fin=open(ref)
    refLen=len(fin.readlines()[1].strip())
    fin.close()
    #print(refLen)
    baseDict={}
    baseDict['A']={}
    baseDict['C']={}
    baseDict['G']={}
    baseDict['T']={}
    baseDict['N']={}
    for i in ['A','C','G','T','N']:
        for j in range(refLen):
            baseDict[i][j+1]=0
    #print(len(mapLine.keys()))
    tmpSam=sam+'.'+str(0)
    for i in mapLine.keys():
        #print(i)
        eachMap=headLine+''.join(mapLine[i])
        fout=open(tmpSam,'w')
        fout.write(eachMap)
        fout.close()
        #t1=time.time()
        os.system('/home/lzl/miniconda3/bin/samtools mpileup -v -u -B -I -Q 1 -t AD -f '+ref +' '+tmpSam+' >'+tmpSam+'.vcf 2>/dev/null')
        #t2=time.time()
        #print('t1 time:%f' % (t2-t1))
        arrAD=getVCFLine(tmpSam+'.vcf')
        #t3=time.time()
        #print('t2 time:%f' % (t3-t2))
        for j in arrAD.keys():
            baseDict[arrAD[j]][j]+=1
        #t4=time.time()
        #print('t3 time:%f' % (t4-t3))
    os.remove(tmpSam)
    os.remove(tmpSam+'.vcf')
    return(baseDict)

def AD2ED(baseDict,ref):
    ad=pd.DataFrame(baseDict)
    fin=open(ref)
    refSeq=fin.readlines()[1].strip()
    fin.close()
    edDict={}
    edDict['ref']={}
    edDict['alt']={}
    edDict['ref_count']={}
    edDict['alt_count']={}
    for i in edDict.keys():
        for j in range(len(refSeq)):
            edDict[i][j+1]={}
    for i in ad.index:
        refBase=refSeq[i-1]
        edDict['ref'][i]=refBase
        eachLine=ad.loc[i,:][0:4]
        if refSeq[i-1] not in ['A','C','G','T']:
            edDict['alt'][i]='*'
            edDict['ref_count'][i]=0
            edDict['alt_count'][i]='0'
        else:
            edDict['ref_count'][i]=eachLine[refBase]
            if sum(eachLine)==eachLine[refBase]:
                edDict['alt'][i]='*'
                edDict['alt_count'][i]='0'
            else:
                otherCount={}
                for j in ['A','C','G','T']:
                    if j !=refBase:
                        if eachLine[j]!=0:
                            otherCount[j]=eachLine[j]
                otherCountSort=sorted(zip(otherCount.values(), otherCount.keys()))[::-1]
                edDict['alt'][i]=','.join([i[1] for i in otherCountSort])
                edDict['alt_count'][i]=','.join([str(i[0]) for i in otherCountSort])
    return(pd.DataFrame(edDict))

def getFreq(refCount,altCount):
    eachAD=[refCount]
    eachAD.extend([int(i) for i in altCount.split(',')])
    if len(eachAD)<2 or sum(eachAD)<5:
        return(False)
    #if eachAD[0]<2:
    #    return(False)
    if len(eachAD)>2:
        if eachAD[2]*2>eachAD[1]:
            return(False)
    if eachAD[1]/sum(eachAD)>0.05:# and eachAD[1]>2
        if len(eachAD)>2:
            if eachAD[2]/sum(eachAD)>0.05:
                return(False)
        return(True)
    else:
        return(False)

def formED(df,chrNum,left,seqLen,circKey):
    edPos=list(df.index)
    for i in edPos:
        if (i<=seqLen) and ((i+seqLen) in edPos):
            if df.loc[i,'ref_count']<df.loc[i+seqLen,'ref_count']:
                df.loc[i,:]=df.loc[i+seqLen,:]
    newDF=df.loc[[i<=seqLen for i in edPos],:].copy()
    newDF.loc[:,'alt']=newDF.loc[:,'alt'].apply(lambda x: x.split(',')[0])
    newDF.loc[:,'alt_count']=newDF.loc[:,'alt_count'].apply(lambda x: int(x.split(',')[0]))
    newDF.loc[:,'chr']=[chrNum for i in range(newDF.shape[0])]
    newDF.loc[:,'pos']=[i+left-1 for i in list(newDF.index)]
    newDF.loc[:,'ratio']=newDF.loc[:,'alt_count']/(newDF.loc[:,'alt_count']+newDF.loc[:,'ref_count'])
    newDF.loc[:,'key']=[circKey for i in range(newDF.shape[0])]
    return(newDF)

def getEd(circKey):
    IDlist=list(FLdf.loc[circKey,'ID'].values)
    fastqOutput=''
    for i in IDlist:
        [r1,r2,r3,r4]=fastqList[i]
        currentID=r1.split(' ')[0][1:]
        #tmpLen=int(len(r2)/2)
        #r2=r2[:tmpLen]+'A'+r2[(tmpLen+1):]
        fastqOutput+="%s%s%s%s\n" % (r1,r2,r3,r4)
        
    eachFaPos=faPos[currentID].split(";")
    eachChr=eachFaPos[0]
    eachLeft=int(eachFaPos[1])
    eachRight=int(eachFaPos[2])
    eachLen=int(eachFaPos[3])
    eachFaSeq=faOutput[currentID]
    foFastq=open(outPrefixTmp+"seq1_"+currentID+".fastq",'w')
    foFastq.write(fastqOutput)
    foFastq.close()
    foFa=open(outPrefixTmp+"seq2_"+currentID+".fa",'w')
    foFa.write(eachFaSeq)
    foFa.close()
    tmpSAM=outPrefixTmp+currentID+".sam"
    cmd1="minimap2 -ax splice -N 1 -k14 "+outPrefixTmp+"seq2_"+currentID+".fa "+outPrefixTmp+"seq1_"+currentID+".fastq>" +tmpSAM+" 2>/dev/null"
    os.system(cmd1)
    ref=outPrefixTmp+"seq2_"+currentID+".fa"
    tmpAD=getEachVCF(tmpSAM,ref)
    tmpED=AD2ED(tmpAD.copy(),ref)
    passED=tmpED.loc[tmpED.apply(lambda x:getFreq(x[2],x[3]),axis=1),:].copy()
    edDF=formED(passED.copy(),eachChr,eachLeft,eachLen,circKey)
    if os.path.exists(outPrefixTmp+"seq1_"+currentID+".fastq"):
        os.remove(outPrefixTmp+"seq1_"+currentID+".fastq")
    if os.path.exists(outPrefixTmp+"seq2_"+currentID+".fa"):    
        os.remove(outPrefixTmp+"seq2_"+currentID+".fa")
    if os.path.exists(outPrefixTmp+currentID+".sam"):    
        os.remove(outPrefixTmp+currentID+".sam")
    if os.path.exists(outPrefixTmp+"seq2_"+currentID+".fa.fai"):
        os.remove(outPrefixTmp+"seq2_"+currentID+".fa.fai")
    if edDF.shape[0]==0:
        return([])
    edList=np.array(edDF.copy())
    edList=[i for i in edList if len(i)>0]
    return(np.ravel(edList).tolist())




def getVar(genomeFile,RG_out,outPrefix,fastqFile,thread):
    global outPrefixTmp,FLdf,fastqList,faPos,faOutput
    outPrefixTmp=outPrefix+'tmp/'
    genome = pyfasta.Fasta(genomeFile)


    FLdf=pd.read_csv(RG_out+'result_Normal.txt',sep='\t')


    FLdf['key']=FLdf['chr']+'|'+FLdf['exon_start']+'|'+FLdf['exon_end']
    edThreshold=4
    pass_key=list(FLdf['key'].value_counts()[FLdf['key'].value_counts()>edThreshold].index)
    FLdf.index=FLdf['key']

    fastqFile=open(fastqFile)
    targetID=list(FLdf['ID'])
    targetKey=dict(zip(targetID,[1 for i in range(len(targetID))]))
    fastqList={}
    while True:
        r1=fastqFile.readline()
        if r1:  
            r2=fastqFile.readline()
            r3=fastqFile.readline()
            r4=fastqFile.readline().strip('\n')
            currentID=r1.split(' ')[0][1:]
            if targetKey.__contains__(currentID):
                fastqList[currentID]=[r1,r2,r3,r4]
        else:
            break
    fastqFile.close()

    hangLen=0
    faOutput={}
    faPos={}
    for i in range(FLdf.shape[0]):
        each=FLdf.iloc[i]
        chr=each['chr']
        exonS=[int(j) for j in each['exon_start'].split(',')]
        exonE=[int(j) for j in each['exon_end'].split(',')]
        maxFa=genome.sequence({'chr': chr, 'start':exonS[0]+1-hangLen, 'stop':exonE[-1]+hangLen})
        faBS=''
        for j in range(2):
            faBS+=maxFa
        name=">ref"
        faOutput[each['ID']]="%s\n%s" % (name, faBS)
        faPos[each['ID']]=chr+";"+str(exonS[0]+1-hangLen)+";"+str(exonE[-1])+";"+str(exonE[-1]-exonS[0]+2*hangLen)

    pool=Pool(processes=thread)
    poolResult=[i for i in pool.map(getEd,pass_key) if len(i)>0]
    #print(poolResult)
    #print('result')
    result=[]
    for i in poolResult:
        result.extend(i)
    result=np.array(result)
    #print(result)
    pool.close()
    pool.join()
    result_res=result.reshape([int(len(result)/8),8])
    resultDF=pd.DataFrame(result_res)
    resultDF.columns=['ref','alt','refCount','altCount','chr','pos','ratio','isoID']
    resultDF=resultDF[['isoID','chr','pos','ref','alt','refCount','altCount','ratio']]
    resultDF.to_csv(outPrefix+'circ_Var.txt',sep='\t',index=False)

if __name__=='__main__':
    main()