import pysam,pandas as pd,sys,os,numpy as np,pyfasta,mappy as mp
from interval import Interval
from progressbar import *
from multiprocessing import Pool
from .RG_circFL_output import getSize,FL2bed
from .RG_detectBS import getReadInfo
filterTh=0.1
BSfilterTh=0.5
FSfilterTh=0.5
unmapped_readLen_Th=50
start_end_diff_Th1=10
start_end_diff_Th2=5
def getRef(isoID):
    chr,start,end=isoID.split('|')
    start=[int(i) for i in start.split(',')]
    end=[int(i) for i in end.split(',')]
    seq=''
    for i in range(len(start)):
        seq+=genome.sequence({'chr':chr,'start':start[i],'stop':end[i]}).upper()
    return(seq)
def ref_pos2Type(eachMap):
    typeList=[]
    n=0
    for cigar in eachMap.cigar:
        if cigar[1]==0:
            typeList.extend([0 for i in range(cigar[0])])
            n+=cigar[0]
        elif cigar[1]==2:
            typeList.extend([2 for i in range(cigar[0])])
            n+=cigar[0]
        elif cigar[1]==3:
            typeList.extend([3 for i in range(cigar[0])])
            n+=cigar[0]
    
    return(typeList)

def findReadPos(eachMap,targetPos):
    ref_len=0
    read_len=0
    readType=[]
    for cigar in eachMap.cigar:
        if cigar[1]==0:
            read_len+=cigar[0]
            ref_len+=cigar[0]
        elif cigar[1]==1:
            read_len+=cigar[0]
        elif cigar[1]==2:
            ref_len+=cigar[0]
        elif cigar[1]==3:
            ref_len+=cigar[0]
        elif cigar[1]==4:
            read_len+=cigar[0]
        elif cigar[1]==5:
            read_len+=cigar[0]
        else:
            continue
        if ref_len in targetPos:
            readType.append(cigar[1])
    return(readType)

def judgeRef(eachMap,targetPos,errorLen=4):
    pRef=ref_pos2Type(eachMap)
    for i in targetPos:
        pRef_bs=pRef[(i-errorLen):(i+errorLen)]
        if sum(pRef_bs)==0:
            return(True)
    return(False)

def judgeRead(eachMap,targetPos,errorLen=4):
    for i in targetPos:
        tmpPos= range(i-errorLen,i+errorLen)
        readType=findReadPos(eachMap,tmpPos)
        if sum(readType)==0:
            return(True)
    return(False)
def judgeJun(eachMap,targetPos):
    if judgeRef(eachMap,targetPos) and judgeRead(eachMap,targetPos):
        return(True)
    else:
        return(False)
def evalCirc(i):
    tmp=FL_noDup.iloc[i,:]
    refSeq=getRef(tmp['isoID'])
    refLen=len(refSeq)
    readID=tmp['readID'].split(',')
    countBS=0
    count=0
    exon_start=[int(i) for i in tmp['exon_start'].split(',')]
    exon_start=[i-exon_start[0] for i in exon_start]
    
    count_SE=0 # the potential chimeric reads
    for j in range(len(readID)):
        eachRead=readID[j]
        copyN=readLen[eachRead]//refLen +2
        eachRef=refSeq * int(copyN)
        eachAlign = mp.Aligner(seq=eachRef)
        readSeq=readDict[eachRead]
        isBSPass=False
        for eachMap in eachAlign.map(readSeq):
            if eachMap.is_primary:
                #print(j)
                #print(eachMap.r_st)
                #print(eachMap.r_en)
                #print(len(eachRef))
                targetPos=[i*refLen -eachMap.r_st for i in range(copyN) if (i*refLen>=eachMap.r_st) and (i*refLen<=eachMap.r_en)]
                #print(targetPos)
                if len(targetPos)>0:
                    start_diff=min(targetPos[0],abs(refLen-targetPos[0]))
                    end_diff=eachMap.r_en % refLen
                    end_diff=min(end_diff,abs(refLen-end_diff))
                    SE_diff=min(start_diff,end_diff)
                    read_diff=readLen[eachRead]-(eachMap.r_en-eachMap.r_st)
                    if (SE_diff<start_end_diff_Th1) & (read_diff>unmapped_readLen_Th):
                        count_SE+=-1
                    else:
                        count_SE+=1
                    '''
                    if (SE_diff<start_end_diff_Th2):
                        count_SE+=-1
                    else:
                        if  (SE_diff<start_end_diff_Th1) & (read_diff>unmapped_readLen_Th):
                            count_SE+=-1
                        else:
                            count_SE+=1
                    '''
                
                isPass=True
                for q in exon_start:
                    newTargetPos=[m+q for m in targetPos]
                    if judgeJun(eachMap,newTargetPos):
                        if q==0 and isBSPass!=True:
                            countBS+=1
                            isBSPass=True
                        continue
                    else:
                        isPass=False
                        break
                if isPass==True:
                    count+=1
                    break                
    return([countBS,count,count_SE])
def calReadSplice(i):
    tmp=FL_noDup.iloc[i,:]
    readArr=[]
    targetI=Interval(tmp['start']-1,tmp['end'])
    tmpRegionS=tmp['chr']+':'+str(tmp['start'])+'-'+str(tmp['start'])
    tmpRegionE=tmp['chr']+':'+str(tmp['end'])+'-'+str(tmp['end'])
    tmpFileS=RG_tmp_outPrefix+str(i)+'_FO_s.sam'
    tmpFileE=RG_tmp_outPrefix+str(i)+'_FO_e.sam'
    CMDS='samtools view -h %s %s>%s' % (bamFile,tmpRegionS,tmpFileS)
    CMDE='samtools view -h %s %s>%s' % (bamFile,tmpRegionE,tmpFileE)
    tmpSAM_s=''
    tmpSAM_e=''
    isPass=True
    if not os.system(CMDS):
        tmpSAM_s=pysam.AlignmentFile(tmpFileS)
    else:
        isPass=False
    if not os.system(CMDE):
        tmpSAM_e=pysam.AlignmentFile(tmpFileE)
    else:
        isPass=False
    if not isPass:
        return([1,1])
    readDict_tmp={}
    for read in tmpSAM_s.fetch():
        readDict_tmp[read.qname]=1
        if read.flag & 260 >0:
            continue
        if dict_ID2type.__contains__(read.qname):
            continue
        readInfo=getReadInfo(read)
        for i in range(len(readInfo[0])):
            readI=Interval(readInfo[0][i],readInfo[1][i])
            if readI.overlaps(targetI):
                readArr.append(readInfo)
    for read in tmpSAM_e.fetch():
        if readDict_tmp.__contains__(read.qname):
            continue
        if read.flag & 260 >0:
            continue
        if dict_ID2type.__contains__(read.qname):
            continue
        readInfo=getReadInfo(read)
        for i in range(len(readInfo[0])):
            readI=Interval(readInfo[0][i],readInfo[1][i])
            if readI.overlaps(targetI):
                readArr.append(readInfo)
    count_s=0
    count_i=0
    for readInfo in readArr:
        tmpS=0
        tmpI=0
        if len(readInfo[0])>1:
            for j in range(len(readInfo[0])-1):
                readI=Interval(readInfo[1][j],readInfo[0][j+1])
                if readI.overlaps(targetI):
                    tmpS=1
                    #print(readInfo)
                    break
        else:
            tmpI=1
        count_s+=tmpS
        count_i+=tmpI
    os.system('rm %s' % tmpFileS)
    os.system('rm %s' % tmpFileE)
    return([count_s,count_i])

def filterOut(readPrefix,outPrefix,fastqFile,genomeFile,thread=1,rmskFile=False):
    global bam,FL_noDup,dict_ID2type,RG_tmp_outPrefix,bamFile,genome,readLen,readDict
    fastq=open(fastqFile)
    
    genome=pyfasta.Fasta(genomeFile)
    bamFile=readPrefix+'test.minimap2.bam'
    RG_tmp_outPrefix=outPrefix+'tmp/'
    FL_ID2type=pd.read_csv(readPrefix+'explainFL_ID2Type.txt',sep='\t')
    dict_ID2type=dict(zip(FL_ID2type['ID'],FL_ID2type['type']))
    bam=pysam.AlignmentFile(bamFile,'r')
    FL_noDup=pd.read_csv(outPrefix+'circFL_Normal.txt',sep='\t')
    totalCount=FL_noDup.loc[:,['circID','readCount']].groupby('circID').agg('sum')
    readDict={}
    readLen={}
    line1=fastq.readline()
    while True:
        if line1:
            line2=fastq.readline().strip()
            line3=fastq.readline()
            line4=fastq.readline()
            ID=line1.strip().split()[0][1:]
            readDict[ID]=line2
            readLen[ID]=len(line2)
            line1=fastq.readline()
        else:
            break
    FL_fo=FL_noDup.loc[:,['circID','isoID','readCount']].copy()
    FL_fo['totalCount']=totalCount.loc[FL_noDup.circID,'readCount'].tolist()
    readSplice=[]
    perfectMatch=[]
    pool=Pool(processes=thread)
    widgets = [ 'CalSpliceRatio: ' , Percentage () , ' ' , Bar ( marker = RotatingMarker ( ) ) ,' ' , ETA ( ) , ' '  ]
    bar=progressbar.ProgressBar(widgets=widgets,maxval=FL_noDup.shape[0]).start()
    eachBin=20*thread
    for i in range(0,FL_noDup.shape[0],eachBin):
        tmpRS=pool.map(calReadSplice,range(i,min(i+eachBin,FL_noDup.shape[0])))
        readSplice.extend(tmpRS)
        tmpBM=pool.map(evalCirc,range(i,min(i+eachBin,FL_noDup.shape[0])))
        perfectMatch.extend(tmpBM)
        bar.update(i)
    pool.close()
    pool.join()
    
    
    FL_fo['FSJ']= [i[0] for i in readSplice] # linear junction
    FL_fo['FN']= [i[1] for i in readSplice] # non junction
    FL_fo['spliceRatio']=(FL_fo['FSJ']+FL_fo['totalCount'])/(FL_fo['FSJ']+FL_fo['totalCount']+FL_fo['FN'])
    FL_fo['perfectM_BS']=[i[0] for i in perfectMatch]
    FL_fo['perfectM_FS']=[i[1] for i in perfectMatch]
    FL_fo['perfectM_SEdiff']=[i[2] for i in perfectMatch]

    totalM=FL_fo.loc[:,['circID','perfectM_BS']].groupby('circID').agg('sum')
    FL_fo['totalM_BS']=totalM.loc[FL_fo.circID,'perfectM_BS'].tolist()
    FL_fo['matchBSratio']=FL_fo['totalM_BS']/FL_fo['totalCount']
    FL_fo['matchFSratio']=FL_fo['perfectM_FS']/FL_fo['readCount']
    

    if not rmskFile:
        FL_fo.to_csv(outPrefix+'circFL_Normal_filterOut.detail',sep='\t',header=True,index=None)
        FL_noDup=FL_noDup.loc[(FL_fo['spliceRatio']>filterTh) & (FL_fo['matchBSratio']>BSfilterTh)& (FL_fo['matchFSratio']>FSfilterTh)& (FL_fo['perfectM_SEdiff']>0),:]
        FL_noDup.to_csv(outPrefix+'circFL_Normal_pass.txt',sep='\t',header=True,index=None)
        FL_bed=FL2bed(FL_noDup)
        FL_bed.to_csv(outPrefix+'circFL_Normal_pass.bed',sep='\t',header=None,index=None)
    else:
        start_bed=FL_noDup.loc[:,['chr','start','start','isoID']].copy()
        start_bed.columns=['chr','start','end','isoID']
        end_bed=FL_noDup.loc[:,['chr','end','end','isoID']].copy()
        end_bed.columns=['chr','start','end','isoID']
        start_bed.iloc[:,1]=start_bed.iloc[:,1]-30
        start_bed.iloc[:,2]=start_bed.iloc[:,2]+30
        end_bed.iloc[:,1]=end_bed.iloc[:,1]-30
        end_bed.iloc[:,2]=end_bed.iloc[:,2]+30
        start_bed.to_csv(outPrefix+'tmp_s.bed',sep='\t',header=None,index=None)
        end_bed.to_csv(outPrefix+'tmp_e.bed',sep='\t',header=None,index=None)
        rmsk_s=outPrefix+'rmsk_s.bed'
        rmsk_e=outPrefix+'rmsk_e.bed'
        os.system('bedtools intersect -wo -a %s -b %s >%s 2>/dev/null' % (outPrefix+'tmp_s.bed',rmskFile, rmsk_s))
        os.system('bedtools intersect -wo -a %s -b %s >%s 2>/dev/null' % (outPrefix+'tmp_e.bed',rmskFile, rmsk_e))
        s_repeat=[]
        with open(rmsk_s) as fin:
            line=fin.readline()
            while line:
                arr=line.split('\t')
                s_repeat.append(arr[3])
                line=fin.readline()
        e_repeat=[]
        with open(rmsk_e) as fin:
            line=fin.readline()
            while line:
                arr=line.split('\t')
                e_repeat.append(arr[3])
                line=fin.readline()

        boudary_repeat=set(s_repeat) & set(e_repeat)
        FL_fo['boudary_repeat']=0
        FL_fo.loc[FL_fo.isoID.isin(boudary_repeat),'boudary_repeat']=1
        FL_fo.to_csv(outPrefix+'circFL_Normal_filterOut.detail',sep='\t',header=True,index=None)
        FL_noDup=FL_noDup.loc[(FL_fo['spliceRatio']>filterTh) & (FL_fo['matchBSratio']>BSfilterTh)& (FL_fo['matchFSratio']>FSfilterTh)& (FL_fo['perfectM_SEdiff']>0) & (FL_fo['boudary_repeat']==0),:]
        FL_noDup.to_csv(outPrefix+'circFL_Normal_pass.txt',sep='\t',header=True,index=None)
        FL_bed=FL2bed(FL_noDup)
        FL_bed.to_csv(outPrefix+'circFL_Normal_pass.bed',sep='\t',header=None,index=None)
        os.remove(rmsk_s)
        os.remove(rmsk_e)
        os.remove(outPrefix+'tmp_s.bed')
        os.remove(outPrefix+'tmp_e.bed')