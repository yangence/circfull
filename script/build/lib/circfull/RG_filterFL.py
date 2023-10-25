import numpy as np,pandas as pd,sys,pyfasta,os,pysam
from interval import Interval # please 'pip uninstall pyinterval' and install this package by 'pip install interval'
from multiprocessing import Pool
from progressbar import *
import mappy as mp


# F1: fusion in different chromosome # F2: fusion in same chromosome
# FC1: fusion candidate in different chromosome # FC2: fusion candiate in same chromosome
# C1: chimeric in different chromosome # C2: chimeric in same chromsome
# N: normal # NC: normal candidate # M: multiple chromosome  # U: Unknown or unsure type

def getSeq(each):
    exonStart=[int(i) for i in each['exon_start'].split(',')]
    exonEnd=[int(i) for i in each['exon_end'].split(',')]
    chr=each['chr']
    seq=''
    for i in range(len(exonStart)):
        seq+=genome.sequence({'chr': chr, 'start':exonStart[i], 'stop':exonEnd[i]})
    return(seq)

def calRead(read_list,refLen):
    align_start=[]
    align_end=[]
    for read in read_list:
        qS=read.query_alignment_start
        qE=read.query_alignment_end
        readcigar=read.cigar
        if readcigar[0][0] == 5:
            qS+=readcigar[0][1]
            qE+=readcigar[0][1]
        align_start.append(qS)
        align_end.append(qE)
    query_pos=sorted(zip(align_start,align_end))
    if len(query_pos)==1:
        return(0)
    threshold_long=0.2 * refLen
    threshold_short=0.2 * refLen
    count_long=0
    count_short=0
    for i in range(1,len(query_pos)):
        #if (query_pos[i][0]-query_pos[i-1][1]) < -threshold_long:  # judge if ref longer than real
            #count_long+=1
        if (query_pos[i][0]-query_pos[i-1][1]) > threshold_short: # judge if ref shorter than real
            count_short+=1
    if (count_short/(len(query_pos)-1)>0.5):
        return(0)
    return(1)
def mp_calRead(mapAll,refLen):
    align_start=[]
    align_end=[]
    for read in mapAll:
        qS=read.q_st
        qE=read.q_en
        readcigar=read.cigar
        if readcigar[0][1] == 5:
            qS+=readcigar[0][0]
            qE+=readcigar[0][0]
        align_start.append(qS)
        align_end.append(qE)
    query_pos=sorted(zip(align_start,align_end))
    if len(query_pos)==1:
        return(0)
    threshold_long=0.2 * refLen
    threshold_short=0.2 * refLen
    count_long=0
    count_short=0
    for i in range(1,len(query_pos)):
        #if (query_pos[i][0]-query_pos[i-1][1]) < -threshold_long:  # judge if ref longer than real
            #count_long+=1
        if (query_pos[i][0]-query_pos[i-1][1]) > threshold_short: # judge if ref shorter than real
            count_short+=1
    if (count_short/(len(query_pos)-1)>0.5):
        return(0)
    return(1)

def mp_read2ref(tmp,type):
    tmpNrow=tmp.shape[0]
    currentID=tmp.index.tolist()[0]
    readSeq=fq_dict[currentID]
    for i in range(tmpNrow):
        each=tmp.iloc[i,:]
        ref_seq=getSeq(each)
        eachAlign=mp.Aligner(seq=ref_seq,k=10,n_threads=1,preset='map-ont')
        mapAll=eachAlign.map(readSeq)
        read_result=mp_calRead(mapAll,len(ref_seq))
        mapAll.close()
        if read_result:
            return('N',i)
    return(type,-1)
    
def read2ref(tmp,type):
    tmpNrow=tmp.shape[0]
    currentID=tmp.index.tolist()[0]
    fqOutput=fq_dict[currentID]
    foFastq=open(outPrefixTmp+"query_"+currentID+".fq",'w')
    foFastq.write(fqOutput)
    foFastq.close()
    for i in range(tmpNrow):
        each=tmp.iloc[i,:]
        ref_seq=getSeq(each)
        faOutput='>%s\n%s\n' % (currentID,ref_seq)
        foFa=open(outPrefixTmp+"ref_"+currentID+".fa",'w')
        foFa.write(faOutput)
        foFa.close()
        cmd="minimap2 -ax map-ont  "+" -k10 "+outPrefixTmp+"ref_"+currentID+".fa "+outPrefixTmp+"query_"+currentID+".fq >"+outPrefixTmp+currentID+".sam 2>/dev/null"
        os.system(cmd)
        os.remove(outPrefixTmp+"ref_"+currentID+".fa")
        if os.path.exists(outPrefixTmp+currentID+".sam"):
            samfile=pysam.AlignmentFile(outPrefixTmp+currentID+".sam","r")
            read_list=[read  for read in samfile.fetch() if read.flag & 4!=4]
            read_result=calRead(read_list,len(ref_seq))
            os.remove(outPrefixTmp+currentID+".sam")
            if read_result:
                os.remove(outPrefixTmp+"query_"+currentID+".fq")
                return('N',i)
    os.remove(outPrefixTmp+"query_"+currentID+".fq")
    return(type,-1)
        
def judgeFusion_diff(tmp):
    type,rowID=mp_read2ref(tmp,'m')
    if type=='N':
        return(type,rowID)
    tmpChr=list(tmp['chr'])
    tmpChr_unique=list(set(tmpChr))
    tmpNrow=tmp.shape[0]
    tmpChr2Num=[]
    tmpChr_len=len(tmpChr_unique)
    tmp_query_start=tmp['query_start'].values
    tmp_query_end=tmp['query_end'].values
    # in same chromosome but different strand
    tmpStrand=np.array(tmp['strand'])
    if len(set(tmpStrand[np.where(np.array(tmpChr)==tmpChr_unique[0])[0]]))>1 or len(set(tmpStrand[np.where(np.array(tmpChr)==tmpChr_unique[1])[0]]))>1:
        return('U',-1)
    for j in range(tmpNrow):
        if tmpChr[j] == tmpChr_unique[0]:
            tmpChr2Num.append(0)
        else:
            tmpChr2Num.append(1)
    if tmpNrow >=4:
        targetNum=[0,1]
        for j in range(tmpNrow):
            targetNum.extend([0,1])
        targetNum=np.array(targetNum[(tmpChr2Num[0]):(tmpNrow+tmpChr2Num[0])])
        if sum(targetNum == tmpChr2Num) == tmpNrow:
            if (tmp_query_start[2]-tmp_query_end[0])>40 and (tmp_query_start[3]-tmp_query_end[1])>40 and (tmp_query_end[1]-tmp_query_end[0])>40 and (tmp_query_end[2]-tmp_query_end[1])>40 and (tmp_query_end[3]-tmp_query_end[2])>40 and (tmp_query_start[1]-tmp_query_end[0])> -40 and (tmp_query_start[2]-tmp_query_end[1])> -40 and (tmp_query_start[3]-tmp_query_end[2])> -40 :
                return('F1',-1)
    return('FC1',-1)
def judgeFusion_same(tmp):
    tmpStrand=tmp['strand']
    tmpChr=list(tmp['chr'])
    tmpChr_unique=list(set(tmpChr))
    tmpChr2Num=[]
    tmpChr_len=len(tmpChr_unique)
    tmpNrow=tmp.shape[0]
    tmp_query_start=tmp['query_start'].values
    tmp_query_end=tmp['query_end'].values
    if tmpNrow == 2:
        #if len(set(tmpStrand))>1:
            #return('U')
        r1S=tmp['reference_start'].values[0]
        r1E=tmp['reference_end'].values[0]
        r2S=tmp['reference_start'].values[tmpNrow-1]
        r2E=tmp['reference_end'].values[tmpNrow-1]
        r1=Interval(r1S,r1E,lower_closed=False)
        r2=Interval(r2S,r2E,lower_closed=False)
        if r1.overlaps(r2):
            #return('NC',-1)
            return('N',-1) # revised for more sensitive
        type,rowID=mp_read2ref(tmp,'U')
        return(type,rowID)
        
    elif tmpNrow == 3:
        #if len(set(tmpStrand))>1:
            #return('U')
        r1S=tmp['reference_start'].values[0]
        r1E=tmp['reference_end'].values[0]
        r2S=tmp['reference_start'].values[tmpNrow-1]
        r2E=tmp['reference_end'].values[tmpNrow-1]
        r1=Interval(r1S,r1E,lower_closed=False)
        r2=Interval(r2S,r2E,lower_closed=False)
        arr1=np.zeros(tmpNrow-1)
        arr2=np.zeros(tmpNrow-1)
        for j in range(1,tmpNrow):
            cS1=tmp['reference_start'].values[j]
            cE1=tmp['reference_end'].values[j]
            cS2=tmp['reference_start'].values[j-1]
            cE2=tmp['reference_end'].values[j-1]
            c1=Interval(cS1,cE1,lower_closed=False)
            c2=Interval(cS2,cE2,lower_closed=False)
            arr1[j-1]=r1.overlaps(c1)
            arr2[j-1]=r2.overlaps(c2)
        if sum(arr1)>0 and sum(arr2)>0:
            if arr1[0]:
                return('N',-1)
            else:
                type,rowID=mp_read2ref(tmp,'FC2')
                return(type,rowID)
                #return('FC2')
        else:
            type,rowID=mp_read2ref(tmp,'C2')
            return(type,rowID)
    elif tmpNrow>3:
        r1S=tmp['reference_start'].values[np.where(tmp['exon_length'].values==max(tmp['exon_length'].values))[0][0]]
        r1E=tmp['reference_end'].values[np.where(tmp['exon_length'].values==max(tmp['exon_length'].values))[0][0]]
        r1=Interval(r1S,r1E,lower_closed=False)
        arr1=np.zeros(tmpNrow)
        for j in range(0,tmpNrow):
            cS1=tmp['reference_start'].values[j]
            cE1=tmp['reference_end'].values[j]
            c1=Interval(cS1,cE1,lower_closed=False)
            arr1[j]=r1.overlaps(c1)
        #if len(set(tmpStrand[np.array(arr1)==0]))>1 or len(set(tmpStrand[np.array(arr1)==1]))>1:
        #        return('U')
        if sum(arr1 == 1)>(tmpNrow-2):
                return('N',-1)
        else:
            targetNum=[0,1]
            for j in range(2*tmpNrow):
                targetNum.extend([0,1])
            arr1=[int(i) for i in arr1]
            targetNum=np.array(targetNum[arr1[0]:(arr1[0]+tmpNrow)])
            if sum(targetNum == arr1) == tmpNrow:
                if (tmp_query_start[2]-tmp_query_end[0])>40 and (tmp_query_start[3]-tmp_query_end[1])>40 and (tmp_query_end[1]-tmp_query_end[0])>40 and (tmp_query_end[2]-tmp_query_end[1])>40 and (tmp_query_end[3]-tmp_query_end[2])>40 and (tmp_query_start[1]-tmp_query_end[0])> -40 and (tmp_query_start[2]-tmp_query_end[1])> -40 and (tmp_query_start[3]-tmp_query_end[2])> -40 :
                    type,rowID=mp_read2ref(tmp,'F2')
                    return(type,rowID)
                else:
                    type,rowID=mp_read2ref(tmp,'FC2')
                    return(type,rowID)
            else:
                type,rowID=mp_read2ref(tmp,'C2')
                return(type,rowID)
    else:
        return('U',-1)

def getType(i):
    tmp=explainFL_pass1.loc[passID[i]].copy()
    tmpChr=list(tmp['chr'])
    tmpChr_unique=list(set(tmpChr))
    tmpChr_len=len(tmpChr_unique)
    tmpNrow=tmp.shape[0]
    tmp_query_start=tmp['query_start'].values
    tmp_query_end=tmp['query_end'].values
    queryESdiff=float('inf')
    rowID=-1
    for j in range(tmpNrow-1):
        queryESdiff=min(abs(tmp_query_start[j+1]-tmp_query_end[j]),queryESdiff)
    if tmpChr_len == 2:
        tmpChrDict={}
        tmpChrDict[tmpChr_unique[0]]=0
        tmpChrDict[tmpChr_unique[1]]=1
        diff_chr=[]
        for j in range(1,len(tmpChr)):
            diff_chr.append(tmpChrDict[tmpChr[j]]-tmpChrDict[tmpChr[j-1]])
        if sum(np.array(diff_chr)!=0)>1 :
            re,rowID=judgeFusion_diff(tmp)
        else:
            re='C1'
    elif tmpChr_len >2:
        re='M'
    else:
        re,rowID=judgeFusion_same(tmp)
    return([passID[i],re,rowID,queryESdiff])

def filterFL(genomeFile,outPrefix,fastqFile,thread):
    if thread>20:
        thread=20
    global outPrefixTmp, genome,passID,explainFL_pass1,fq_dict
    outPrefixTmp=outPrefix+'tmp/'
    genome=pyfasta.Fasta(genomeFile)
    oldFq=open(fastqFile)
    fileName=outPrefix+'explainFL.txt'
    
    explainFL=pd.read_csv(fileName,header=None,names=["ID","chr","reference_start","reference_end","query_start","query_end","exon_start","exon_end","exon_length","strand","leftSeq","rightSeq"],sep='\t')

    tmp1=explainFL[explainFL['ID'].duplicated()]
    #tmp2=tmp1['ID'].duplicated()
    #tmp3=tmp1[tmp2]
    #passID=list(tmp3['ID'].drop_duplicates())
    passID=list(tmp1['ID'].drop_duplicates())
    explainFL_pass=explainFL[explainFL['ID'].isin(passID)]
    explainFL_pass=explainFL_pass.sort_values(by=["ID","query_start"])
    explainFL_pass1=explainFL_pass.set_index('ID')
    fq_dict={}
    readID_dict={}
    for i in passID:
        readID_dict[i]=1
    while True:
        fq1=oldFq.readline()
        if fq1:
            fq1=fq1.strip()
            fq2=oldFq.readline().strip()
            fq3=oldFq.readline()
            fq4=oldFq.readline().strip()
            ID=fq1.split('\t')[0]
            ID=ID.split(' ')[0][1:]
            if readID_dict.__contains__(ID):
                #fq_dict[ID]=fq1+'\n'+fq2+'\n'+fq3+fq4+'\n'
                fq_dict[ID]=fq2
        else:
            break
    oldFq.close()

    pool=Pool(processes=thread)
    widgets = [ 'getType: ' , Percentage () , ' ' , Bar ( marker = RotatingMarker ( ) ) ,' ' , ETA ( ) , ' '  ]
    bar=progressbar.ProgressBar(widgets=widgets,maxval=len(passID)).start()
    eachBin=20*thread
    allType_list=[]
    for i in range(0,len(passID),eachBin):
        allType_list.extend(pool.map(getType,range(i,min(i+eachBin,len(passID)))))
        bar.update(i)
    
    passID_df=pd.DataFrame(allType_list,columns=['ID','type','rowID','diff'])
    pool.close()
    pool.join()
    passID_df=passID_df.set_index('type')
    passID_df.to_csv(outPrefix+"explainFL_ID2Type.txt",sep="\t")
    passID_df_type=list(set(passID_df.index))
    if 'F1' in passID_df_type:
        explainFL_Fusion1=explainFL_pass1.loc[passID_df.loc[['F1'],'ID']]
        explainFL_Fusion1.to_csv(outPrefix+"explainFL_Fusion1.txt",sep="\t",header=True)
    if 'F2' in passID_df_type:
        explainFL_Fusion2=explainFL_pass1.loc[passID_df.loc[['F2'],'ID']]
        explainFL_Fusion2.to_csv(outPrefix+"explainFL_Fusion2.txt",sep="\t",header=True)
    if 'N' in passID_df_type:
        explainFL_Normal=explainFL_pass1.loc[passID_df.loc[['N'],'ID']]
        passID_dict={}
        geneID=''
        count=0
        n=0
        for i in explainFL_Normal.index.tolist():
            if n==0:
                count=0
            else:
                if i == old:
                    count+=1
                else:
                    count=0
            passID_dict[i+'-'+str(count)]=n
            old=i
            n+=1
        passID_df_n=passID_df.loc[['N'],:].copy()
        N_pass_ID=passID_df_n.loc[passID_df_n.rowID==-1 ,'ID'].tolist()
        passID_df_n_adjust=passID_df_n.loc[passID_df_n.rowID!=-1 ,:]
        FLdf_adj=explainFL_Normal.iloc[[passID_dict[i] for i in (passID_df_n_adjust.ID+'-'+passID_df_n_adjust.rowID.map(str)).tolist()],:]
        explainFL_Normal=pd.concat([explainFL_Normal.loc[N_pass_ID,:],FLdf_adj])
        explainFL_Normal.to_csv(outPrefix+"explainFL_Normal.txt",sep="\t",header=True)
    if 'NC' in passID_df_type:
        explainFL_NC=explainFL_pass1.loc[passID_df.loc['NC','ID']]
        explainFL_NC.to_csv(outPrefix+"explainFL_NC.txt",sep="\t",header=True)

