'''
Usage: circfull partial -f fastq -g genome -c circ -j junc [-t threads] [-o output]

Options:
    -h --help                   Show help message.
    -v --version                Show version.
    -f fastq                    circFL-seq fastq file.
    -g genome                   Fasta file of genome.
    -c circ                     circFL full-length file.
    -j junc                     explainFL file.         
    -t threads                  Number of threads [default: 20].
    -o output                   Output dir [default: circFL_out].
'''
from .genericFun import *
import sys,time,pandas as pd,numpy as np,docopt,os,time,pysam
from .RG_circFL_output import FL2bed

def getTargetRead(explainFL_file):
    explainFL_in=open(explainFL_file)
    target_read=[]
    read2pos={}
    each=explainFL_in.readline()
    while each:
        each_arr=each.strip().split('\t')
        each_next=explainFL_in.readline()
        if not each_next:
            break
        else:
            each_next_arr=each_next.strip().split('\t')
            if each_arr[0]==each_next_arr[0] and each_arr[1]==each_next_arr[1]:
                target_read.append(each_arr[0])                
            each=each_next
            if each_arr[0] not in read2pos:
                read2pos[each_arr[0]]=[[each_arr[1],int(each_arr[2]),int(each_arr[3])]]# additional code to record the orgional mapping position 03/16/24
            else:
                read2pos[each_arr[0]].append([each_arr[1],int(each_arr[2]),int(each_arr[3])])
    if each_arr[0] not in read2pos: # the last one
                read2pos[each_arr[0]]=[[each_arr[1],int(each_arr[2]),int(each_arr[3])]]# additional code to record the orgional mapping position 03/16/24
    else:
        read2pos[each_arr[0]].append([each_arr[1],int(each_arr[2]),int(each_arr[3])])
    explainFL_in.close()
    target_read=list(set(target_read))
    return(target_read,read2pos)

def getTargetFastq(target_read,fastq_file,target_fastq):
    target_read_dict=dict(zip(target_read,[1]*len(target_read)))
    fq=open(fastq_file)
    target_fq=open(target_fastq,'w')
    line1=fq.readline()
    line2=fq.readline()
    line3=fq.readline()
    line4=fq.readline()
    while line1:
        readID=line1.split(' ')[0][1:]
        if target_read_dict.__contains__(readID):
            target_fq.write(line1+line2+line3+line4)
        line1=fq.readline()
        line2=fq.readline()
        line3=fq.readline()
        line4=fq.readline()
    target_fq.close()
    fq.close()
def filterCircFL(df):
    pass1=df.loc[:,'start']>0
    pass2=df.loc[:,'end']>0
    pass3=pass1 & pass2
    return(df.loc[pass3,:])

def fa2twofa(ref_file,out_file):
    ref=open(ref_file)
    faout=open(out_file,'w')
    line1=ref.readline()
    line2=ref.readline()
    while line1:
        faout.write(line1+line2.strip()+line2)
        line1=ref.readline()
        line2=ref.readline()
    faout.close()
    ref.close()
def getReadInfo(read):
    rS=read.reference_start
    rE=read.reference_end
    qS=read.query_alignment_start
    qE=read.query_alignment_end
    readcigar=read.cigar
    if readcigar[0][0] == 5:
        qS+=readcigar[0][1]
        qE+=readcigar[0][1]
    exonS=[rS]
    exonE=[]
    currentPos=rS-1
    i=0
    for each in readcigar:
        i=i+1
        if i == len(readcigar):
            if each[0] in (0,2):
                currentPos+=each[1]
            exonE.append(currentPos+1)
        else:
            if each[0] ==0:
                currentPos+=each[1]
            elif each[0] == 2:
                if each[1]<10: # if insertation more than 10, treat it as N
                    currentPos+=each[1]
                else:
                    exonE.append(currentPos+1)
                    currentPos+=each[1]
                    exonS.append(currentPos+1)
            elif each[0] == 3:
                exonE.append(currentPos+1)
                currentPos+=each[1]
                exonS.append(currentPos+1)
    return([exonS,exonE])


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
        
def judgeIsoOverlap(isoID,pos):
    isoID_arr=isoID.split('|')
    start1=int(isoID_arr[1].split(',')[0])
    start2=int(pos[1])
    end1=int(isoID_arr[2].split(',')[-1])
    end2=int(pos[2])
    if isoID_arr[0]==pos[0]:
        if start1 > end2 or start2 > end1:
            return False
        return True
    else:
        return False

def judgeIsoOverlap_list(isoID,pos_list):
    n=len(pos_list)
    if n<2:
        return(False)
    overlap_judge=[]
    for i in range(n):
        overlap_judge.append(judgeIsoOverlap(isoID,pos_list[i]))
    if sum(overlap_judge)>=2:
        return(True)
    else:
        return(False)
    
def countPartial(sam,circ_df,read2pos,hangLen=5):
    samfile=pysam.AlignmentFile(sam,"r")
    iso2reads_dict={}
    read2iso_dict={}
    readMaxLen_dict={}
    iso2len=dict(zip(circ_df.isoID,circ_df.len))
    for read in samfile.fetch():
        if (read.flag & 4 != 4):
            readInfo=getReadInfo(read)
        else:
            continue
        ExonS=readInfo[0]
        ExonE=readInfo[1]
        iso_len=iso2len[read.reference_name]
        nExon=len(ExonS)
        aligned_len=sum([(ExonE[i]-ExonS[i]) for i in range(nExon)])
        isPass=False
        for i in range(nExon):
            if (ExonS[i]<(iso_len-hangLen)) and (ExonE[i]>(iso_len+hangLen)) :
                targetPos=[iso_len-read.reference_start]
                if judgeJun(read,targetPos):
                    isPass=True           
        if isPass:
            if readMaxLen_dict.__contains__(read.query_name):
                if readMaxLen_dict[read.query_name].__contains__(read.reference_name):
                    readMaxLen_dict[read.query_name][read.reference_name]=max([readMaxLen_dict[read.query_name][read.reference_name],aligned_len/nExon])
                else:
                    readMaxLen_dict[read.query_name][read.reference_name]=aligned_len/nExon
            else:
                readMaxLen_dict[read.query_name]={}
                readMaxLen_dict[read.query_name][read.reference_name]=aligned_len/nExon
    for read in readMaxLen_dict.keys():
        value=readMaxLen_dict[read]
        maxLen=0
        target_ref=''
        for ref in value.keys():
            eachLen=value[ref]
            if eachLen>maxLen:
                target_ref=ref
                maxLen=eachLen
        read2iso_dict[read]=target_ref
        if iso2reads_dict.__contains__(target_ref):
            iso2reads_dict[target_ref].append(read)
        else:
            iso2reads_dict[target_ref]=[read]
    circ_df_read2iso={}
    for i in range(circ_df.shape[0]):
        each=circ_df.iloc[i,:]
        readID=each['readID'].split(',')
        for j in readID:
            circ_df_read2iso[j]=each['isoID']
    for i in read2iso_dict.keys():
        if circ_df_read2iso.__contains__(i):
            continue
        else:
            if i in read2pos:
                read_refpos=read2pos[i] # additional code 03/16/2024
                isoID=read2iso_dict[i] #
                if judgeIsoOverlap_list(isoID,read_refpos): #
                    circ_df_read2iso[i]=read2iso_dict[i]
    iso2read_circ_df={}
    for read,iso in circ_df_read2iso.items():
        if iso2read_circ_df.__contains__(iso):
            iso2read_circ_df[iso].append(read)
        else:
            iso2read_circ_df[iso]=[read]
    new_circ_df=circ_df.copy()
    new_reads=[]
    new_readCount=[]
    for i in new_circ_df.isoID.tolist():
        tread=iso2read_circ_df[i]
        new_readCount.append(len(tread))
        new_reads.append(','.join(tread))
    new_circ_df.readCount=new_readCount
    new_circ_df.readID=new_reads
    return(new_circ_df)

def partial(options):
    circ_file=options['-c']
    explainFL_file=options['-j']
    fastq=options['-f']
    genome=options['-g']
    thread=int(options['-t'])
    outDir=options['-o']
    plog('Check circFL file')
    fileCheck(circ_file)
    plog('Check circFL file')
    fileCheck(explainFL_file)
    plog('Check fastq file')
    checkFastq(fastq)
    plog('Check genome file')
    readFaFile(genome)
    if not os.path.exists(outDir):
        os.mkdir(outDir)
    if outDir[-1]!='/':
        outPrefix=outDir+'/'
    else:
        outPrefix=outDir
    partial_outPrefix=outPrefix+'partial/'
    partial_tmp_outPrefix=partial_outPrefix+'tmp/'
    createDir(outPrefix);createDir(partial_outPrefix);createDir(partial_tmp_outPrefix)
    sam=partial_tmp_outPrefix+'target_reads.sam'
    plog("Read circFL file")
    circ_df=pd.read_csv(circ_file,sep='\t')
    circ_df=filterCircFL(circ_df)
    plog("Get target reads")
    target_read,read2pos=getTargetRead(explainFL_file)
    getTargetFastq(target_read,fastq,partial_tmp_outPrefix+'target_reads.fastq')
    plog("BED to Fasta file")
    circ_df_bed=FL2bed(circ_df)
    circ_df_bed.to_csv(partial_outPrefix+'circFL_Normal.bed',sep='\t',header=None,index=None)
    faFile=partial_tmp_outPrefix+'circFL_Normal.fa'
    os.system("bedtools getfasta  -split -nameOnly -fi %s -bed %s -fo %s" % (genome,partial_outPrefix+'circFL_Normal.bed',faFile))
    fileCheck(faFile)
    twoFL_file=partial_tmp_outPrefix+'twoFL.fa'
    fa2twofa(partial_tmp_outPrefix+'circFL_Normal.fa',twoFL_file)
    plog("Align target reads")
    os.system("minimap2 -ax splice -ub -k14 -w 4 -t %d %s %s >%s" % (thread,twoFL_file,partial_tmp_outPrefix+'target_reads.fastq',sam))
    plog("Count partial reads")
    new_df=countPartial(sam,circ_df,read2pos)
    new_df.to_csv(partial_outPrefix+'circFL_Normal.txt',sep='\t',index=None)
    plog('All done!!!')
    
