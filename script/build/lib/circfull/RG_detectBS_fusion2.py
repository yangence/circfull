import os, sys, pyfasta, pysam, pandas as pd, numpy as np
from multiprocessing import Pool
hangLen=150


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
                currentPos+=each[1]
            elif each[0] == 3:
                exonE.append(currentPos+1)
                currentPos+=each[1]
                exonS.append(currentPos+1)
    exonLen=sum([exonE[i]-exonS[i] for i in range(len(exonS))])
    return([exonS,exonE])
def getEachPos(rlist,type='type1',strand=True):
    [r1,r2,r3,r4]=rlist
    currentID=r1.split(' ')[0][1:]
    fastqOutput="%s%s%s%s" % (r1,r2,r3,r4)
    eachFaPos=faPos[currentID+'_'+type]
    eachChr1=eachFaPos[0]
    eachChr2=eachFaPos[1]
    eachLeft=int(eachFaPos[2])
    eachRight=int(eachFaPos[3])
    eachLen=int(eachFaPos[4])
    eachFaSeq=faOutput[currentID+'_'+type]
    fqFile=outPrefixTmp+"seq1_"+type+'_'+currentID+".fastq"
    foFastq=open(fqFile,'w')
    foFastq.write(fastqOutput)
    foFastq.close()
    faFile=outPrefixTmp+"seq2_"+type+'_'+currentID+".fa"
    foFa=open(faFile,'w')
    foFa.write(eachFaSeq)
    foFa.close()
    samFile=outPrefixTmp+currentID+"_"+type+".sam"
    cmd="minimap2 -ax splice -k14 " +isSecond+" "+faFile+" "+fqFile+" >"+samFile+" 2>/dev/null"
    #cmd="minimap2 -ax splice -k14 "+faFile+" "+fqFile+" >"+samFile+" 2>/dev/null"
    os.system(cmd)
    samfile=pysam.AlignmentFile(samFile,"r")
    BSright=[]
    Mright=[]
    BSleft=[]
    Mleft=[]
    for read in samfile.fetch():
        if read.flag & 4 != 4:
            readInfo=getReadInfo(read)
        else:
            continue
        ExonS=readInfo[0]
        ExonE=readInfo[1]
        ExonS_diff=abs(np.array(ExonS)-eachLen-hangLen)
        ExonE_diff=abs(np.array(ExonE)-eachLen+hangLen)
        ExonS_diff_idx=np.where(ExonS_diff==min(ExonS_diff))[0]
        ExonE_diff_idx=np.where(ExonE_diff==min(ExonE_diff))[0]
        commonIdx=set(ExonS_diff_idx-1) & set(ExonE_diff_idx)
        if len(commonIdx)==0:
            if len(ExonS_diff_idx)==1 and len(ExonE_diff_idx)==1:
                if ExonS_diff_idx[0] == ExonE_diff_idx[0]:
                    if ExonS[ExonS_diff_idx[0]]<eachLen:
                        if  ExonS_diff_idx[0]<len(ExonS)-1:
                            if ExonS[ExonS_diff_idx[0]+1]>eachLen:
                                commonIdx=[ExonE_diff_idx[0]]
                    else:
                        if  ExonE_diff_idx[0]>0:
                            if ExonE[ExonE_diff_idx[0]-1]<eachLen:
                                commonIdx=[ExonE_diff_idx[0]-1]
        for index in commonIdx:
        #1 based position
            if strand:
                tmpright=eachLeft+ExonE[index]-1
                tmpleft=eachRight+ExonS[index+1]-eachLen
                BSright.append(tmpright)
                BSleft.append(tmpleft)
                Mright.append(genome.sequence({'chr': eachChr1, 'start':tmpright+1, 'stop':tmpright+2}).upper())
                Mleft.append(genome.sequence({'chr': eachChr2, 'start':tmpleft-2, 'stop':tmpleft-1}).upper())
            else:
                if type == 'type1':
                    tmpright=eachLeft+ExonE[index]-1
                    tmpleft=eachRight-(ExonS[index+1]-eachLen)
                    BSright.append(tmpright)
                    BSleft.append(tmpleft)
                    Mright.append(genome.sequence({'chr': eachChr1, 'start':tmpright+1, 'stop':tmpright+2}).upper())
                    Mleft.append(genome.sequence({'chr': eachChr2, 'start':tmpleft+1, 'stop':tmpleft+2,'strand':'-'}).upper())
                else:
                    tmpright=eachLeft-ExonE[index]
                    tmpleft=eachRight+ExonS[index+1]-eachLen
                    BSright.append(tmpright)
                    BSleft.append(tmpleft)
                    Mright.append(genome.sequence({'chr': eachChr1, 'start':tmpright-2, 'stop':tmpright-1,'strand':'-'}).upper())
                    Mleft.append(genome.sequence({'chr': eachChr2, 'start':tmpleft-2, 'stop':tmpleft-1}).upper())
    samfile.close()
    os.remove(fqFile)
    os.remove(faFile)
    os.remove(samFile)
    #print([eachChr1,eachChr2,BSleft,BSright,Mleft,Mright])
    return([eachChr1,eachChr2,BSleft,BSright,Mleft,Mright])

def getBs_fusion1(rlist):
    [r1,r2,r3,r4]=rlist
    currentID=r1.split(' ')[0][1:]
    [eachChr1_first,eachChr2_first,BSleft_first,BSright_first,Mleft_first,Mright_first]=getEachPos(rlist,'type1',faStrand[currentID])
    [eachChr1_second,eachChr2_second,BSleft_second,BSright_second,Mleft_second,Mright_second]=getEachPos(rlist,'type2',faStrand[currentID])
    tmpOut="%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (currentID,
    eachChr1_first,
    eachChr1_second,
    ','.join([str(i) for i in BSright_first]),
    ','.join([str(i) for i in BSleft_first]),
    ','.join([str(i) for i in Mright_first]),
    ','.join([str(i) for i in Mleft_first]),
    eachChr1_second,
    eachChr2_second,
    ','.join([str(i) for i in BSright_second]),
    ','.join([str(i) for i in BSleft_second]),
    ','.join([str(i) for i in Mright_second]),
    ','.join([str(i) for i in Mleft_second]),
    str(faStrand[currentID]))
    return(tmpOut)
    
def detectBS_fusion2(options):
    global isSecond,faPos,faStrand,faOutput,outPrefixTmp,genome
    genomeFile=options[0]
    outPrefix=options[1]
    thread=options[2]
    outPrefixTmp=outPrefix+'tmp/'
    isSecond=''
    if len(options)>3:
        isSecond='-uf'
    fastqFile=open(outPrefix+'fusion1.fq')
    FLdf_fusion2=pd.read_csv(outPrefix+"explainFL_Fusion2.txt",sep='\t')
    FLdf_fusion2=FLdf_fusion2.set_index('ID')
    FLdf_fusion2['exon_start']=FLdf_fusion2['exon_start'].map(str)
    FLdf_fusion2['exon_end']=FLdf_fusion2['exon_end'].map(str)
    FLdf_fusion2_ID=list(set(FLdf_fusion2.index))
    genome = pyfasta.Fasta(genomeFile)

    faOutput={}
    faPos={}
    faStrand={}
    for i in range(len(FLdf_fusion2_ID)):
        each_df=FLdf_fusion2.loc[FLdf_fusion2_ID[i]]
        # one chromosome
        tmpIdx1=0
        tmpIdx2=1
        if each_df.iloc[0,2]>each_df.iloc[1,2]:
            tmpIdx1=1
            tmpIdx2=0
        each=each_df.iloc[range(tmpIdx1,each_df.shape[0],2),:].sort_values('exon_length',ascending=False).iloc[0,]
        chr_first=each['chr']
        exonS_first=[int(j) for j in each['exon_start'].split(',')]
        exonE_first=[int(j) for j in each['exon_end'].split(',')]
        strand_first=each['strand']
        # another chromosome
        each=each_df.iloc[range(tmpIdx2,each_df.shape[0],2),:].sort_values('exon_length',ascending=False).iloc[0,]
        chr_second=each['chr']
        exonS_second=[int(j) for j in each['exon_start'].split(',')]
        exonE_second=[int(j) for j in each['exon_end'].split(',')]
        strand_second=each['strand']
        type1_first=genome.sequence({'chr': chr_first, 'start':exonS_first[0]+1-hangLen, 'stop':exonE_first[-1]+hangLen})
        type2_first=type1_first
        if strand_first == strand_second:
            faStrand[FLdf_fusion2_ID[i]]=True
            # Type 1 first end + second start
            type1_second=genome.sequence({'chr': chr_second, 'start':exonS_second[0]+1-hangLen, 'stop':exonE_second[-1]+hangLen})
            type1_BS=type1_first+type1_second
            type1_name=">type1"
            faOutput[FLdf_fusion2_ID[i]+'_type1']="%s\n%s" % (type1_name, type1_BS)
            faPos[FLdf_fusion2_ID[i]+'_type1']= [chr_first,chr_second,exonS_first[0]+1-hangLen,exonS_second[0]+1-hangLen,exonE_first[-1]-exonS_first[0]+2*hangLen]
            # Type 2 second end + first start
            type2_second=type1_second
            type2_BS=type2_second + type2_first
            type2_name=">type2"
            faOutput[FLdf_fusion2_ID[i]+'_type2']="%s\n%s" % (type2_name, type2_BS)
            faPos[FLdf_fusion2_ID[i]+'_type2']= [chr_second,chr_first,exonS_second[0]+1-hangLen,exonS_first[0]+1-hangLen,exonE_second[-1]-exonS_second[0]+2*hangLen]
        else:
            faStrand[FLdf_fusion2_ID[i]]=False
            # Type 1 first end + second end
            type1_second=genome.sequence({'chr': chr_second,  'start':exonS_second[0]+1-hangLen, 'stop':exonE_second[-1]+hangLen,'strand': '-'})
            type1_BS=type1_first+type1_second
            type1_name=">type1"
            faOutput[FLdf_fusion2_ID[i]+'_type1']="%s\n%s" % (type1_name, type1_BS)
            faPos[FLdf_fusion2_ID[i]+'_type1']= [chr_first,chr_second,exonS_first[0]+1-hangLen,exonE_second[-1]+hangLen,exonE_first[-1]-exonS_first[0]+2*hangLen]
            # Type 2 second start + first start
            type2_second=type1_second
            type2_BS=type2_second + type2_first
            type2_name=">type2"
            faOutput[FLdf_fusion2_ID[i]+'_type2']="%s\n%s" % (type2_name, type2_BS)
            faPos[FLdf_fusion2_ID[i]+'_type2']= [chr_second,chr_first,exonE_second[-1]+hangLen,exonS_first[0]+1-hangLen,exonE_second[-1]-exonS_second[0]+2*hangLen]

    targetKey=dict(zip(FLdf_fusion2_ID,[1 for i in range(len(FLdf_fusion2_ID))]))
    fastqList=[]
    while True:
        r1=fastqFile.readline()
        if r1:  
            r2=fastqFile.readline()
            r3=fastqFile.readline()
            r4=fastqFile.readline().strip('\n')
            currentID=r1.split(' ')[0][1:]
            if targetKey.__contains__(currentID):
                fastqList.append([r1,r2,r3,r4])
        else:
            break
    fastqFile.close()
    pool=Pool(processes=thread)
    result=pool.map(getBs_fusion1,fastqList)
    pool.close()
    pool.join()
    if len(result)>1:
        fout=open(outPrefix+'BS_Fusion2.txt','w')
        fout.write("\n".join(result)+'\n')
        fout.close()
        return(True)
    else:
        return(False)
