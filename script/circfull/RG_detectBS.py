import os, sys, pyfasta, pysam, pandas as pd, numpy as np
from multiprocessing import Pool
from progressbar import *
hangLen=300


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

def getBs(rlist):
    [r1,r2,r3,r4]=rlist
    currentID=r1.strip().split(' ')[0][1:]
    fastqOutput="%s%s%s%s" % (r1,r2,r3,r4)
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
    cmd="minimap2 -ax splice "+strandFastq+" -k14 "+outPrefixTmp+"seq2_"+currentID+".fa "+outPrefixTmp+"seq1_"+currentID+".fastq >"+outPrefixTmp+currentID+".sam 2>/dev/null"
    os.system(cmd)
    if not os.path.exists(outPrefixTmp+currentID+".sam"):
        return('')
    samfile=pysam.AlignmentFile(outPrefixTmp+currentID+".sam","r")
    BSright=[]
    Mright=[]
    BSleft=[]
    Mleft=[]
    for read in samfile.fetch():
        if (read.flag & 4 != 4):
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
            tmpright=eachLeft+ExonE[index]-1
            tmpleft=eachLeft+ExonS[index+1]-eachLen
            BSright.append(tmpright)
            BSleft.append(tmpleft)
            Mright.append(genome.sequence({'chr': eachChr, 'start':tmpright+1, 'stop':tmpright+2}).upper())
            Mleft.append(genome.sequence({'chr': eachChr, 'start':tmpleft-2, 'stop':tmpleft-1}).upper())
    samfile.close()
    os.remove(outPrefixTmp+"seq1_"+currentID+".fastq")
    os.remove(outPrefixTmp+"seq2_"+currentID+".fa")
    os.remove(outPrefixTmp+currentID+".sam")
    if len(BSleft)==0:
        return('')
    return(currentID+"\t"+eachChr+"\t"+','.join([str(i) for i in BSleft])+"\t"+','.join([str(i) for i in BSright])+"\t"+','.join([str(i) for i in Mleft])+"\t"+','.join([str(i) for i in Mright]))

def detectBS(options):
    global strandFastq,outPrefixTmp,genome,faPos,faOutput
    genomeFile=options[0]
    outPrefix=options[1]
    fastqFile=options[2]
    thread=int(options[3])
    strandFastq=''
    if len(options)>4:
        strandFastq='-uf'
    outPrefixTmp=outPrefix+'tmp/'
    FLdf=pd.read_csv(outPrefix+"explainFL_Normal_adj.txt",sep='\t')
    #FLdf_NC=pd.read_csv(outPrefix+"explainFL_NC.txt",sep='\t')
    #FLdf=pd.concat([FLdf,FLdf_NC]).copy() # include candidate to increase sensitive
    FLdf=FLdf.sort_values(by=["ID","exon_length"],ascending=False)
    FLdf_counts=FLdf['ID'].value_counts()
    FLdf=FLdf.drop_duplicates('ID')

    genome = pyfasta.Fasta(genomeFile)
    fastqFile=open(fastqFile)
    fout=open(outPrefix+'BS_Normal.txt','w')

    targetID=list(FLdf['ID'])

    targetKey=dict(zip(targetID,[1 for i in range(len(targetID))]))

    fastqList=[]
    while True:
        r1=fastqFile.readline()
        if r1:  
            r2=fastqFile.readline()
            r3=fastqFile.readline()
            r4=fastqFile.readline().strip('\n')
            currentID=r1.strip().split(' ')[0][1:]
            if targetKey.__contains__(currentID):
                fastqList.append([r1,r2,r3,r4])
        else:
            break
    fastqFile.close()
############################
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
    widgets = [ 'getBs: ' , Percentage () , ' ' , Bar ( marker = RotatingMarker ( ) ) ,' ' , ETA ( ) , ' '  ]
    bar=progressbar.ProgressBar(widgets=widgets,maxval=len(fastqList)).start()
    eachBin=20*thread
    allResult_list=[]
    for i in range(0,len(fastqList),eachBin):
        allResult_list.extend(pool.map(getBs,fastqList[i:min(i+eachBin,len(fastqList))]))
        bar.update(i)
    result=np.array(allResult_list)
    result=result[np.where(result!='')[0]]
    pool.close()
    pool.join()
    for i in result:
        fout.write(i+"\n")
    fout.close()
