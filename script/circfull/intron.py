'''
Usage: circfull intron -f fastq -g genome -a anno -c circ -i intron [-o output]

Options:
    -h --help                   Show help message.
    -v --version                Show version.
    -f fastq                    circFL-seq fastq file.
    -g genome                   Fasta file of genome.
    -a anno                     Tabix indexed gtf file of gene annotation.
    -c circ                     circFL full-length file.
    -i intron                   Intron candiates, output of circfull anno.
    -o output                   Output dir [default: circFL_out].
'''
from .genericFun import *
import sys,pandas as pd,pysam,pyfasta, mappy as mp,numpy as np

def noDup_list(x):
    newList=[]
    xDict=dict()
    for i in x:
        tmp=str(i[0])+'_'+str(i[1])
        if xDict.__contains__(tmp):
            xDict[tmp]+=1
        else:
            xDict[tmp]=1
            newList.append(i)
    return(newList)

class geneInform:
     def __init__(self, gtfFile):
        tabixfile=pysam.TabixFile(gtfFile)
        gene2strand_dict={}
        gene2trans_dict={}
        gene2exon_dict={}
        trans2exon_dict={}
        trans2gene_dict={}
        gene2status_dict={}
        gene2class_dict={}
        trans2class_dict={}  
        gene2chr_dict={}

        for gtf in tabixfile.fetch(parser=pysam.asGTF()):
            if gtf.feature=='exon':
                current_geneID=gtf.gene_id
                current_transID=gtf.transcript_id
                if gene2exon_dict.__contains__(current_geneID):
                    gene2exon_dict[current_geneID].append([gtf.start,gtf.end])
                else:
                    gene2exon_dict[current_geneID]=[[gtf.start,gtf.end]]
                if trans2exon_dict.__contains__(current_transID):
                    trans2exon_dict[current_transID].append([gtf.start,gtf.end])
                else:
                    trans2exon_dict[current_transID]=[[gtf.start,gtf.end]]
            if gtf.feature=='transcript':
                current_geneID=gtf.gene_id
                current_transID=gtf.transcript_id
                if 'transcript_status' in gtf.keys():
                    trans2class_dict[current_transID]=gtf.transcript_status
                else:
                    trans2class_dict[current_transID]='NOVEL'
                if gene2trans_dict.__contains__(current_geneID):
                    gene2trans_dict[current_geneID].append(current_transID)
                else:
                    gene2trans_dict[current_geneID]=[current_transID]
                trans2gene_dict[current_transID]=current_geneID
            if gtf.feature=='gene':
                current_geneID=gtf.gene_id
                gene2strand_dict[current_geneID]=gtf.strand
                gene2chr_dict[current_geneID]=gtf.contig
                if 'gene_status' in gtf.keys():
                    gene2status_dict[current_geneID]=gtf.gene_status
                else:
                    gene2status_dict[current_geneID]='NOVEL'
                gene2class_dict[current_geneID]=gtf.gene_name

        class2gene_dict=dict(zip(gene2class_dict.values(),gene2class_dict.keys()))

        for i in gene2exon_dict.keys():
            gene2exon_dict[i]=noDup_list(gene2exon_dict[i])
            
        self.gene2strand_dict=gene2strand_dict
        self.gene2trans_dict=gene2trans_dict
        self.gene2exon_dict=gene2exon_dict
        self.trans2exon_dict=trans2exon_dict
        self.trans2gene_dict=trans2gene_dict
        self.gene2status_dict=gene2status_dict
        self.gene2class_dict=gene2class_dict
        self.trans2class_dict=trans2class_dict
        self.gene2chr_dict=gene2chr_dict

def createIntrondict(gene2trans,trans2exon):
    def find_intron(exonlist):
        exonNum=len(exonlist)
        intronlist=[]
        if exonNum>1:
            for i in range(exonNum-1):
                start=exonlist[i][1]+1 # 1-based position
                end=exonlist[i+1][0] # 1-based position
                intronlist.append([start,end])
            return(intronlist)
        else:
            return(intronlist)
    gene2intron_dict={}
    for geneID in gene2trans.keys():
        transID_arr=gene2trans[geneID]
        for transID in transID_arr:
            exon_list=trans2exon[transID]
            intron_list=find_intron(exon_list)
            if len(intron_list)>0:
                if gene2intron_dict.__contains__(geneID):            
                    gene2intron_dict[geneID].extend(intron_list)
                else:
                    gene2intron_dict[geneID]=intron_list
    for i in gene2intron_dict.keys():
             gene2intron_dict[i]=noDup_list(gene2intron_dict[i])
    return(gene2intron_dict)


def intronOverlap(gene_df,targetDf,hangLen=20):
    def judgeOverlap(intronarr,start,end,geneStrand,hangLen):
        overlaplist=[]
        for each in intronarr:
            if geneStrand=='+':
                if (start <= each[0]+hangLen) and (start>=each[0]-hangLen):
                    overlaplist.append(each)
            else:
                if (end<= each[1]+hangLen) and (end>=each[1]-hangLen):
                    overlaplist.append(each)
        return(overlaplist)
    
    intron2isoID={}
    for geneID, intron_arr in gene_df.gene2intron.items():
        geneName=gene_df.gene2class_dict[geneID]
        geneStrand=gene_df.gene2strand_dict[geneID]
        subdf=targetDf.loc[targetDf.geneName==geneName,:]
        if subdf.shape[0]==0:
            continue
        #print(subdf.shape[0])
        for i in range(subdf.shape[0]):
            each_entry=subdf.iloc[i,:]
            each_exon_start=int(each_entry['exon_start'])
            each_exon_end=int(each_entry['exon_end'])
            chro=each_entry['chr']
            isoID=each_entry['isoID']
            overlap_intron=judgeOverlap(intron_arr,each_exon_start,each_exon_end,geneStrand,hangLen)
            for j in overlap_intron:
                intron_key=chro+'|'+str(j[0])+'|'+str(j[1]) #0-based
                if intron2isoID.__contains__(intron_key):
                    intron2isoID[intron_key].append(isoID)
                else:
                    intron2isoID[intron_key]=[isoID]
    return(intron2isoID)

def readfastq(RG_df,fastq_file):
    read2fq={}
    iso2readID={}
    for i in range(RG_df.shape[0]):
        each=RG_df.iloc[i,:]
        isoID=each['isoID']
        readname=each['readID'].split(',')
        iso2readID[isoID]=readname
    require_reads=[]
    for i in  iso2readID.values():
        for j in i:
            require_reads.append(j)
    #print(require_reads)
    with open(fastq_file) as fin:
        line1=fin.readline()
        while line1:
            line2=fin.readline()
            line3=fin.readline()
            line4=fin.readline()
            readID=line1.split(' ')[0][1:]
            if readID in require_reads:
                read2fq[readID]=line2.strip()
            line1=fin.readline()
    return(iso2readID,read2fq)

def mp_calRead_intron(mapAll,refLen):
    def checkInsertion(insertionArr,refLen):
        for i in insertionArr:
            if refLen+5>=i[0] and refLen-5<=i[0]:
                return(False)
        return(True)
    isPass=False
    isLariat=False
    adjustPos=[0,0]
    for read in mapAll:
        rS=read.r_st# 0-based
        rE=read.r_en
        qS=read.q_st
        qE=read.q_en
        readcigar=read.cigar
        #print(read.cigar_str)
        if readcigar[0][0] == 5:
            qS+=readcigar[0][1]
            qE+=readcigar[0][1]
        posS=[rS]
        posE=[]
        currentPos=rS
        insertPos=[]
        delPos=[]
        i=0
        lastCIGAR=0
        for each in readcigar:
            i=i+1
            if i == len(readcigar):
                if each[1] ==0:
                    currentPos+=each[0]
                    if len(posS)>len(posE):
                        posE.append(currentPos)
                    else:
                        posE[-1]=currentPos
            else:
                if each[1] ==0:
                    currentPos+=each[0]
                    if len(posS)>len(posE):
                        posE.append(currentPos)
                    else:
                        posE[-1]=currentPos
                    lastCIGAR=0
                elif each[1] ==2:
                    currentPos+=each[0]
                    posS.append(currentPos)
                    lastCIGAR=2
                    delPos.append([currentPos,each[0]])
                elif each[1]==1:
                    lastCIGAR=1
                    insertPos.append([currentPos,each[0]])
                elif each[1] == 3:
                    posS=[]
                    posE=[]
                    continue
        #print([posS,posE])
        #print(insertPos)
        #print(delPos)
        matchIntron=False
        for i in range(len(posS)):
            each_S=posS[i]
            each_E=posE[i]
            if (each_S+1+5<refLen) and (each_E-5>refLen):
                matchIntron=True
                #isPass=True
                if checkInsertion(insertPos,refLen):
                    isPass=True
                    break
                else:
                    break
        isLariat=False
        adjustPos=[0,0]
        if not matchIntron:
            for j in range(len(delPos)):
                delPos_site=delPos[j][0]
                delPos_len=delPos[j][1]
                delPos_site_lastEnd=delPos_site-delPos_len
                if delPos_site_lastEnd==refLen:
                    adjustPos=[delPos_len,0] #ajdust left site
                    isLariat=True
                    return([isPass,isLariat,adjustPos])
                elif delPos_site==refLen:
                    adjustPos=[0,delPos_len] #adjust right site
                    isLariat=True
                    return([isPass,isLariat,adjustPos])
        if isPass:
            return([isPass,isLariat,[0,0]])
            
    return([isPass,isLariat,adjustPos])
    

def getSeq(each):
    exonStart=[int(i) for i in each['exon_start'].split(',')]
    exonEnd=[int(i) for i in each['exon_end'].split(',')]
    chr=each['chr']
    seq=''
    for i in range(len(exonStart)):
        seq+=genome.sequence({'chr': chr, 'start':exonStart[i], 'stop':exonEnd[i]})
    return(seq)

def revCom(seq):
    seqnew=seq.lower()
    seqnew=seqnew.replace('a','T').replace('t','A').replace('g','C').replace('c','G')[::-1]
    return(seqnew)

def mp_read2intronref(readID,read2fastq,intronKey):
    fastqlist=[read2fastq[i] for i in readID]
    intronpos=intronKey.split('|')
    introndict={'chr':intronpos[0],'exon_start':intronpos[1],'exon_end':intronpos[2]}
    ref_seq=getSeq(introndict)
    #print(len(ref_seq))
    #print(ref_seq)
    ref_seq+=ref_seq
    judge_dict={}
    for i in range(len(fastqlist)):
        readSeq=fastqlist[i]
        eachAlign=mp.Aligner(seq=ref_seq,k=10,n_threads=1,preset='map-ont')
        mapAll=eachAlign.map(readSeq)
        #print(readID[i])
        #print(readSeq)
        read_result=mp_calRead_intron(mapAll,len(ref_seq)/2)
        mapAll.close()
        judge_dict[readID[i]]=read_result
    return(judge_dict)


def combineJudge(readJudge):
    merge_iso2intron={}
    merge_iso2num={}
    def ligateorlariat(intronKey,arr):
        if arr[0]:
            return(['ligate',intronKey,intronKey,0])
        elif arr[1]:
            intronlist=intronKey.split('|')
            chro=intronlist[0]
            start=int(intronlist[1])+arr[2][0]
            end=int(intronlist[2])-arr[2][1]
            lariatKey="{}|{}|{}".format(chro,start,end)
            return(['lariat',lariatKey,intronKey,arr[2][0]+arr[2][1]])
        else:
            return(['','',intronKey,0])
    readID2type={}    
    for isoID, detail in readJudge.items():
        for intronKey, readResult in detail.items():
            for readID, readDetail in readResult.items():
                if readID2type.__contains__(readID):
                    if readID2type[readID][0]!='ligate':
                        tmp_read2type=ligateorlariat(intronKey,readDetail)
                        if tmp_read2type[0]=='ligate':
                            readID2type[readID]==tmp_read2type
                        elif tmp_read2type[0]=='lariat':
                            if tmp_read2type[3]<readID2type[readID][3]:
                                readID2type[readID]==tmp_read2type
                else:
                    tmp_read2type=ligateorlariat(intronKey,readDetail)
                    if tmp_read2type[0]!='':
                        readID2type[readID]=ligateorlariat(intronKey,readDetail)
            
    return(readID2type)
                
def prepareFormat(df):
    iso2reads={}
    iso2num={}
    for i in range(df.shape[0]):
        tmp=df.iloc[i,:]
        if iso2reads.__contains__(tmp[1]):
            iso2reads[tmp[1]].append(tmp.name)
        else:
            iso2reads[tmp[1]]=[tmp.name]
    for i in iso2reads.keys():
        iso2num[i]=len(iso2reads[i])
        iso2reads[i]=','.join(iso2reads[i])
    tmp=df.loc[~df.iloc[:,1].duplicated(),:]
    tmp.index=tmp.iloc[:,1]
    subtmp=tmp.loc[iso2num.keys(),:]
    ndf=pd.DataFrame({'isoID':iso2num.keys(),'intronID':subtmp.iloc[:,2],'type':subtmp.iloc[:,0],'distance':subtmp.iloc[:,3],'readCount':iso2num.values(),'readID':iso2reads.values()})   
    ndf=ndf.sort_values('readCount',ascending=False)
    return(ndf)

def intronFormat2RGformat(df,gene_df):
    circID=df.isoID
    circID_arr=[i.split('|') for i in circID]
    chro=[i[0]  for i in circID_arr]
    start=[int(i[1])  for i in circID_arr]
    end=[int(i[2])  for i in circID_arr]
    circlen=list(np.array(end)-np.array(end)+1)
    exonNum=[1] * len(start)
    leftSeq_list=[]
    for i in range(len(start)):
        leftSeq_list.append(genome.sequence({'chr': chro[i], 'start':start[i]-2, 'stop':start[i]-1}))
    rightSeq_list=[]
    for i in range(len(start)):
        rightSeq_list.append(genome.sequence({'chr': chro[i], 'start':end[i]+1, 'stop':end[i]+2}))
    motif_list=[]
    for i in range(len(start)):
        motif_list.append(leftSeq_list[i]+rightSeq_list[i])
    geneID=[gene_df.intron2gene[i] for i in df.intronID.tolist()]
    geneName=[gene_df.gene2class_dict[i] for i in geneID]
    strand=[gene_df.gene2strand_dict[i] for i in geneID]
    RG_format_df=pd.DataFrame({'circID':circID,'isoID':circID,'chr':chro,'start':start,'end':end,'len':circlen,'exonNum':exonNum,
                             'exon_start':start,'exon_end':end,'motif':motif_list,'leftSeq':leftSeq_list,'rightSeq':rightSeq_list,
                             'exon_leftSeq':leftSeq_list,'exon_rightSeq':rightSeq_list,'strand':strand,'geneName':geneName,
                              'readCount':df.readCount.tolist(),'readID':df.readID.tolist()})
    return(RG_format_df)
            
def intron(options):
    global genome
    circ_file=options['-c']
    fastq=options['-f']
    genomeFile=options['-g']
    anno=options['-a']
    intron_canFile=options['-i']
    outDir=options['-o']
    plog('Check circFL file')
    fileCheck(circ_file)
    plog('Check fastq file')
    checkFastq(fastq)
    plog('Check genome file')
    readFaFile(genomeFile)
    genome=pyfasta.Fasta(genomeFile)
    plog('Check anno file')
    readGTFfile(anno)
    if not os.path.exists(outDir):
        os.mkdir(outDir)
    if outDir[-1]!='/':
        outPrefix=outDir+'/'
    else:
        outPrefix=outDir
    intron_outPrefix=outPrefix+'intron/'
    createDir(outPrefix);createDir(intron_outPrefix)
    plog("Read circFL file")
    RG_df=pd.read_csv(circ_file,sep='\t')
    plog("Prepare intron candidates")
    intron_can=pd.read_csv(intron_canFile,sep='\t')
    gene_df=geneInform(anno)
    gene2intron=createIntrondict(gene_df.gene2trans_dict,gene_df.trans2exon_dict)

    gene_df.gene2intron=gene2intron

    intron2gene={}

    for geneID,intron_arr in gene2intron.items():
        chro=gene_df.gene2chr_dict[geneID]
        for i in intron_arr:
            intronKey=chro+'|'+str(i[0])+'|'+str(i[1])
            intron2gene[intronKey]=geneID

    gene_df.intron2gene=intron2gene
    intron_can.exon_start=intron_can.exon_start.map(str)
    intron_can.exon_end=intron_can.exon_end.map(str)
    intron_can_pass=intron_can.loc[intron_can.exon_start.map(lambda x: len(x.split(',')))==1,:] # 1-based start; only single exon
    #intron_can_pass=intron_can.loc[intron_can.exon_start.map(lambda x: len(x.split(',')))>0,:] # 1-based start; multiple exon
    target_intron2isoID=intronOverlap(gene_df,intron_can_pass,hangLen=20)

    target_iso2intron={}

    for intron,i in target_intron2isoID.items():
        for j in i:
            if target_iso2intron.__contains__(j):
                target_iso2intron[j].append(intron)
            else:
                target_iso2intron[j]=[intron]

    RG_df_sub=RG_df.loc[RG_df.isoID.isin(target_iso2intron.keys()),:]
    
    plog("Read fastq file")
    RG_iso2readID,RG_read2fastq=readfastq(RG_df_sub,fastq)
    
    plog("Intron-aware alignment")
    readsJudge_iso2intron={}
    for isoID, intron_arr in target_iso2intron.items():
        if RG_iso2readID.__contains__(isoID):
            readsJudge_iso2intron[isoID]={}
            readID_arr=RG_iso2readID[isoID]
            for each_intron in intron_arr:
                #chr10|63599276|63599665': {'chr10|63599284|63599669
                readsJudge_iso2intron[isoID][each_intron]=mp_read2intronref(readID_arr,RG_read2fastq,each_intron)
    
    plog("Prepare output files")          
    merge_reads2type=combineJudge(readsJudge_iso2intron)
    merge_reads2type_df=pd.DataFrame(merge_reads2type).T
    merge_reads2type_df_dict=dict(zip(merge_reads2type_df.index,merge_reads2type_df.iloc[:,1]))    
    merge_reads2type_df_format=prepareFormat(merge_reads2type_df)

    RG_df_sub_pass_idx=[]
    for i in range(RG_df_sub.shape[0]):
        each=RG_df_sub.iloc[i,:]
        readID=each['readID'].split(',')
        if len(set(readID) & set(merge_reads2type_df_dict.keys()))>0:
            RG_df_sub_pass_idx.append(i)

    RG_df_sub_pass=RG_df_sub.iloc[RG_df_sub_pass_idx,:]
    RG_df_1=RG_df.loc[~RG_df.isoID.isin(RG_df_sub_pass.isoID),:]
    RG_df_2=intronFormat2RGformat(merge_reads2type_df_format,gene_df)
    RG_df_2_detail=RG_df_2.copy()
    RG_df_2_detail.loc[:,'intronID']=merge_reads2type_df_format.loc[:,'intronID'].tolist()
    RG_df_2_detail.loc[:,'type']=merge_reads2type_df_format.loc[:,'type'].tolist()
    RG_df_2_detail.loc[:,'distance']=merge_reads2type_df_format.loc[:,'distance'].tolist()

    RG_df_3=pd.concat([RG_df_1.iloc[:,0:18],RG_df_2])
    RG_df_3.to_csv(intron_outPrefix+'circFL_Normal.txt',sep='\t',index=None)
    RG_df_2_detail.to_csv(intron_outPrefix+'circFL_intron.txt',sep='\t',index=None)
    plog('All done!!!')              
