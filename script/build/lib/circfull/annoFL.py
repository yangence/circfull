'''
Usage: circfull anno -b bed -a anno [-o output]

Options:
    -h --help                   Show help message.
    -v --version                Show version.
    -b bed                      Full-length circRNA position information in BED format (12 columns).
    -a anno                     Tabix indexed gtf file of gene annotation.
    -o output                   Output file [default: circ_anno.txt].
'''

import pandas as pd
import pysam
from .genericFun import *

def bed2FL(x):
    chr=x[0]
    start=x[1]+1
    end=x[2]
    exon_len=[int(i) for i in x[10].split(',')]
    exon_start_pos=[int(i) for i in x[11].split(',')]
    exon_start=[start+i for i in exon_start_pos]
    exon_end=[exon_start[i]+exon_len[i]-1 for i in range(len(exon_len))]
    exon_start=[str(i) for i in exon_start]
    exon_end=[str(i) for i in exon_end]
    isoID=x[3]
    strand=x[5]
    arr=[chr,start,end,isoID,strand,','.join(exon_start),','.join(exon_end),sum(exon_len)]
    return(arr)
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
    
def comparePartExon2gene(geneExon,exon_start,exon_end):
    start_type,end_type,exon_score,inexon_score=compareExon2gene(geneExon,exon_start,exon_end)
    ref_start=[i[0]+1 for i in geneExon]
    ref_end=[i[1] for i in geneExon]
    gene_start=ref_start[0]
    gene_end=ref_end[-1]
    num=len(exon_start)
    overlapIdx_start=[]
    overlapIdx_end=[]
    for i in range(num):
        if (exon_start[i]>=gene_start) and (exon_start[i]<=gene_end):
            overlapIdx_start.append(1)
        else:
            overlapIdx_start.append(0)
        if (exon_end[i]>=gene_start) and (exon_end[i]<=gene_end):
            overlapIdx_end.append(1)
        else:
            overlapIdx_end.append(0)
    start_type=start_type.split(',')
    end_type=end_type.split(',')
    for i in range(num):
        if overlapIdx_start[i]==0:
            start_type[i]='intergenic'
        if overlapIdx_end[i]==0:
            end_type[i]='intergenic'
    start_type=','.join(start_type)
    end_type=','.join(end_type)
    return(start_type,end_type,exon_score,inexon_score)
def locatePos(geneExon,pos):
    ref_start=[i[0]+1 for i in geneExon]
    ref_end=[i[1] for i in geneExon]
    ref_num=len(ref_start)
    judgeBody='intron'
    for i in range(ref_num):
        if (ref_start[i]<pos) and (ref_end[i]>pos):
            judgeBody='inE'
    return(judgeBody)
def judgeIntron(geneExon,s,e):
    ref_start=[i[0]+1 for i in geneExon]
    ref_end=[i[1] for i in geneExon]
    ref_num=len(ref_start)
    for i in range(ref_num):
        if s<ref_start[i] and e>ref_end[i]:
            return(True)
    return(False)
    
def compareExon2gene(geneExon,exon_start,exon_end):
    ref_start=[i[0]+1 for i in geneExon]
    ref_end=[i[1] for i in geneExon]
    num=len(exon_start)
    start_match=[i in ref_start for i in exon_start]
    end_match=[i in ref_end for i in exon_end]
    in_start_match=0
    in_end_match=0
    for i in range(num):
        for j in range(len(ref_start)):
            if (exon_start[i]>=ref_start[j]) & (exon_start[i]<=ref_end[j]):
                in_start_match+=1
            if (exon_end[i]>=ref_start[j]) & (exon_end[i]<=ref_end[j]):
                 in_end_match+=1
    
    start_pos=[]
    end_pos=[]
    for i in range(num):
        if start_match[i]:
            start_pos.append('m')
        else:
            start_pos.append(locatePos(geneExon,exon_start[i]))
        if end_match[i]:
            end_pos.append('m')
        else:
            end_pos.append(locatePos(geneExon,exon_end[i]))
        if start_pos[-1]=='intron':
            if end_pos[-1] !='intron':
                start_pos[-1]='outE'
            else:
                if judgeIntron(geneExon,exon_start[i],exon_end[i]):
                    start_pos[-1]='outE'
                    end_pos[-1]='outE'
        if end_pos[-1]=='intron':
            if start_pos[-1] !='intron':
                end_pos[-1]='outE'
    match_score=sum(start_match)+sum(end_match)        
    start_pos_all=','.join(start_pos)
    end_pos_all=','.join(end_pos)
    return(start_pos_all,end_pos_all,match_score,in_start_match+in_end_match)  

def getCircType(tmp):
    exon_start=[int(i) for i in tmp['exon_start'].split(',')]
    exon_end=[int(i) for i in tmp['exon_end'].split(',')]

    canGene_full=[] # all in the gene with same sense
    canGene_part=[] # only a part in the gene with same sense
    canGene_anti=[] # antisense
    for gtf in tabixfile.fetch(tmp['chr'], tmp['start']-1, tmp['end']-1,parser=pysam.asGTF()):
        if gtf.feature=='gene':
            if gtf.strand==tmp['strand']:
                if (gtf.start<=tmp['start']) and (gtf.end>=tmp['end']):
                    canGene_full.append(gtf.gene_id)
                else:
                    canGene_part.append(gtf.gene_id)
            else:
                canGene_anti.append(gtf.gene_id)                
    start_type=''
    end_type=''
    circ_type=''
    geneName=''
    if len(canGene_full)>0:
        circ_type='full'
        exon_score=-1
        inexon_score=-1
        geneName=''
        for i in canGene_full:
            geneExon=gene2exon_dict[i]
            start_type_tmp,end_type_tmp,exon_score_tmp,inexon_score_tmp=compareExon2gene(geneExon,exon_start,exon_end)
            if exon_score_tmp>exon_score:
                start_type=start_type_tmp
                end_type=end_type_tmp
                exon_score=exon_score_tmp
                inexon_score=inexon_score_tmp
                geneName=gene2class_dict[i]
            elif exon_score_tmp==exon_score:
                if inexon_score_tmp>inexon_score:
                    start_type=start_type_tmp
                    end_type=end_type_tmp
                    inexon_score=inexon_score_tmp
                    geneName=gene2class_dict[i]
            else:
                continue
    elif len(canGene_part)>0:
        if len(canGene_part)>1:
            min_start,max_end,geneName,geneID_minmax=getMinMax(canGene_part)
            if (min_start<tmp['start']) and (max_end>=tmp['end']):
                circ_type='read through'
                geneExon_1=gene2exon_dict[geneID_minmax[0]]
                start_type_tmp_1,end_type_tmp_1,exon_score_tmp_1,inexon_score_tmp1=comparePartExon2gene(geneExon_1,exon_start,exon_end)   
                geneExon_2=gene2exon_dict[geneID_minmax[1]]
                start_type_tmp_2,end_type_tmp_2,exon_score_tmp_2,inexon_score_tmp2=comparePartExon2gene(geneExon_2,exon_start,exon_end)   
                return(combin2type(start_type_tmp_1,start_type_tmp_2),combin2type(end_type_tmp_1,end_type_tmp_2),circ_type,geneName)
        circ_type='part'
        exon_score=-1
        inexon_score=-1
        geneName=''
        for i in canGene_part:
            geneExon=gene2exon_dict[i]
            start_type_tmp,end_type_tmp,exon_score_tmp,inexon_score_tmp=comparePartExon2gene(geneExon,exon_start,exon_end)
            if exon_score_tmp>exon_score:
                start_type=start_type_tmp
                end_type=end_type_tmp
                exon_score=exon_score_tmp
                inexon_score=inexon_score_tmp
                geneName=gene2class_dict[i]
            elif exon_score_tmp==exon_score:
                if inexon_score_tmp>inexon_score:
                    start_type=start_type_tmp
                    end_type=end_type_tmp
                    inexon_score=inexon_score_tmp
                    geneName=gene2class_dict[i]
            else:
                continue  
    elif len(canGene_anti)>0:
        circ_type='antisense'
        geneName=';'.join([gene2class_dict[i]+'-AS1' for i in canGene_anti])
        start_type=','.join(['a']* len(exon_start))
        end_type=','.join(['a']* len(exon_end))
    else:
        circ_type='intergenic'
        start_type=','.join(['intergenic']* len(exon_start))
        end_type=','.join(['intergenic']* len(exon_end))
    return(start_type,end_type,circ_type,geneName)

def combin2type(t1,t2):
    t1=t1.split(',')
    t2=t2.split(',')
    num=len(t1)
    tc=[]
    for i in range(num):
        if t1[i]!='intergenic' or t2[i]=='intergenic':
            tc.append(t1[i])
        elif t1[i]=='intergenic' or t2[i]!='intergenic':
            tc.append(t2[i])
        elif t1[i]=='intergenic' or t2[i]=='intergenic':
            tc.append('intergenic')
        else:
            if 'm' in [t1[i],t2[i]]:
                tc.append('m')
            else:
                tc.append(t1[i])
    return(','.join(tc))
    
def getMinMax(canGene_part):
    min_start=float('Inf')
    max_end=0
    geneName_min=''
    geneName_max=''
    geneID_min=''
    geneID_max=''
    for i in canGene_part:
        if gene2status_dict[i] in ['NOVEL']: # remove sense_overlapping misc_RNA
            continue
        fail=True
        for j in gene2trans_dict[i]:
            if trans2class_dict[j] not in ['NOVEL']:
                fail=False
                break
        if fail:
            continue
        geneExon=gene2exon_dict[i]
        min_start_tmp=min([i[0] for i in geneExon])
        max_end_tmp=max([i[1] for i in geneExon])
        if min_start_tmp< min_start:
            min_start=min_start_tmp
            geneName_min=gene2class_dict[i]
            geneID_min=i
        if max_end_tmp>max_end:
            max_end=max_end_tmp
            geneName_max=gene2class_dict[i]
            geneID_max=i
            
    return(min_start,max_end,geneName_min+';'+geneName_max,[geneID_min,geneID_max])

def bed2df(bedFile,gtfFile):
    global tabixfile,gene2strand_dict,gene2trans_dict,gene2exon_dict,trans2exon_dict,trans2gene_dict,gene2status_dict,gene2class_dict,trans2class_dict
    FL_bed=pd.read_csv(bedFile,sep='\t',header=None)
    FL_bed.iloc[:,10]=FL_bed.iloc[:,10].map(str)
    FL_bed.iloc[:,11]=FL_bed.iloc[:,11].map(str)
    FL=pd.DataFrame(FL_bed.apply(lambda x: bed2FL(x),axis=1).tolist(),columns=['chr','start','end','isoID','strand','exon_start','exon_end','len'])
    tabixfile=pysam.TabixFile(gtfFile)
    gene2strand_dict={}
    gene2trans_dict={}
    gene2exon_dict={}
    trans2exon_dict={}
    trans2gene_dict={}
    gene2status_dict={}
    gene2class_dict={}
    trans2class_dict={}
    
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
            if 'gene_status' in gtf.keys():
                gene2status_dict[current_geneID]=gtf.gene_status
            else:
                gene2status_dict[current_geneID]='NOVEL'
            gene2class_dict[current_geneID]=gtf.gene_name
    for i in gene2exon_dict.keys():
        gene2exon_dict[i]=noDup_list(gene2exon_dict[i])
        
    start_type_list=[]
    end_type_list=[]
    circ_type_list=[]
    geneName_list=[]
    for i in range(FL.shape[0]):
        tmp=FL.iloc[i,:]
        start_type,end_type,circ_type,geneName=getCircType(tmp)
        start_type_list.append(start_type)
        end_type_list.append(end_type)
        circ_type_list.append(circ_type)
        geneName_list.append(geneName)
        
    FL['start_type']=start_type_list
    FL['end_type']=end_type_list
    FL['geneName']=geneName_list
    
    BSJ_type_list=[]
    for i in range(FL.shape[0]):
        circ_type=circ_type_list[i]
        start_type=start_type_list[i].split(',')
        end_type=end_type_list[i].split(',')
        if circ_type in ['intergenic','antisense','read through']:
            BSJ_type=circ_type
        elif circ_type=='part':
            if start_type[0] =='intergenic' and end_type[-1]=='intergenic':
                BSJ_type='novel UTR5;3'
            elif start_type[0] =='intergenic':
                if FL.iloc[i,:]['strand']=='+':
                    BSJ_type='novel UTR5'
                else:
                    BSJ_type='novel UTR3'
            else:
                if FL.iloc[i,:]['strand']=='+':
                    BSJ_type='novel UTR3'
                else:
                    BSJ_type='novel UTR5'
        else:
            num=len(start_type)
            # judge BSJ type
            if start_type[0]=='m':
                if end_type[-1]=='m':
                    BSJ_type='m'
                elif end_type[-1] in ['inE','outE']:
                    if FL.iloc[i,:]['strand']=='+':
                        BSJ_type='AS5'
                    else:
                        BSJ_type='AS3'
                elif end_type[-1]=='intron':
                    BSJ_type='Intron retention'
                else:
                    print('1')
            elif start_type[0] in ['inE','outE']:
                if end_type[-1]=='m':
                    if FL.iloc[i,:]['strand']=='+':
                        BSJ_type='AS5'
                    else:
                        BSJ_type='AS3'
                elif end_type[-1] in ['inE','outE']:
                    BSJ_type='AS5,AS3'
                elif end_type[-1]=='intron':
                    BSJ_type='Intron retention'
                else:
                     print('1')
            elif start_type[0] =='intron':
                BSJ_type='Intron retention'
            else:
                 print('1')
                            
        BSJ_type_list.append(BSJ_type)
        
    FSJ_type_list=[]
    for i in range(FL.shape[0]):
        FSJ_type=[]
        circ_type=circ_type_list[i]
        start_type=start_type_list[i].split(',')
        end_type=end_type_list[i].split(',')
        if circ_type in ['intergenic','antisense','read through']:       
            FSJ_type=circ_type
        else:
            num=len(start_type)
        #judge FSJ type
            if num==1:
                FSJ_type='1'
            else:
                for j in range(num-1):
                    donar=end_type[j]
                    acceptor=start_type[j+1]
                    if donar=='m':
                        if acceptor=='m':
                            FSJ_type.append('m')
                        elif acceptor in ['inE','outE']:
                            if FL.iloc[i,:]['strand']=='+':
                                FSJ_type.append('AS3')
                            else:
                                FSJ_type.append('AS5')
                        elif acceptor=='intron':
                            FSJ_type.append('Intron retention')
                        else:
                            if FL.iloc[i,:]['strand']=='+':
                                FSJ_type.append('novel UTR3')
                            else:
                                FSJ_type.append('novel UTR5')
                    elif donar  in ['inE','outE']:
                        if acceptor=='m':
                            if FL.iloc[i,:]['strand']=='+':
                                FSJ_type.append('AS5')
                            else:
                                FSJ_type.append('AS3')
                        elif acceptor in ['inE','outE']:
                            FSJ_type.append('AS5,AS3')
                        elif acceptor=='intron':
                            FSJ_type.append('Intron retention')
                        else:
                            if FL.iloc[i,:]['strand']=='+':
                                FSJ_type.append('novel UTR3')
                            else:
                                FSJ_type.append('novel UTR5')
                    elif donar=='intron':
                        if acceptor=='intergenic':
                            if FL.iloc[i,:]['strand']=='+':
                                FSJ_type.append('novel UTR3')
                            else:
                                FSJ_type.append('novel UTR5')
                        else:
                            FSJ_type.append('Intron retention')
                    else:
                        if acceptor=='intergenic':
                            FSJ_type.append('novel UTR5;3')
                        else:
                            if FL.iloc[i,:]['strand']=='+':
                                FSJ_type.append('novel UTR5')
                            else:
                                FSJ_type.append('novel UTR3')
                FSJ_type=','.join(FSJ_type)          
        FSJ_type_list.append(FSJ_type)
        
    detail_BSJ_type=[]
    com_BSJ_type=[]
    for i in BSJ_type_list:
        if i=='AS3':
            n='N3SS'
            m='NSS'
        elif i=='AS5':
            n='N5SS'
            m='NSS'
        elif i=='AS5,AS3':
            n='N5SS,N3SS'
            m='NSS'
        elif i=='Intron retention':
            n='intronic'
            m='intronic'
        elif i=='antisense':
            n='antisense'
            m='antisense'
        elif i=='intergenic':
            n='intergenic'
            m='intergenic'
        elif i=='m':
            n='exonic'
            m='exonic'
        elif i=='novel UTR3':
            n='novel UTR3'
            m='novel UTR'
        elif i=='novel UTR5':
            n='novel UTR5'
            m='novel UTR'
        elif i=='novel UTR5;3':
            n='novel UTR5,UTR3'
            m='novel UTR'
        elif i=='read through':
            n='read through'
            m='read through'
        else:
            n='unknown'
            m='unknown'
        detail_BSJ_type.append(n)
        com_BSJ_type.append(m)
    FL['detail_type']=detail_BSJ_type
    FL['type']=com_BSJ_type
    return(FL)
    
def annoFL(options):
    bed=options['-b']
    anno=options['-a']
    outFile=options['-o']
    plog('Check anno file')
    readGTFfile(anno)
    plog('Check BED file')
    fileCheck(anno)
    plog('Begin annotation')
    FL=bed2df(bed,anno)
    FL.to_csv(outFile,header=True,index=None,sep='\t')
    plog('Done!!!')