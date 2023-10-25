import pysam,pandas as pd,sys,os,numpy as np
RTthre=1000000
def FL2bed(FL_noDup):
    FL_bed=FL_noDup.loc[:,['chr','start','end','isoID']]
    FL_bed['start']=FL_bed['start']-1
    FL_bed['score']=0
    FL_bed['strand']=FL_noDup['strand']
    FL_bed['thickStart']=0
    FL_bed['thickEnd']=0
    FL_bed['itemRgb']=0
    exonNum=FL_noDup['exon_start'].map(lambda x: len(str(x).split(','))).tolist()
    FL_bed['blockCount ']=exonNum
    exonSize=FL_noDup.apply(lambda x: getSize(x),axis=1).tolist()
    FL_bed['blockSizes']=exonSize
    startPos=FL_noDup.apply(lambda x: getStart(x),axis=1).tolist()
    FL_bed['blockStarts']=startPos
    FL_bed.loc[FL_bed['strand']=='U','strand']=''
    return(FL_bed)
    
def fus2toNormal(x):
    ID=x['ID']
    circID=x['chr_first']+'|'+str(x['start_first'])+'|'+str(x['end_second'])
    tmpChr=x['chr_first']
    start=x['start_first']
    end=x['end_second']
    tmpLen=x['len_first']+x['len_second']
    exonNum=x['exonNum_first']+x['exonNum_second']
    exons_start=str(x['exon_start_first'])+','+str(x['exon_start_second'])
    exons_end=str(x['exon_end_first'])+','+str(x['exon_end_second'])
    motif=x['exon_leftSeq_first'][0:2]+x['exon_rightSeq_second'][-2:]
    leftSeq=x['exon_leftSeq_first'][0:2]
    rightSeq=x['exon_rightSeq_second'][-2:]
    exon_leftSeq=x['exon_leftSeq_first']+','+x['exon_leftSeq_second']
    exon_rightSeq=x['exon_rightSeq_first']+','+x['exon_rightSeq_second']
    strand=x['strand_first']
    geneName=[]
    if type(x['geneName_first'])==str:
        geneName.append(x['geneName_first'])
    if type(x['geneName_second'])==str:
        geneName.append(x['geneName_second'])
    geneName=','.join(list(set(geneName)))
    return(ID,circID,tmpChr,start,end,tmpLen,exonNum,exons_start,exons_end,motif,leftSeq,rightSeq,exon_leftSeq,exon_rightSeq,strand,geneName)
def getSize(x):
    start=[int(i) for i in str(x['exon_start']).split(',')]
    end=[int(i) for i in str(x['exon_end']).split(',')]
    size=[]
    for i in range(len(start)):
        size.append(str(end[i]-start[i]+1))
    return(','.join(size))
def getStart(x):
    start=[int(i) for i in str(x['exon_start']).split(',')]
    s0=x['start']
    sarry=[]
    for i in range(len(start)):
        sarry.append(str(start[i]-s0))
    return(','.join(sarry))

def circFL_output(outPrefix):
    global FL_noDup

    FL=pd.read_csv(outPrefix+'result_Normal.txt',sep='\t')
    if os.path.exists(outPrefix+'fusion/result_fusion1.txt'):
        FL_fus1=pd.read_csv(outPrefix+'fusion/result_fusion1.txt',sep='\t')
    else:
        FL_fus1=pd.DataFrame(columns=['ID', 'circID', 'chr_first', 'start_first', 'end_first', 'len_first',
           'exonNum_first', 'exon_start_first', 'exon_end_first',
           'exon_leftSeq_first', 'exon_rightSeq_first', 'strand_first',
           'geneName_first', 'chr_second', 'start_second', 'end_second',
           'len_second', 'exonNum_second', 'exon_start_second', 'exon_end_second',
           'exon_leftSeq_second', 'exon_rightSeq_second', 'strand_second',
           'geneName_second'])
    if os.path.exists(outPrefix+'fusion/result_fusion2.txt'):
        FL_fus2=pd.read_csv(outPrefix+'fusion/result_fusion2.txt',sep='\t')
    else:
        FL_fus2=pd.DataFrame(columns=['ID', 'circID', 'chr_first', 'start_first', 'end_first', 'len_first',
           'exonNum_first', 'exon_start_first', 'exon_end_first',
           'exon_leftSeq_first', 'exon_rightSeq_first', 'strand_first',
           'geneName_first', 'chr_second', 'start_second', 'end_second',
           'len_second', 'exonNum_second', 'exon_start_second', 'exon_end_second',
           'exon_leftSeq_second', 'exon_rightSeq_second', 'strand_second',
           'geneName_second'])

    FL_fus2_same1=FL_fus2[(FL_fus2['geneName_first']==FL_fus2['geneName_second']) & (FL_fus2['strand_first']==FL_fus2['strand_second'])]
    FL_fus2_diff1=FL_fus2[((FL_fus2['geneName_first']==FL_fus2['geneName_second']) & (FL_fus2['strand_first']==FL_fus2['strand_second']))==False]
    
    FL_fus2_same2=FL_fus2_diff1[(FL_fus2_diff1['end_second']-FL_fus2_diff1['start_first'])<=RTthre]

    comID=FL_fus2_same1['ID'].tolist()+FL_fus2_same2['ID'].tolist()
    FL_fus2_same=FL_fus2[FL_fus2['ID'].isin(comID)]
    FL_fus2_diff=FL_fus2[[not i for i in FL_fus2['ID'].isin(FL_fus2_same['ID'])]]
    
    FL_fus1=pd.concat([FL_fus1,FL_fus2_diff],axis=0)
    if FL_fus1.shape[0]>0:
        FL_fus1['isoID']=FL_fus1['chr_first']+'|'+FL_fus1['exon_start_first'].map(str)+'|'+FL_fus1['exon_end_first'].map(str)+'|'+FL_fus1['chr_second']+'|'+FL_fus1['exon_start_second'].map(str)+'|'+FL_fus1['exon_end_second'].map(str)
        FL_fus1_isocount=FL_fus1.loc[:,'isoID'].value_counts()
        FL_fus1.index=FL_fus1['isoID']
        FL_fus1_noDup=FL_fus1.iloc[:,1:].drop_duplicates('isoID').copy()
        FL_fus1_noDup.loc[FL_fus1_isocount.index.tolist(),'readCount']=FL_fus1_isocount.values
        FL_fus1_noDup=FL_fus1_noDup.sort_values('readCount',ascending =False)
        readID_dict={}
        for i in range(FL_fus1.shape[0]):
            tmp=FL_fus1.iloc[i,:]
            if readID_dict.__contains__(tmp['isoID']):
                readID_dict[tmp['isoID']].append(tmp['ID'])
            else:
                readID_dict[tmp['isoID']]=[tmp['ID']]
        FL_fus1_noDup['readID']=FL_fus1_noDup['isoID'].map(lambda x: ','.join(readID_dict[x]))
        FL_fus1_noDup.loc[FL_fus1_noDup['geneName_first'].isnull(),'geneName_first']=''
        FL_fus1_noDup.loc[FL_fus1_noDup['geneName_second'].isnull(),'geneName_second']=''
        FL_fus1_noDup=FL_fus1_noDup.loc[:,['circID', 'isoID','chr_first', 'start_first', 'end_first', 'len_first',
           'exonNum_first', 'exon_start_first', 'exon_end_first',
           'exon_leftSeq_first', 'exon_rightSeq_first', 'strand_first',
           'geneName_first', 'chr_second', 'start_second', 'end_second',
           'len_second', 'exonNum_second', 'exon_start_second', 'exon_end_second',
           'exon_leftSeq_second', 'exon_rightSeq_second', 'strand_second',
           'geneName_second',  'readCount', 'readID']]
        FL_fus1_noDup.to_csv(outPrefix+'circFL_Fusion.txt',sep='\t',header=True,index=None)


        

    if FL_fus2_same.shape[0]>0:
        FL_fus2_same=pd.DataFrame(FL_fus2_same.apply(fus2toNormal,axis=1).tolist(),columns=FL.columns)
        FL=pd.concat([FL,FL_fus2_same])
    FL['exon_start']=FL['exon_start'].map(str).tolist()
    FL['exon_end']=FL['exon_end'].map(str).tolist()
    FL['isoID']=FL['chr']+'|'+FL['exon_start']+'|'+FL['exon_end']
    FL_isocount=FL.loc[:,'isoID'].value_counts()
    FL.index=FL['isoID']

    FL_noDup=FL.iloc[:,1:].drop_duplicates('isoID').copy()
    FL_noDup.loc[FL_isocount.index.tolist(),'readCount']=FL_isocount.values
    FL_noDup=FL_noDup.sort_values('readCount',ascending =False)

    readID_dict={}
    for i in range(FL.shape[0]):
        tmp=FL.iloc[i,:]
        if readID_dict.__contains__(tmp['isoID']):
            readID_dict[tmp['isoID']].append(tmp['ID'])
        else:
            readID_dict[tmp['isoID']]=[tmp['ID']]
            
    FL_noDup['readID']=FL_noDup['isoID'].map(lambda x: ','.join(readID_dict[x]))
    FL_noDup.loc[FL_noDup['geneName'].isnull(),'geneName']=''
    FL_noDup=FL_noDup.loc[:,['circID','isoID','chr','start','end','len','exonNum','exon_start','exon_end','motif','leftSeq','rightSeq','exon_leftSeq','exon_rightSeq','strand','geneName','readCount','readID']]
   
    FL_noDup.to_csv(outPrefix+'circFL_Normal.txt',sep='\t',header=True,index=None)
    FL_noDup['exon_start']=FL_noDup['exon_start'].map(str).tolist()
    FL_noDup['exon_end']=FL_noDup['exon_end'].map(str).tolist()
    FL_bed=FL2bed(FL_noDup)
    FL_bed.to_csv(outPrefix+'circFL_Normal.bed',sep='\t',header=None,index=None)