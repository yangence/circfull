import pandas as pd,sys,os
from .RG_circFL_output import FL2bed
def getR2I(dfa,dfb=None):
    read2iso={}
    if type(dfb)!=type(None):
        for i in range(dfb.shape[0]):
            tmp=dfb.iloc[i,]
            reads=tmp['readID'].split(',')
            isoID=tmp['isoID']
            for j in reads:
                read2iso[j]=isoID
    
    for i in range(dfa.shape[0]):
        tmp=dfa.iloc[i,]
        reads=tmp['readID'].split(',')
        isoID=tmp['isoID']
        for j in reads:
            read2iso[j]=isoID
    return(read2iso)
    
def getI2R(read2iso):
    iso2read={}
    for read in read2iso.keys():
        isoID=read2iso[read]
        if iso2read.__contains__(isoID):
            iso2read[isoID].append(read)
        else:
            iso2read[isoID]=[read]
    return(iso2read)
def getMdf(dfa,dfb=None,iso2read=None):
    if type(dfb)==type(None):
        ndfA=dfa
    else:
        ndfA=pd.concat([dfa,dfb])
    ndfA=ndfA.drop_duplicates('isoID')
    ndfA.index=ndfA.isoID
    ndfA=ndfA.loc[set(iso2read.keys()),:]
    allReadList=[]
    allReadCount=[]
    for i in range(ndfA.shape[0]):
        tmp=ndfA.iloc[i,:]
        isoID=tmp['isoID']
        readList=iso2read[isoID]
        allReadList.append(','.join(readList))
        allReadCount.append(len(readList))
    ndfA['readID']=allReadList
    ndfA['readCount']=allReadCount
    ndfA=ndfA.sort_values('readCount',ascending=False)
    return(ndfA)

def read2strand(df):
    r2s={}
    for i in range(df.shape[0]):
        tmp=df.iloc[i,]
        reads=tmp['readID'].split(',')
        strand=tmp['strand']
        for j in reads:
            r2s[j]=strand
    return(r2s)

def changeStrand(df,r2s,r2i,dfs):
    ndfA=pd.concat([df,dfs])
    ndfA=ndfA.drop_duplicates('isoID')
    ndfA.index=ndfA.isoID
    isoID_list=[]
    for i in range(df.shape[0]):
        tmp=df.iloc[i,]
        reads=tmp['readID'].split(',')
        strand=tmp['strand']
        strand_score=0
        isoID_dict={}
        for r in reads:
            if r2s.__contains__(r):
                if r2s[r]==strand:
                    strand_score+=1
                else:
                    strand_score+=-1
                    tmp_iso=r2i[r]
                    if isoID_dict.__contains__(tmp_iso):
                        isoID_dict[tmp_iso]+=1
                    else:
                        isoID_dict[tmp_iso]=1
        if strand_score>=0:
            isoID_list.append(tmp['isoID'])
        else:
            isoID=''
            tmpNum=0
            for j in isoID_dict.keys():
                if isoID_dict[j]>tmpNum:
                    isoID=j
                    tmpNum=isoID_dict[j]
            isoID_list.append(isoID)
    newDf=ndfA.loc[isoID_list,ndfA.columns[:-1]]
    newDf['readID']=df['readID'].tolist()
    n_r2i=getR2I(newDf)
    n_i2r=getI2R(n_r2i)
    newDf=getMdf(dfa=newDf,iso2read=n_i2r)
    return(newDf)
    
def merge(RG,outDir,cRG=False,sRG=False):
    ndf1=pd.read_csv(RG+'circFL_Normal.txt',sep='\t')
    if cRG:
        ndf2=pd.read_csv(cRG+'circFL_Normal.txt',sep='\t')
    if sRG:
        ndf3=pd.read_csv(sRG+'circFL_Normal.txt',sep='\t')
        ndf3_r2s=read2strand(ndf3)
        ndf3_r2i=getR2I(ndf3)
    if cRG:
        read2iso=getR2I(ndf1,ndf2)
        iso2read=getI2R(read2iso)
        ndf1_2=getMdf(ndf1,ndf2,iso2read)
        if sRG:
            ndf1_2_strand=ndf1_2.loc[ndf1_2.strand!='U',:]
            read2iso=getR2I(ndf1_2_strand,ndf3)
            iso2read=getI2R(read2iso)
            ndf1_2_3=getMdf(ndf1_2_strand,ndf3,iso2read)
            ndf1_2_3=changeStrand(ndf1_2_3,ndf3_r2s,ndf3_r2i,ndf3)
            ndf1_2_3.to_csv(outDir+'circFL_Normal.txt',sep='\t',index=None)
            ndf1_2_3_bed=FL2bed(ndf1_2_3)
            ndf1_2_3_bed.to_csv(outDir+'circFL_Normal.bed',sep='\t',header=None,index=None)
        else:
            ndf1_2.to_csv(outDir+'circFL_Normal.txt',sep='\t',index=None)
            ndf1_2_bed=FL2bed(ndf1_2)
            ndf1_2_bed.to_csv(outDir+'circFL_Normal.bed',sep='\t',header=None,index=None)
    else:
        ndf1_strand=ndf1.loc[ndf1.strand!='U',:]
        read2iso=getR2I(ndf1_strand,ndf3)
        iso2read=getI2R(read2iso)
        ndf1_3=getMdf(ndf1_strand,ndf3,iso2read)
        ndf1_3=changeStrand(ndf1_3,ndf3_r2s,ndf3_r2i,ndf3)
        ndf1_3.to_csv(outDir+'circFL_Normal.txt',sep='\t',index=None)
        ndf1_3_bed=FL2bed(ndf1_3)
        ndf1_3_bed.to_csv(outDir+'circFL_Normal.bed',sep='\t',header=None,index=None)
    '''
    fdf1=dir1+'circFL_Fusion.txt'
    fdf2=dir2+'circFL_Fusion.txt'
    if not os.path.exists(fdf1):
        print('ERROR: %s not exit!!!' % fdf1)
        return()
    if not os.path.exists(fdf2):
        print('ERROR: %s not exit!!!' % fdf2)
        return()
    '''