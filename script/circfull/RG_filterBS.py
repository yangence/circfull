import pandas as pd, numpy as np,pyfasta, sys,os,warnings
from multiprocessing import Pool
warnings.filterwarnings("ignore")
errorLen=40

            
def mapExon(x):
    chr=x['chr']
    start=x['start'].split(',')
    end=x['end'].split(',')
    leftSeq=x['leftSeq'].split(',')
    rightSeq=x['rightSeq'].split(',')
    startKeyList=[chr+'_'+i for i in start]
    endKeyList=[chr+'_'+i for i in end]
    num=len(start)
    adjS=[]
    adjE=[]
    seqS=[]
    seqE=[]
    for i in range(num):
        startValue=[]
        endValue=[]
        startKey=startKeyList[i]
        endKey=endKeyList[i]
        if ExonSdict.__contains__(startKey):
            startValue=ExonSdict[startKey]
        if ExonEdict.__contains__(endKey):
            endValue=ExonEdict[endKey]
        cValueCan=list(set(startValue) & set(endValue))
        if isSecond and len(cValueCan)>0:
            strandScore=strandFile.loc[x['ID'],'score']
            cValue=[]
            for j in cValueCan:
                if (strandDict[j]=='+' and strandScore>=0) or (strandDict[j]=='-' and strandScore<=0):
                    cValue.append(j)
        else:
            cValue=cValueCan
        if len(cValue)==0:
            adjS.append(int(start[i]))
            adjE.append(int(end[i]))
        elif len(cValue)==1:
            tmp_gtf=gtf_exon.loc[cValue[0]]
            if len(tmp_gtf.shape)==1:
                adjS.append(tmp_gtf['start'])
                adjE.append(tmp_gtf['end'])
            else:
                tmp_gtf_start=np.array(tmp_gtf['start'])
                tmp_gtf_start_diff=abs(tmp_gtf_start-int(start[i]))
                adjS.append(tmp_gtf_start[np.where(tmp_gtf_start_diff==min(tmp_gtf_start_diff))[0][0]])            
                tmp_gtf_end=np.array(tmp_gtf['end'])
                tmp_gtf_end_diff=abs(tmp_gtf_end-int(end[i]))
                adjE.append(tmp_gtf_end[np.where(tmp_gtf_end_diff==min(tmp_gtf_end_diff))[0][0]])
        else:
            mSlist=[]
            mElist=[]
            mDifflist=[]
            for j in range(len(cValue)):
                tmp_gtf=gtf_exon.loc[cValue[j]]
                if len(tmp_gtf.shape)==1:
                    mSlist.append(tmp_gtf['start'])
                    mElist.append(tmp_gtf['end'])
                    mDifflist.append(abs(tmp_gtf['start']-int(start[i]))+abs(tmp_gtf['end']-int(end[i])))
                else:
                    tmp_gtf_start=np.array(tmp_gtf['start'])
                    tmp_gtf_start_diff=abs(tmp_gtf_start-int(start[i]))
                    mSlist.append(tmp_gtf_start[np.where(tmp_gtf_start_diff==min(tmp_gtf_start_diff))[0][0]])            
                    tmp_gtf_end=np.array(tmp_gtf['end'])
                    tmp_gtf_end_diff=abs(tmp_gtf_end-int(end[i]))
                    mElist.append(tmp_gtf_end[np.where(tmp_gtf_end_diff==min(tmp_gtf_end_diff))[0][0]])
                    mDifflist.append(min(tmp_gtf_start_diff)+min(tmp_gtf_end_diff))
            mDifflist=np.array(mDifflist)
            mIdx=np.where(mDifflist==min(mDifflist))[0][0]
            adjS.append(mSlist[mIdx])
            adjE.append(mElist[mIdx])
    for i in range(num):
        seqS.append(genome.sequence({'chr': chr, 'start':adjS[i]-2, 'stop':adjS[i]-1}).upper())
        seqE.append(genome.sequence({'chr': chr, 'start':adjE[i]+1, 'stop':adjE[i]+2}).upper())
        adjS[i]=str(adjS[i])
        adjE[i]=str(adjE[i])
    '''
    for i in range(num):
        motif=seqS[i]+seqE[i]
        if motif in ['ACCT','AGGT']:
            if leftSeq[i]+rightSeq[i] in ['ACCT','AGGT']:
                adjS[i]=start[i]
                adjE[i]=end[i]
                seqS[i]=leftSeq[i]
                seqE[i]=rightSeq[i]
    '''
    x['start']=','.join(adjS)
    x['end']=','.join(adjE)
    x['leftSeq']=','.join(seqS)
    x['rightSeq']=','.join(seqE)
    return(x)



def adjustBS(i):
    tmp=BS_Normal.iloc[i].copy()
    tmp=mapExon(tmp)
    chr=tmp['chr']
    tmpID=tmp['ID']
    tmpPS=np.array([int(i) for i in tmp['start'].split(',')])
    tmpPE=np.array([int(i) for i in tmp['end'].split(',')])
    tmpPlen=np.array(tmpPE)-np.array(tmpPS)
    tmpPlenMax=max(tmpPlen)
    tmpPBoth=[str(i[0])+'|'+str(i[1]) for i in zip(tmpPS,tmpPE)]
    tmpSS=tmp['leftSeq'].split(',')
    tmpSE=tmp['rightSeq'].split(',')
    tmpSBoth=[''.join(i) for i in zip(tmpSS,tmpSE)]
    tmpSBothDF=pd.DataFrame(zip(tmpSS,tmpSE,tmpSBoth,list(tmpPlen)),index=tmpPBoth,columns=['start','end','both','len'])
    tmpPBothDF=pd.DataFrame(zip(tmpPS,tmpPE),index=tmpPBoth,columns=['start','end'])    
    bsNum=len(tmpPS)
    baseRightMotif={'GT':'AG','CT':'AC'}
    baseLeftMotif={'AG':'GT','AC':'CT'}
    if True:
        if bsNum == 1:
            pos_start=tmpPS[0]
            adj_start = 'all_bsNum1'
            pos_end=tmpPE[0]
            adj_end = 'all_bsNum1'
        else:
            tmpPBoth_table = pd.Series(index=tmpPBoth).index.value_counts()
            tmpPBoth_table_max=max(tmpPBoth_table)
            tmpPBoth_table_can=list(tmpPBoth_table[tmpPBoth_table==tmpPBoth_table_max].index)
            tmpSBothDF=tmpSBothDF.loc[tmpPBoth_table_can]
            tmpPBoth_table_can=list(tmpSBothDF.index)
            tmpSBoth_table=list(tmpSBothDF['both'])
            tmpSS_table=list(tmpSBothDF['start'])
            tmpSE_table=list(tmpSBothDF['end'])
            tmpPlen_table=np.array(list(tmpSBothDF['len']))
            idx_potential_tmpPBoth=[]
            if isSecond:
                if len(set(tmpSBoth_table) & set(setMotif1[tmpID]))>0:
                    for i in range(len(tmpSBoth_table)):
                        if tmpSBoth_table[i] in ['AGGT','ACCT']:
                            idx_potential_tmpPBoth.append(i)
                    potential_tmpPlen=tmpPlen_table[idx_potential_tmpPBoth]
                    potential_tmpPlen_max=max(potential_tmpPlen)
                    idx_potential_tmpPBoth_only=tmpPBoth_table_can[idx_potential_tmpPBoth[np.where(potential_tmpPlen==potential_tmpPlen_max)[0][0]]] 
                    adj_start='all_bsNumM_BothMotif_BothMax'
                    adj_end='all_bsNumM_BothMotif_BothMax'
                elif len(set(tmpSS_table) & set(setMotif2[tmpID]))>0:
                    for i in range(len(tmpSS_table)):
                        if tmpSS_table[i] in ['AG','AC']:
                            idx_potential_tmpPBoth.append(i)
                    potential_tmpPlen=tmpPlen_table[idx_potential_tmpPBoth]
                    potential_tmpPlen_max=max(potential_tmpPlen)
                    idx_potential_tmpPBoth_only=tmpPBoth_table_can[idx_potential_tmpPBoth[np.where(potential_tmpPlen==potential_tmpPlen_max)[0][0]]]
                    adj_start='all_bsNumM_LeftMotif_BothMax'
                    adj_end='all_bsNumM_LeftMotif_RightnoMotif_BothMax'
                elif len(set(tmpSE_table) & set(setMotif3[tmpID]))>0:
                    for i in range(len(tmpSE_table)):
                        if tmpSE_table[i] in ['CT','GT']:
                            idx_potential_tmpPBoth.append(i)
                    potential_tmpPlen=tmpPlen_table[idx_potential_tmpPBoth]
                    potential_tmpPlen_max=max(potential_tmpPlen)
                    idx_potential_tmpPBoth_only=tmpPBoth_table_can[idx_potential_tmpPBoth[np.where(potential_tmpPlen==potential_tmpPlen_max)[0][0]]]
                    adj_start='all_bsNumM_LeftnoMotif_BothMax'
                    adj_end='all_bsNumM_LeftnoMotif_RightMotif_BothMax'
                else:
                    potential_tmpPlen_max=max(tmpPlen_table)
                    idx_potential_tmpPBoth_only=tmpPBoth_table_can[np.where(tmpPlen_table==potential_tmpPlen_max)[0][0]]
                    adj_start='all_bsNumM_noMotif_BothMax'
                    adj_end='all_bsNumM_noMotif_BothMax'
            else:
                if len(set(tmpSBoth_table) & set(['AGGT','ACCT']))>0:
                    for i in range(len(tmpSBoth_table)):
                        if tmpSBoth_table[i] in ['AGGT','ACCT']:
                            idx_potential_tmpPBoth.append(i)
                    potential_tmpPlen=tmpPlen_table[idx_potential_tmpPBoth]
                    potential_tmpPlen_max=max(potential_tmpPlen)
                    idx_potential_tmpPBoth_only=tmpPBoth_table_can[idx_potential_tmpPBoth[np.where(potential_tmpPlen==potential_tmpPlen_max)[0][0]]] 
                    adj_start='all_bsNumM_BothMotif_BothMax'
                    adj_end='all_bsNumM_BothMotif_BothMax'
                elif len(set(tmpSS_table) & set(['AG','AC']))>0:
                    for i in range(len(tmpSS_table)):
                        if tmpSS_table[i] in ['AG','AC']:
                            idx_potential_tmpPBoth.append(i)
                    potential_tmpPlen=tmpPlen_table[idx_potential_tmpPBoth]
                    potential_tmpPlen_max=max(potential_tmpPlen)
                    idx_potential_tmpPBoth_only=tmpPBoth_table_can[idx_potential_tmpPBoth[np.where(potential_tmpPlen==potential_tmpPlen_max)[0][0]]]
                    adj_start='all_bsNumM_LeftMotif_BothMax'
                    adj_end='all_bsNumM_LeftMotif_RightnoMotif_BothMax'
                elif len(set(tmpSE_table) & set(['CT','GT']))>0:
                    for i in range(len(tmpSE_table)):
                        if tmpSE_table[i] in ['CT','GT']:
                            idx_potential_tmpPBoth.append(i)
                    potential_tmpPlen=tmpPlen_table[idx_potential_tmpPBoth]
                    potential_tmpPlen_max=max(potential_tmpPlen)
                    idx_potential_tmpPBoth_only=tmpPBoth_table_can[idx_potential_tmpPBoth[np.where(potential_tmpPlen==potential_tmpPlen_max)[0][0]]]
                    adj_start='all_bsNumM_LeftnoMotif_BothMax'
                    adj_end='all_bsNumM_LeftnoMotif_RightMotif_BothMax'
                else:
                    potential_tmpPlen_max=max(tmpPlen_table)
                    idx_potential_tmpPBoth_only=tmpPBoth_table_can[np.where(tmpPlen_table==potential_tmpPlen_max)[0][0]]
                    adj_start='all_bsNumM_noMotif_BothMax'
                    adj_end='all_bsNumM_noMotif_BothMax'
            tmpPBothDF_potentail=tmpPBothDF.loc[idx_potential_tmpPBoth_only]
            if len(tmpPBothDF_potentail.shape) == 1:
                pos_start=tmpPBothDF_potentail['start']
                pos_end=tmpPBothDF_potentail['end']
            else:
                pos_start=tmpPBothDF_potentail['start'].iloc[0]
                pos_end=tmpPBothDF_potentail['end'].iloc[0]
    seq_start=genome.sequence({'chr': chr, 'start':pos_start-2, 'stop':pos_start-1}).upper()
    seq_end=genome.sequence({'chr': chr, 'start':pos_end+1, 'stop':pos_end+2}).upper()
    return('\t'.join([tmp['ID'],chr,str(pos_start),str(pos_end),seq_start,seq_end,str(adj_start),str(adj_end)]))

def filterBS(options):
    global genome,ExonSdict,ExonEdict,strandDict,BS_Normal,isSecond,gtf_exon,setMotif1,setMotif2,setMotif3,strandFile
    genomeFile=options[0]
    genome = pyfasta.Fasta(genomeFile)
    gtfFile=options[1]
    outPrefix=options[2]
    thread=options[3]
    isSecond=False

    BSFile=outPrefix+"BS_Normal.txt"
    fout=open(outPrefix+'BS_Normal_adj.txt','w')
    BS_Normal=pd.read_csv(BSFile,names=['ID','chr','start','end','leftSeq','rightSeq'],sep='\t')
    BS_Normal=BS_Normal.dropna(axis = 0)
    if len(options)>4:
        isSecond=True
        strandFile=options[4]
        strandFile.index=strandFile['ID']
        setMotif1={}
        setMotif2={}
        setMotif3={}
        for i in BS_Normal['ID'].values:
            tmpScore=strandFile.loc[i,'score']
            if tmpScore>0:
                setMotif1[i]=['AGGT']
                setMotif2[i]=['AG']
                setMotif3[i]=['GT']
            elif tmpScore<0:
                setMotif1[i]=['ACCT']
                setMotif2[i]=['AC']
                setMotif3[i]=['CT']
            else:
                setMotif1[i]=['AGGT','ACCT']
                setMotif2[i]=['AG','AC']
                setMotif3[i]=['GT','CT']

    gtf = pd.read_csv(gtfFile,sep="\t",names=['chr','source','type','start','end','score','strand','phase','attributes'],comment='#')
    gtf_exon=gtf[gtf['type']=='exon'].sort_values(by=['chr','start'])
    gtf_exon['key']=gtf_exon['chr'].values+'|'+gtf_exon['start'].map(str).values+gtf_exon['strand'].values+gtf_exon['end'].map(str).values
    gtf_exon['gene']=gtf_exon.apply(lambda x: x['attributes'].split(';')[0].split(' ')[1][1:-1],axis=1)
    gtf_exon['geneKey']=gtf_exon['gene']+'|'+gtf_exon['key']
    gtf_exon=gtf_exon.set_index('gene')
    gtf_exon=gtf_exon.sort_index()

    if os.path.exists(os.path.dirname(gtfFile)+'/ExonSdict.npy'):
        ExonSdict=np.load(os.path.dirname(gtfFile)+'/ExonSdict.npy',allow_pickle=True).item()
        ExonEdict=np.load(os.path.dirname(gtfFile)+'/ExonEdict.npy',allow_pickle=True).item()
        strandDict=np.load(os.path.dirname(gtfFile)+'/strandDict.npy',allow_pickle=True).item()
    else:
        if not isSecond and not (os.path.exists(outPrefix+'ExonSdict.npy') and os.path.exists(outPrefix+'ExonEdict.npy') and os.path.exists(outPrefix+'strandDict.npy')):
            ExonSdict={}
            ExonEdict={}
            strandDict={}
            for i in range(gtf_exon.shape[0]):
                tmp=gtf_exon.iloc[i]
                tmpChr=tmp['chr']
                tmpStart=tmp['start']
                tmpEnd=tmp['end']
                tmpValue=tmp.name
                strandDict[tmpValue]=tmp['strand']
                for j in range(errorLen*2+1):
                    startKey=tmpChr+'_'+str(tmpStart-errorLen+j)
                    endKey=tmpChr+'_'+str(tmpEnd-errorLen+j)
                    if ExonSdict.__contains__(startKey):
                        ExonSdict[startKey].append(tmpValue)
                    else:
                        ExonSdict[startKey]=[tmpValue]
                    if ExonEdict.__contains__(endKey):
                        ExonEdict[endKey].append(tmpValue)
                    else:
                        ExonEdict[endKey]=[tmpValue]
            np.save(outPrefix+'ExonSdict.npy',ExonSdict)
            np.save(outPrefix+'ExonEdict.npy',ExonEdict)
            np.save(outPrefix+'strandDict.npy',strandDict)
        else:
                if isSecond:
                    path=outPrefix+'../'
                else:
                    path=outPrefix
                
                ExonSdict=np.load(path+'ExonSdict.npy',allow_pickle=True).item()
                ExonEdict=np.load(path+'ExonEdict.npy',allow_pickle=True).item()
                strandDict=np.load(path+'strandDict.npy',allow_pickle=True).item()
            
    pool=Pool(processes=thread)
    result=pool.map(adjustBS,[i for i in range(BS_Normal.shape[0])])
    pool.close()
    pool.join()
    for i in result:
        fout.write(i+'\n')
    fout.close()
