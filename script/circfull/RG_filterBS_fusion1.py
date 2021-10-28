import pandas as pd, numpy as np,pyfasta, sys,warnings,os
warnings.filterwarnings("ignore")
errorLen=40

def mapExon(x):
    chr1=x['leftSchr']
    chr2=x['leftEchr']
    leftSpos=x['leftSpos'].split(',')
    leftEpos=x['leftEpos'].split(',')
    rightSpos=x['rightSpos'].split(',')
    rightEpos=x['rightEpos'].split(',')
    leftSseq=x['leftSseq'].split(',')
    leftEseq=x['leftEseq'].split(',')
    rightSseq=x['rightSseq'].split(',')
    rightEseq=x['rightEseq'].split(',')
    strand=x['strand']
    isAdj=[]
    if strand:
        for i in range(len(leftSpos)):
            
            leftSposKey=chr1+'_'+leftSpos[i]
            if ExonEdict.__contains__(leftSposKey):
                tmpExonEdict=ExonEdict[leftSposKey]
                tmpleftSpos=int(leftSpos[i])
                tmpDiff=abs(np.array(tmpExonEdict)-tmpleftSpos)
                leftSpos[i]=int(tmpExonEdict[np.where(np.array(tmpDiff)==min(tmpDiff))[0][0]])
                leftSseq[i]=genome.sequence({'chr':chr1,'start':leftSpos[i]+1,'stop':leftSpos[i]+2})
                isAdj.append('leftS')
            leftEposKey=chr2+'_'+leftEpos[i]
            if ExonSdict.__contains__(leftEposKey):
                tmpExonSdict=ExonSdict[leftEposKey]
                tmpleftEpos=int(leftEpos[i])
                tmpDiff=abs(np.array(tmpExonSdict)-tmpleftEpos)
                leftEpos[i]=tmpExonSdict[np.where(np.array(tmpDiff)==min(tmpDiff))[0][0]]
                leftEseq[i]=genome.sequence({'chr':chr2,'start':leftEpos[i]-2,'stop':leftEpos[i]-1})
                isAdj.append('leftE')
        for i in range(len(rightSpos)):
            rightSposKey=chr2+'_'+rightSpos[i]
            if ExonEdict.__contains__(rightSposKey):
                tmpExonEdict=ExonEdict[rightSposKey]
                tmprightSpos=int(rightSpos[i])
                tmpDiff=abs(np.array(tmpExonEdict)-tmprightSpos)
                rightSpos[i]=tmpExonEdict[np.where(np.array(tmpDiff)==min(tmpDiff))[0][0]]
                rightSseq[i]=genome.sequence({'chr':chr2,'start':rightSpos[i]+1,'stop':rightSpos[i]+2})
                isAdj.append('rightS')
            rightEposKey=chr1+'_'+rightEpos[i]
            if ExonSdict.__contains__(rightEposKey):
                tmpExonSdict=ExonSdict[rightEposKey]
                tmprightEpos=int(rightEpos[i])
                tmpDiff=abs(np.array(tmpExonSdict)-tmprightEpos)
                rightEpos[i]=tmpExonSdict[np.where(np.array(tmpDiff)==min(tmpDiff))[0][0]]
                rightEseq[i]=genome.sequence({'chr':chr1,'start':rightEpos[i]-2,'stop':rightEpos[i]-1})
                isAdj.append('rightE')
    else:
        for i in range(len(leftSpos)):
            leftSposKey=chr1+'_'+leftSpos[i]
            if ExonEdict.__contains__(leftSposKey):
                tmpExonEdict=ExonEdict[leftSposKey]
                tmpleftSpos=int(leftSpos[i])
                tmpDiff=abs(np.array(tmpExonEdict)-tmpleftSpos)
                leftSpos[i]=tmpExonEdict[np.where(np.array(tmpDiff)==min(tmpDiff))[0][0]]
                leftSseq[i]=genome.sequence({'chr':chr1,'start':leftSpos[i]+1,'stop':leftSpos[i]+2})
                isAdj.append('leftS')
            leftEposKey=chr2+'_'+leftEpos[i]
            if ExonEdict.__contains__(leftEposKey):
                tmpExonEdict=ExonEdict[leftEposKey]
                tmpleftEpos=int(leftEpos[i])
                tmpDiff=abs(np.array(tmpExonEdict)-tmpleftEpos)
                leftEpos[i]=tmpExonEdict[np.where(np.array(tmpDiff)==min(tmpDiff))[0][0]]
                leftEseq[i]=genome.sequence({'chr':chr2,'start':leftEpos[i]+1,'stop':leftEpos[i]+2,'strand':'-'})
                isAdj.append('leftE')
        for i in range(len(rightSpos)):
            rightSposKey=chr2+'_'+rightSpos[i]
            if ExonSdict.__contains__(rightSposKey):
                tmpExonSdict=ExonSdict[rightSposKey]
                tmprightSpos=int(rightSpos[i])
                tmpDiff=abs(np.array(tmpExonSdict)-tmprightSpos)
                rightSpos[i]=tmpExonSdict[np.where(np.array(tmpDiff)==min(tmpDiff))[0][0]]
                rightSseq[i]=genome.sequence({'chr':chr2,'start':rightSpos[i]-2,'stop':rightSpos[i]-1,'strand':'-'})
                isAdj.append('rightS')
            rightEposKey=chr1+'_'+rightEpos[i]
            if ExonSdict.__contains__(rightEposKey):
                tmpExonSdict=ExonSdict[rightEposKey]
                tmprightEpos=int(rightEpos[i])
                tmpDiff=abs(np.array(tmpExonSdict)-tmprightEpos)
                rightEpos[i]=tmpExonSdict[np.where(np.array(tmpDiff)==min(tmpDiff))[0][0]]
                rightEseq[i]=genome.sequence({'chr':chr1,'start':rightEpos[i]-2,'stop':rightEpos[i]-1})
                isAdj.append('rightE')
    x['leftSpos']=','.join(map(str,leftSpos))
    x['leftEpos']=','.join(map(str,leftEpos))
    x['rightSpos']=','.join(map(str,rightSpos))
    x['rightEpos']=','.join(map(str,rightEpos))
    x['leftSseq']=','.join(leftSseq)
    x['leftEseq']=','.join(leftEseq)
    x['rightSseq']=','.join(rightSseq)
    x['rightEseq']=','.join(rightEseq)
    return([x,len(set(isAdj))])
        
def adjustBS(tmp,type='type1'):
    tmp=tmp.copy()
    chr1=tmp['chr1']
    chr2=tmp['chr2']
    strand=tmp['strand']
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
            if len(set(tmpSBoth_table) & set(['GTAG','CTAC']))>0:
                for i in range(len(tmpSBoth_table)):
                    if tmpSBoth_table[i] in ['GTAG','CTAC']:
                        idx_potential_tmpPBoth.append(i)
                potential_tmpPlen=tmpPlen_table[idx_potential_tmpPBoth]
                potential_tmpPlen_max=max(potential_tmpPlen)
                idx_potential_tmpPBoth_only=tmpPBoth_table_can[idx_potential_tmpPBoth[np.where(potential_tmpPlen==potential_tmpPlen_max)[0][0]]] 
            elif len(set(tmpSS_table) & set(['GT','CT']))>0:
                for i in range(len(tmpSS_table)):
                    if tmpSS_table[i] in ['GT','CT']:
                        idx_potential_tmpPBoth.append(i)
                potential_tmpPlen=tmpPlen_table[idx_potential_tmpPBoth]
                potential_tmpPlen_max=max(potential_tmpPlen)
                idx_potential_tmpPBoth_only=tmpPBoth_table_can[idx_potential_tmpPBoth[np.where(potential_tmpPlen==potential_tmpPlen_max)[0][0]]]
            elif len(set(tmpSE_table) & set(['AG','AC']))>0:
                for i in range(len(tmpSE_table)):
                    if tmpSE_table[i] in ['AG','AC']:
                        idx_potential_tmpPBoth.append(i)
                potential_tmpPlen=tmpPlen_table[idx_potential_tmpPBoth]
                potential_tmpPlen_max=max(potential_tmpPlen)
                idx_potential_tmpPBoth_only=tmpPBoth_table_can[idx_potential_tmpPBoth[np.where(potential_tmpPlen==potential_tmpPlen_max)[0][0]]]
            else:
                potential_tmpPlen_max=max(tmpPlen_table)
                idx_potential_tmpPBoth_only=tmpPBoth_table_can[np.where(tmpPlen_table==potential_tmpPlen_max)[0][0]]
            tmpPBothDF_potentail=tmpPBothDF.loc[idx_potential_tmpPBoth_only]
            if len(tmpPBothDF_potentail.shape) == 1:
                pos_start=tmpPBothDF_potentail['start']
                pos_end=tmpPBothDF_potentail['end']
            else:
                pos_start=tmpPBothDF_potentail['start'].iloc[0]
                pos_end=tmpPBothDF_potentail['end'].iloc[0]
    if strand:
        seq_start=genome.sequence({'chr': chr1, 'start':pos_start+1, 'stop':pos_start+2}).upper()
        seq_end=genome.sequence({'chr': chr2, 'start':pos_end-2, 'stop':pos_end-1}).upper()
    else:
        if type == 'type1':
            seq_start=genome.sequence({'chr': chr1, 'start':pos_start+1, 'stop':pos_start+2}).upper()
            seq_end=genome.sequence({'chr': chr2, 'start':pos_end+1, 'stop':pos_end+2,'strand':'-'}).upper()
        else:
            seq_start=genome.sequence({'chr': chr1, 'start':pos_start-2, 'stop':pos_start-1,'strand':'-'}).upper()
            seq_end=genome.sequence({'chr': chr2, 'start':pos_end-2, 'stop':pos_end-1}).upper()
    if type == 'type1':
        return('\t'.join([tmp['ID'],chr1,chr2,str(pos_start),str(pos_end),seq_start,seq_end]))
    else:
        return('\t'.join([str(pos_start),str(pos_end),seq_start,seq_end,str(strand)]))
        

def filterBS_fusion1(genomeFile,gtfFile,outPrefix):
    global genome,ExonSdict,ExonEdict
    genome = pyfasta.Fasta(genomeFile)
    BSFile=outPrefix+"BS_Fusion1.txt"
    BS_Fusion1=pd.read_csv(BSFile,names=['ID','leftSchr','leftEchr','leftSpos','leftEpos','leftSseq','leftEseq',
    'rightSchr','rightEchr','rightSpos','rightEpos','rightSseq','rightEseq','strand'],sep='\t',dtype=object)
    BS_Fusion1=BS_Fusion1.dropna(axis = 0)
    if BS_Fusion1.shape[0]<1:
        return(False)
    
    ExonSdict={}
    ExonEdict={}

    gtf = pd.read_csv(gtfFile,sep="\t",names=['chr','source','type','start','end','score','strand','phase','attributes'],comment='#')
    gtf_exon=gtf[gtf['type']=='exon'].sort_values(by=['chr','start']).copy()

    gtf_exon['key']=gtf_exon['chr'].values+'|'+gtf_exon['start'].map(str).values+gtf_exon['strand'].values+gtf_exon['end'].map(str).values
    gtf_exon=gtf_exon.drop_duplicates('key')
    gtf_exon_start=gtf_exon.drop_duplicates('start')
    gtf_exon_end=gtf_exon.drop_duplicates('end')
    gtf_exon_start.index=range(gtf_exon_start.shape[0])
    gtf_exon_end.index=range(gtf_exon_end.shape[0])
    for i in range(gtf_exon_start.shape[0]):
        tmp=gtf_exon_start.iloc[i]
        tmpChr=tmp['chr']
        tmpStart=tmp['start']
        for j in range(tmpStart-errorLen,tmpStart+errorLen+1):
            eachKey=tmpChr+'_'+str(j)
            if ExonSdict.__contains__(eachKey):
                ExonSdict[eachKey].append(tmpStart)
            else:
                ExonSdict[eachKey]=[tmpStart]
    for i in range(gtf_exon_end.shape[0]):
        tmp=gtf_exon_end.iloc[i]
        tmpChr=tmp['chr']
        tmpEnd=tmp['end']
        for j in range(tmpEnd-errorLen,tmpEnd+errorLen+1):
            eachKey=tmpChr+'_'+str(j)
            if ExonEdict.__contains__(eachKey):
                ExonEdict[eachKey].append(tmpEnd)
            else:
                ExonEdict[eachKey]=[tmpEnd]
    
    BS_Fusion1=pd.DataFrame([mapExon(BS_Fusion1.iloc[i].copy())[0] for i in range(BS_Fusion1.shape[0])])
    BS_Fusion1_left=BS_Fusion1[['ID','leftSchr','leftEchr','leftSpos','leftEpos','leftSseq','leftEseq','strand']].copy()
    BS_Fusion1_right=BS_Fusion1[['ID','rightSchr','rightEchr','rightSpos','rightEpos','rightSseq','rightEseq','strand']].copy()
    BS_Fusion1_left.columns=['ID','chr1','chr2','start','end','leftSeq','rightSeq','strand']
    BS_Fusion1_right.columns=['ID','chr1','chr2','start','end','leftSeq','rightSeq','strand']
    BS_Fusion1_adj=pd.DataFrame([(adjustBS(BS_Fusion1_left.iloc[i])+'\t'+adjustBS(BS_Fusion1_right.iloc[i],type='type2')).split('\t') for i in range(BS_Fusion1_left.shape[0])],columns=['ID','leftChr','rightChr','leftSpos','leftEpos','leftSseq','leftEseq','rightSpos','rightEpos','rightSseq','rightEseq','strand'])
    BS_Fusion1_adj=BS_Fusion1_adj.sort_values(['leftChr','rightChr','leftSpos','leftEpos'])
    BS_Fusion1_adj['circID']=BS_Fusion1_adj['leftChr']+'|'+BS_Fusion1_adj['rightChr']+'|'+BS_Fusion1_adj['leftSpos']+'|'+BS_Fusion1_adj['leftEpos']+'|'+BS_Fusion1_adj['rightSpos']+'|'+BS_Fusion1_adj['rightEpos']
    BS_Fusion1_adj.to_csv(outPrefix+'BS_Fusion1_adj.txt',sep='\t',index=None)
    BS_Fusion1_adj_1=BS_Fusion1_adj.loc[:,['ID','leftChr','rightEpos','leftSpos','rightEseq','leftSseq']].copy()
    BS_Fusion1_adj_1.to_csv(outPrefix+'BS_Fusion1_adj_1.txt',sep='\t',index=None,header=None)
    leftList=[]
    rightList=[]
    leftListS=[]
    rightListS=[]
    for i in range(BS_Fusion1_adj.shape[0]):
        tmp=BS_Fusion1_adj.iloc[i]
        if tmp['strand'] == 'True':
            leftList.append(tmp['leftEpos'])
            rightList.append(tmp['rightSpos'])
            leftListS.append(tmp['leftEseq'])
            rightListS.append(tmp['rightSseq'])
        else:
            leftList.append(tmp['rightSpos'])
            rightList.append(tmp['leftEpos'])
            leftListS.append(tmp['rightSseq'])
            rightListS.append(tmp['leftEseq'])
    BS_Fusion1_adj_2=pd.DataFrame({'ID':list(BS_Fusion1_adj['ID'].values),'chr':list(BS_Fusion1_adj['rightChr'].values),'start':leftList,'end':rightList,'leftSeq':leftListS,'rightSeq':rightListS})
    BS_Fusion1_adj_2.to_csv(outPrefix+'BS_Fusion1_adj_2.txt',sep='\t',index=None,header=None)
    return(True)