import pandas as pd, numpy as np,pyfasta,sys

def estimateStrand(x):
    leftSeq=x['exon_leftSeq'].split(',')
    rightSeq=x['exon_rightSeq'].split(',')
    seqC=pd.Series(index=[leftSeq[i]+rightSeq[i] for i in range(len(leftSeq))]).index.value_counts()
    if seqC.index[0] in ['ACCT','AGGT']:
        if seqC.shape[0]>1:
            if (seqC.iloc[0]-seqC.iloc[1])>2:
                return(seqC.index[0])
        else:
            return(seqC.index[0])
    return('')   


# '+' means the same direction as RNA
def getStrand(x):
    motif=x['judgeStrand']
    if motif == '':
        return('')
    map=x['mapStrand']
    if motif == 'ACCT' and map == '+':
        return('-')
    elif motif == 'AGGT' and map == '+':
        return('+')
    elif motif == 'ACCT' and map == '-':
        return('+')
    else:
        return('-')

def judgeStrand(RG_out,strand_out):
    explainFL=pd.read_csv(RG_out+'explainFL_Normal_adj.txt',sep='\t')
    explainFL=explainFL.drop_duplicates('ID')
    explainFL=explainFL.sort_values('ID')
    explainFL.index=explainFL['ID']
    consFL=pd.read_csv(RG_out+"constructFL_Normal_adj.txt",sep='\t')
    consFL.index=consFL['ID']
    consFLHQ=consFL.iloc[[i  in ['ACCT','AGGT']  for i in consFL['motif']]].copy()
    consFLHQ_numM=consFLHQ[consFLHQ['exonNum']>1].copy()
    consFLHQ_num1=consFLHQ[consFLHQ['exonNum']==1].copy()
    consFLHQ_numM['judgeStrand']=consFLHQ_numM.apply(estimateStrand,axis=1)
    consFLHQ_numM['mapStrand']=explainFL.loc[consFLHQ_numM['ID'],'strand'].copy()
    consFLHQ_numM=consFLHQ_numM[consFLHQ_numM['judgeStrand']!=''].copy()

    consFLHQ_numM['seqStrand']=consFLHQ_numM.apply(getStrand,axis=1)

    consFLHQ_num1=consFLHQ_num1.sort_values('circID')
    consFLHQ_num1_table=consFLHQ_num1['circID'].value_counts()
    consFLHQ_num1_PASS_circID=consFLHQ_num1_table.iloc[np.where(consFLHQ_num1_table>3)].index
    if len(consFLHQ_num1_PASS_circID)==0:
        print('Warring: Expression is too low, which could influence strand predication.')
        consFLHQ_num1_PASS_circID=consFLHQ_num1_table.iloc[np.where(consFLHQ_num1_table>2)].index
        if len(consFLHQ_num1_PASS_circID)==0:
            consFLHQ_num1_PASS_circID=consFLHQ_num1_table.iloc[np.where(consFLHQ_num1_table>1)].index
            if len(consFLHQ_num1_PASS_circID)==0:
                consFLHQ_num1_PASS_circID=consFLHQ_num1_table.iloc[np.where(consFLHQ_num1_table>0)].index
    consFLHQ_num1_PASS=consFLHQ_num1.loc[[ i in consFLHQ_num1_PASS_circID for i in consFLHQ_num1['circID'] ]].copy()
    consFLHQ_num1_PASS.loc[consFLHQ_num1_PASS.index,'judgeStrand']=consFLHQ_num1_PASS['motif']
    consFLHQ_num1_PASS.loc[consFLHQ_num1_PASS.index,'mapStrand']=explainFL.loc[consFLHQ_num1_PASS['ID'],'strand'].copy()
    consFLHQ_num1_PASS.loc[consFLHQ_num1_PASS.index,'seqStrand']=consFLHQ_num1_PASS.apply(getStrand,axis=1)
    consFLHQ_PASS=pd.concat([consFLHQ_num1_PASS,consFLHQ_numM],axis=0)
    consFLHQ_PASS_uniq=consFLHQ_PASS.drop_duplicates('circID')

    leftDict={}
    rightDict={}

    for i in range(consFLHQ_PASS_uniq.shape[0]):
        tmp=consFLHQ_PASS_uniq.iloc[i,]
        chr=tmp['chr']
        start=int(tmp['start'])
        end=int(tmp['end'])
        read=tmp['ID']
        for j in range(41):
            leftKey=chr+'|'+str(start-20+j)
            if leftDict.__contains__(leftKey):
                leftDict[leftKey]+=[read]
            else:
                leftDict[leftKey]=[read]
            rightKey=chr+'|'+str(start-20+j)
            if rightDict.__contains__(rightKey):
                rightDict[rightKey]+=[read]
            else:
                rightDict[rightKey]=[read]

    consFLHQ_PASS_circID=list(consFLHQ_PASS_uniq['circID'])
    consFLHQ_PASS_circID_dict={}
    for i in consFLHQ_PASS_circID:
        consFLHQ_PASS_circID_dict[i]=1
        
    consFL_fail=consFL.loc[[not consFLHQ_PASS_circID_dict.__contains__(i) for i in consFL['circID']]].copy()
    consFL_fail_overlapID=[]
    for i in range(len(consFL_fail)):
        tmp=consFL_fail.iloc[i,]
        chr=tmp['chr']
        start=tmp['start']
        end=tmp['end']
        read=tmp['ID']
        leftKey=chr+'|'+str(start)
        rightKey=chr+'|'+str(end)
        if leftDict.__contains__(leftKey) and rightDict.__contains__(rightKey):
            overlapID=list(set(leftDict[leftKey]) and set(rightDict[rightKey]))
            if len(overlapID)>0:
                consFL_fail_overlapID.append(overlapID[0])
            else:
                consFL_fail_overlapID.append('')
        else:
            consFL_fail_overlapID.append('')

    consFL_fail.loc[consFL_fail.index,'judgeStrand']=''
    consFL_fail.loc[np.array(consFL_fail_overlapID)!='','judgeStrand']=list(consFLHQ_PASS_uniq.loc[np.array(consFL_fail_overlapID)[np.where(np.array(consFL_fail_overlapID)!='')],'judgeStrand'])
    consFL_fail.loc[consFL_fail.index,'mapStrand']=explainFL.loc[consFL_fail['ID'],'strand'].copy()
    consFL_fail.loc[consFL_fail.index,'seqStrand']=consFL_fail.apply(getStrand,axis=1)

    consFL_strand=pd.concat([consFLHQ_PASS,consFL_fail],axis=0,sort=False)
    consFL_strand=consFL_strand.set_index('ID')
    consFL_strand.to_csv(strand_out+'consFL_strand.txt',sep='\t')
