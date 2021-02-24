import pysam,pandas as pd,sys,os

def getSize(x):
    start=[int(i) for i in x['exon_start'].split(',')]
    end=[int(i) for i in x['exon_end'].split(',')]
    size=[]
    for i in range(len(start)):
        size.append(str(end[i]-start[i]+1))
    return(','.join(size))
def getStart(x):
    start=[int(i) for i in x['exon_start'].split(',')]
    s0=x['start']
    sarry=[]
    for i in range(len(start)):
        sarry.append(str(start[i]-s0))
    return(','.join(sarry))
def getStrand(x):
    score=scoreDict[x['ID']]
    motif=x['motif']
    exon_leftSeq=x['exon_leftSeq'].split(',')[1:]
    exon_rightSeq=x['exon_rightSeq'].split(',')[:-1]
    exon_motif=[]
    if len(exon_leftSeq)>0:
        for i in range(len(exon_leftSeq)):
            exon_motif.append(exon_leftSeq[i]+exon_rightSeq[i])
    if score>0:
        return('+')
    elif score<0:
        return('-')
    else:
        if motif in ['AGGT','AGGC']:
            return('+')
        elif motif in ['ACCT','GCCT']:
            return('-')
        else:
            if len(exon_motif)>0:
                for i in exon_motif:
                    if i in ['AGGT','AGGC']:
                        return('+')
                    elif i in ['ACCT','GCCT']:
                        return('-')
    return('U')

def annotateNormal(options):
    global scoreDict
    outPrefix=options[0]
    gtfFile=options[1]
    FL=pd.read_csv(outPrefix+'constructFL_Normal_adj.txt',sep='\t')
    FL=FL.dropna(axis=0)
    FL_th=pd.read_csv(outPrefix+'circSeq.th',sep='\t',names=['ID','consN','readLen','start','end','consLen','copyNum','full','consensus'])
    FL_th['usage']=(FL_th.end-FL_th.start)/FL_th.readLen
    FL_th=FL_th.loc[FL_th.usage>0.5]
    FL_th=FL_th.loc[FL_th.copyNum>1.2]
    FL.index=FL['ID']
    if len(list(set(FL.ID)-set(FL_th.ID)))==0:
        sys.exit('No circRNAs were detected!!!')
    FL=FL.loc[list(set(FL.ID)-set(FL_th.ID))]
    scoreDict={}
    if len(options)>2:
        strandFile=options[2]
        isSecond=True
        for i in range(strandFile.shape[0]):
            tmp=strandFile.iloc[i]
            scoreDict[tmp['ID']]=tmp['score']
    else:
        for i in list(FL.ID):
            scoreDict[i]=0
        
    strandList=FL.apply(getStrand,axis=1)

    FL['strand']=strandList
    FL_uniq=FL.drop_duplicates('circID').copy()
    FL_uniq=FL_uniq.set_index('circID')
    tabixfile = pysam.TabixFile(gtfFile)
    geneList=[]
    for i in range(FL_uniq.shape[0]):
        tmp=FL_uniq.iloc[i,]
        geneList1=[]
        geneList2=[]
        if tmp['strand'] == 'U':
            if tmp['chr']  in tabixfile.contigs:
                for gtf in tabixfile.fetch(tmp['chr'], tmp['start']-1, tmp['start'],parser=pysam.asGTF()):
                    if gtf.start==tmp['start']-1:
                        FL_uniq.iloc[i,-1]=gtf.strand
                        tmp=FL_uniq.iloc[i,]
                        #print("gtf.start=%d;FL=%d;strand=%s" % (gtf.start,tmp['start'],FL_uniq.iloc[i,-1]))
                        break
        if tmp['strand'] == 'U':
            if tmp['chr']  in tabixfile.contigs:
                for gtf in tabixfile.fetch(tmp['chr'], tmp['end']-1, tmp['end'],parser=pysam.asGTF()):
                    if gtf.end==tmp['end']-1:
                        FL_uniq.iloc[i,-1]=gtf.strand
                        tmp=FL_uniq.iloc[i,]
                        break
        if tmp['chr']  in tabixfile.contigs:
            for gtf in tabixfile.fetch(tmp['chr'], tmp['start']-1, tmp['start'],parser=pysam.asGTF()):
                if gtf.feature=='gene':
                    gn=gtf.gene_name
                    if gtf.strand==tmp['strand'] or tmp['strand'] == 'U':
                        geneList1.append(gn)
                    else:
                        if gn[-4:]=='-AS1':
                            geneList1.append(gn[0:(len(gn)-4)])
                        else:
                            geneList1.append(gn+'-AS1')
            for gtf in tabixfile.fetch(tmp['chr'], tmp['end']-1, tmp['end'],parser=pysam.asGTF()):
                if gtf.feature=='gene':
                    gn=gtf.gene_name
                    if gtf.strand==tmp['strand'] or tmp['strand'] == 'U':
                        geneList2.append(gn)
                    else:
                        if gn[-4:]=='-AS1':
                            geneList2.append(gn[0:(len(gn)-4)])
                        else:
                            geneList2.append(gn+'-AS1')
        geneCom=list(set(geneList1) & set(geneList2))
        if len(geneCom)>0:
            geneList.append(';'.join(geneCom))
        else:
            geneList1=[]
            if tmp['chr']  in tabixfile.contigs:
                for gtf in tabixfile.fetch(tmp['chr'], tmp['start']-1, tmp['end']-1,parser=pysam.asGTF()):
                    if gtf.feature=='gene':
                        gn=gtf.gene_name
                        if gtf.strand==tmp['strand'] or tmp['strand'] == 'U':
                            geneList1.append(gn)
                        else:
                            if gn[-4:]=='-AS1':
                                geneList1.append(gn[0:(len(gn)-4)])
                            else:
                                geneList1.append(gn+'-AS1')
            geneList.append(';'.join(list(set(geneList1))))
            
    FL_uniq['geneName']=geneList
    FL['geneName']=FL_uniq.loc[FL.circID,'geneName'].values
    FL['strand']=FL_uniq.loc[FL.circID,'strand'].values
    FL.to_csv(outPrefix+'result_Normal.txt',sep="\t",index=None)