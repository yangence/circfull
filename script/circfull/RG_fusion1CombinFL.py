import pandas as pd,sys,pyfasta,sys,os,pysam


def getSeq(chr,exonS,exonE,strand):
    exonStart=[int(i) for i in exonS.split(',')]
    exonEnd=[int(i) for i in exonE.split(',')]
    seq=''
    for i in range(len(exonStart)):
        seq+=genome.sequence({'chr': chr, 'start':exonStart[i], 'stop':exonEnd[i]})
    if strand=='-':
        seq=seq.upper()[::-1].replace('A','t').replace('T','a').replace('C','g').replace('G','c').upper()
    return(seq)

def getStrand(x):
    exon_leftSeq_first=x['exon_leftSeq_first'].split(',')[1:]
    exon_rightSeq_first=x['exon_rightSeq_first'].split(',')[:-1]
    exon_leftSeq_second=x['exon_leftSeq_second'].split(',')[1:]
    exon_rightSeq_second=x['exon_rightSeq_second'].split(',')[:-1]
    exon_motif_left=[x['exon_leftSeq_first'].split(',')[0]+x['exon_rightSeq_first'].split(',')[-1]]
    exon_motif_right=[x['exon_leftSeq_second'].split(',')[0]+x['exon_rightSeq_second'].split(',')[-1]]
    strand_first='U'
    strand_second='U'
    if len(exon_leftSeq_first)>0:
        for i in range(len(exon_leftSeq_first)):
            exon_motif_left.append(exon_leftSeq_first[i]+exon_rightSeq_first[i])
    if len(exon_leftSeq_second)>0:
        for i in range(len(exon_leftSeq_second)):
            exon_motif_right.append(exon_leftSeq_second[i]+exon_rightSeq_second[i])
    for i in exon_motif_left:
        if i in ['AGGT','AGGC']:
            strand_first='+'
            break
        elif i in ['ACCT','GCCT']:
            strand_first='-'
            break
    for i in exon_motif_right:
        if i in ['AGGT','AGGC']:
            strand_second='+'
            break
        elif i in ['ACCT','GCCT']:
            strand_second='-'
            break
            
    if strand_first == 'U':
        if x['chr_first']  in tabixfile.contigs:
            for gtf in tabixfile.fetch(x['chr_first'], x['start_first']-1, x['start_first'],parser=pysam.asGTF()):
                if gtf.start==x['start_first']-1:
                    strand_first=gtf.strand
                    break
    if strand_first == 'U':
        if x['chr_first']  in tabixfile.contigs:
            for gtf in tabixfile.fetch(x['chr_first'], x['end_first']-1, x['end_first'],parser=pysam.asGTF()):
                if gtf.end==x['end_first']-1:
                    strand_first=gtf.strand
                    break
    if strand_second == 'U':
        if x['chr_second']  in tabixfile.contigs:
            for gtf in tabixfile.fetch(x['chr_second'], x['start_second']-1, x['start_second'],parser=pysam.asGTF()):
                if gtf.start==x['start_second']-1:
                    strand_second=gtf.strand
                    break
    if strand_second == 'U':
        if x['chr_second']  in tabixfile.contigs:
            for gtf in tabixfile.fetch(x['chr_second'], x['end_second']-1, x['end_second'],parser=pysam.asGTF()):
                if gtf.end==x['end_second']-1:
                    strand_second=gtf.strand
                    break
    return([strand_first,strand_second])
def getGene(chr,start,end,strand):
    tmp={'chr':chr,'start':start,'end':end,'strand':strand}
    geneList1=[]
    geneList2=[]
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
        geneCom=geneCom
    else:
        geneCom=[]
        for gtf in tabixfile.fetch(tmp['chr'], tmp['start']-1, tmp['end']-1,parser=pysam.asGTF()):
            if gtf.feature=='gene':
                gn=gtf.gene_name
                if gtf.strand==tmp['strand'] or tmp['strand'] == 'U':
                    geneCom.append(gn)
                else:
                    if gn[-4:]=='-AS1':
                        geneCom.append(gn[0:(len(gn)-4)])
                    else:
                        geneCom.append(gn+'-AS1')
    geneCom=list(set(geneCom))
    noAS=[]
    if len(geneCom)>1:
        for j in geneCom:
            if j[-4:]!='-AS1':
                noAS.append(j)
        if len(noAS)>0:
            geneCom=noAS
    return(';'.join(geneCom))

def fusion1CombinFL(options):
    global genome, tabixfile
    genomeFile=options[0]
    genome=pyfasta.Fasta(genomeFile)
    outPrefix=options[1]
    foutFa=open(outPrefix+'fusionSeq.fa','w')
    gtfFile=options[2]

    consFL_1=pd.read_csv(outPrefix+'constructFL_Fusion1_adj_1.txt',sep='\t')
    consFL_2=pd.read_csv(outPrefix+'constructFL_Fusion1_adj_2.txt',sep='\t')
    consFL_1['exon_start']=consFL_1['exon_start'].map(str).tolist()
    consFL_1['exon_end']=consFL_1['exon_end'].map(str).tolist()
    consFL_2['exon_start']=consFL_2['exon_start'].map(str).tolist()
    consFL_2['exon_end']=consFL_2['exon_end'].map(str).tolist()

    BSdf=pd.read_csv(outPrefix+'BS_Fusion1_adj.txt',sep='\t')
    explainF=pd.read_csv(outPrefix+"explainFL_Fusion1.txt",sep='\t')

    commonID=set(BSdf.ID) & set(consFL_1.ID ) & set(consFL_2.ID )

    consFL_1=consFL_1.set_index('ID').loc[commonID]
    consFL_2=consFL_2.set_index('ID').loc[commonID]
    BSdf=BSdf.set_index('ID').loc[commonID]
    explainF=explainF.drop_duplicates('ID').set_index('ID').loc[commonID]

    consFL_all=pd.concat([BSdf[['circID']],consFL_1[['chr','start','end','len','exonNum','exon_start','exon_end','exon_leftSeq','exon_rightSeq']],explainF['strand'],
                          consFL_2[['chr','start','end','len','exonNum','exon_start','exon_end','exon_leftSeq','exon_rightSeq']],BSdf[['strand']]],axis=1).copy()

    consFL_all.columns=['circID','chr_first','start_first','end_first','len_first','exonNum_first','exon_start_first','exon_end_first','exon_leftSeq_first','exon_rightSeq_first','strand_first',
                                 'chr_second','start_second','end_second','len_second','exonNum_second','exon_start_second','exon_end_second','exon_leftSeq_second','exon_rightSeq_second','strand_second']
    consFL_all['circID']=consFL_all['chr_first'].values+'|'+consFL_all['start_first'].map(str).values+'|'+consFL_all['end_first'].map(str).values+'|'+consFL_all['chr_second'].values+'|'+consFL_all['start_second'].map(str).values+'|'+consFL_all['end_second'].map(str).values

    strandDict={'+':'-','-':'+'}
    strand_second=[]
    if len(consFL_all.shape)>0:
        for i in range(consFL_all.shape[0]):
            if consFL_all.iloc[i]['strand_second']:
                strand_second.append(consFL_all.iloc[i]['strand_first'])
            else:
                strand_second.append(strandDict[consFL_all.iloc[i]['strand_first']])

    consFL_all['strand_second']=strand_second



    for i in range(consFL_all.shape[0]):
        tmp=consFL_all.iloc[i]
        seq1=getSeq(tmp['chr_first'],tmp['exon_start_first'],tmp['exon_end_first'],tmp['strand_first']).upper()
        seq2=getSeq(tmp['chr_second'],tmp['exon_start_second'],tmp['exon_end_second'],tmp['strand_second']).upper()
        foutFa.write('>'+tmp.name+'\n'+seq1+seq2+'\n')

    foutFa.close()
    cmd='TideHunter -f 2 -l -c 1 -t 4 '+outPrefix+'fusionSeq.fa>'+outPrefix+'fusionSeq.th 2>/dev/null'
    #print(cmd)
    os.system(cmd)
    fusionTH=pd.read_csv(outPrefix+'fusionSeq.th',sep='\t',names=['ID','consN','readLen','start','end','consLen','copyNum','full','consensus'])
    fusionTH['usage']=(fusionTH.end-fusionTH.start+1)/fusionTH.readLen
    fusionTH=fusionTH.loc[(fusionTH.copyNum>1.2) &(fusionTH.usage>0.5)]
    consFL_all_pass=consFL_all.loc[consFL_all.index.isin(set(consFL_all.index)-set(fusionTH['ID']))].copy()
    tabixfile = pysam.TabixFile(gtfFile)

    if len(options)<=4: # no stranded
        strandList=pd.DataFrame(consFL_all_pass.apply(getStrand,axis=1).to_list())
        consFL_all_pass['strand_first']=strandList.iloc[:,0].to_list()
        consFL_all_pass['strand_second']=strandList.iloc[:,1].to_list()



        
    geneNameList1=[]
    geneNameList2=[]
    for i in range(consFL_all_pass.shape[0]):
        tmp=consFL_all_pass.iloc[i]
        geneNameList1.append(getGene(tmp['chr_first'],tmp['start_first'],tmp['end_first'],tmp['strand_first']))
        geneNameList2.append(getGene(tmp['chr_second'],tmp['start_second'],tmp['end_second'],tmp['strand_second']))
        
    consFL_all_pass['geneName_first']=geneNameList1
    consFL_all_pass['geneName_second']=geneNameList2

    consFL_out=consFL_all_pass[['circID', 'chr_first', 'start_first', 'end_first', 'len_first',
           'exonNum_first', 'exon_start_first', 'exon_end_first',
           'exon_leftSeq_first', 'exon_rightSeq_first', 'strand_first','geneName_first',
           'chr_second', 'start_second', 'end_second', 'len_second',
           'exonNum_second', 'exon_start_second', 'exon_end_second',
           'exon_leftSeq_second', 'exon_rightSeq_second', 'strand_second','geneName_second']]
           
    consFL_out.to_csv(outPrefix+'result_fusion1.txt',sep='\t')