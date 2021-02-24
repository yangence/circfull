import numpy as np,pandas as pd,sys
def createStrandFq(RG_dir,outPrefix,fastqFile,gtfFile):
    predictStrandDf=pd.read_csv(outPrefix+'strandProbability.txt',sep='\t',names=['ID','p'])
    predictStrandDf=predictStrandDf.sort_values('ID').set_index('ID')
    gtf = pd.read_csv(gtfFile,sep="\t",names=['chr','source','type','start','end','score','strand','phase','attributes'],comment='#')
    gtf_exon=gtf[gtf['type']=='exon'].sort_values(by=['chr','start'])
    exonSdict={}
    exonEdict={}
    for i in range(gtf_exon.shape[0]):
        tmp=gtf_exon.iloc[i,]
        tmpS=tmp['chr']+'|'+str(tmp['start'])
        tmpE=tmp['chr']+'|'+str(tmp['end'])
        exonSdict[tmpS]=tmp['strand']
        exonEdict[tmpE]=tmp['strand']
        
    exonDict=dict(exonSdict,**exonEdict)
    explainFL=pd.read_csv(RG_dir+'explainFL_Normal_adj.txt',sep='\t')
    explainFL=explainFL.drop_duplicates('ID')
    explainFL=explainFL.sort_values('ID')
    explainFL.index=explainFL['ID']

    consFL=pd.read_csv(RG_dir+"constructFL_Normal_adj.txt",sep='\t')
    consFL.index=consFL['ID']

    consFL_exonKey={}
    consFL_cluster={}
    consFL_c2exon={}
    for i in range(consFL.shape[0]):
        tmp=consFL.iloc[i]
        exonS=tmp['exon_start'].split(',')
        exonE=tmp['exon_end'].split(',')
        chr=tmp['chr']
        tmpK=[chr+'|'+i for i in exonS+exonE]
        tmpRead=tmp['ID']
        clusterName=[]
        for j in tmpK:
            if consFL_exonKey.__contains__(j):
                clusterName.append(consFL_exonKey[j])
        clusterName=list(set(clusterName))
        if len(clusterName)==0:
            clusterName='C'+str(i)
            consFL_c2exon[clusterName]=[]
            for j in tmpK:
                consFL_exonKey[j]=clusterName
                consFL_c2exon[clusterName].append(j)
            consFL_cluster[clusterName]=[tmpRead]
        elif len(clusterName)==1:
            clusterName=clusterName[0]
            for j in tmpK:
                consFL_exonKey[j]=clusterName
                consFL_c2exon[clusterName].append(j)
            consFL_c2exon[clusterName]=list(set(consFL_c2exon[clusterName]))
            consFL_cluster[clusterName].append(tmpRead)
        else:
            for j in range(1,len(clusterName)):
                consFL_cluster[clusterName[0]].extend(consFL_cluster[clusterName[j]])
                for p in consFL_c2exon[clusterName[j]]:
                    consFL_exonKey[p]=clusterName[0]
                    consFL_c2exon[clusterName[0]].append(p)
                consFL_c2exon[clusterName[j]]=[]
                consFL_cluster[clusterName[j]]=[]
            
            clusterName=clusterName[0]
            for j in tmpK:
                consFL_exonKey[j]=clusterName
                consFL_c2exon[clusterName].append(j)
            consFL_c2exon[clusterName]=list(set(consFL_c2exon[clusterName]))
            consFL_cluster[clusterName].append(tmpRead)
            
    consFL_cluster_name=list(consFL_cluster.keys())

    cluster_strand={}
    for i in range(len(consFL_cluster)):
        clusterName=consFL_cluster_name[i]
        if len(consFL_cluster[clusterName])==0:
            continue
        strandScore=0 # gene strand origin
        clusterRead=consFL_cluster[clusterName]
        clusterExon=consFL_c2exon[clusterName]
        for j in clusterRead:
            mapS=explainFL.loc[j,'strand']
            rawS=predictStrandDf.loc[j.split('_')[0],'p']
            if mapS=='+':
                if rawS>0.5:
                    strandScore+=1
                else:
                    strandScore+=-1
            else:
                if rawS>0.5:
                     strandScore+=-1
                else:
                    strandScore+=1
        for j in clusterExon:
            if exonDict.__contains__(j):
                if exonDict[j]=='+':
                    strandScore=1
                else:
                    strandScore=-1
        cluster_strand[clusterName]=strandScore
        
    zeroCluster=list(np.where(np.array(cluster_strand.values())==0)[0])
    consFL_cluster_name=list(cluster_strand.keys())
    if len(zeroCluster)>0:
        for i in zeroCluster:
            consFL_cluster_name.pop(i)

    circStrand=open(outPrefix+'circOrignStrand.txt','w')
    for i in range(len(consFL_cluster_name)):
        strandScore=cluster_strand[consFL_cluster_name[i]]
        clusterRead=consFL_cluster[consFL_cluster_name[i]]
        for j in clusterRead:
            circStrand.write(j+'\t'+str(strandScore)+'\n')
    circStrand.close()

    readStrand={}
    for i in range(len(consFL_cluster_name)):
        strandScore=cluster_strand[consFL_cluster_name[i]]
        clusterRead=consFL_cluster[consFL_cluster_name[i]]
        for j in clusterRead:
            mapS=explainFL.loc[j,'strand']
            if mapS=='+':
                if strandScore>0:
                    readStrand[j]=1
                else:
                    readStrand[j]=0
            else:
                if strandScore>0:
                    readStrand[j]=0
                else:
                    readStrand[j]=1


    oldFq=open(fastqFile)
    newFq=open(outPrefix+'strandedFq.fastq','w')

    while True:
        fq1=oldFq.readline()
        if fq1:
            fq2=oldFq.readline().strip()
            fq3=oldFq.readline()
            fq4=oldFq.readline().strip()
            ID=fq1.split('\t')[0]
            ID=ID.split(' ')[0].split('_')[0][1:]
            if readStrand.__contains__(ID):
                pro=readStrand[ID]
                if pro==0:
                    fq2=fq2[::-1].replace('A','t').replace('T','a').replace('C','g').replace('G','c').upper()
                    fq4=fq4[::-1]
                newFq.write(fq1+fq2+'\n'+fq3+fq4+'\n')   
        else:
            break

    oldFq.close()
    newFq.close()
