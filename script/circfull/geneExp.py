'''
Usage: circfull geneExp -f fastq -r ref [-t threads] [-o output]

Options:
    -h --help                   Show help message.
    -v --version                Show version.
    -f fastq                    circFL-seq fastq file.
    -r ref                      Gene reference
    -t threads                  Number of threads [default: 20].
    -o output                   Output dir [default: circFL_out].
'''
from .genericFun import *
import os,sys,time
import pandas as pd
def getGeneID(x):
    if x[0:7]=='gene_id':
        geneID=x.split('|')[0].split(',')[1][1:-1]
        geneName=x.split('|')[2].split(':')[0].split(',')[1][1:-1]
        return(geneID)
    else:
        return(x)

def geneName(x):
    if x[0:7]=='gene_id':
        geneName=x.split('|')[2].split(':')[0].split(',')[1][1:-1]
        return(geneName)
    else:
        return(x)

def geneCount(outPrefix):
    geneMap=pd.read_csv(outPrefix+'geneMap.paf',sep='\t',header=None,usecols=[0,5])
    geneMap.columns=['ID','anno']
    geneMap_filter=geneMap.drop_duplicates('ID').copy()
    geneMap_filter['geneID']=geneMap_filter.iloc[:,1].map(getGeneID)
    geneMap_filter['geneName']=geneMap_filter.iloc[:,1].map(geneName)
    geneCount=geneMap_filter['geneID'].value_counts()
    geneMap_filter.index=geneMap_filter.geneID
    geneCountDf=pd.DataFrame(geneCount,columns=['geneID','num'])
    geneCountDf['num']=geneCountDf['geneID']
    geneCountDf['geneID']=geneCountDf.index
    geneID2geneName=dict(zip(geneMap_filter['geneID'],geneMap_filter['geneName']))
    geneCountDf['geneName']=[geneID2geneName[i] for i in geneCountDf.geneID.to_list()]
    geneCountDf['ratio']=geneCountDf['num']/sum(geneCountDf['num'])
    geneCountDf=geneCountDf.loc[:,['geneID','geneName','num','ratio']]
    geneCountDf.to_csv(outPrefix+'geneCountDf.txt',index=None,sep='\t')


def map2ref(ref,fastq,thread,outPrefix):
    os.system("minimap2 -x splice -t %i %s %s >%s 2>/dev/null" % (thread,ref,fastq,outPrefix+'geneMap.paf'))

def geneExp(options):
    fastq=options['-f']
    thread=int(options['-t'])
    ref=options['-r']
    outDir=options['-o']
    plog('Check gene reference')
    fileCheck(ref)
    plog('Check fastq file')
    checkFastq(fastq)
    if not os.path.exists(outDir):
        os.mkdir(outDir)
    if outDir[-1]!='/':
        outPrefix=outDir+'/'
    else:
        outPrefix=outDir
    geneExp_outPrefix=outPrefix+'geneExp/'
    geneExp_tmp_outPrefix=geneExp_outPrefix+'tmp/'
    createDir(geneExp_outPrefix);createDir(geneExp_tmp_outPrefix)
    
    plog('Align fastq to reference: map2ref')
    map2ref(ref,fastq,thread,geneExp_outPrefix)
    
    plog('Count gene expression: geneCount')
    geneCount(geneExp_outPrefix)
    