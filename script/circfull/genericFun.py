import os,sys,time,pysam
import pandas as pd

def plog(x):
    t=time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    print('#### %s #### \n%s'   % (t,x))

def fileCheck(file):
    if not os.path.exists(file):
        sys.exit('ERROR: %s is not exist!!!' % file)
def checkFastq(file):
    if not os.path.exists(file):
        sys.exit('ERROR: %s is not exist!!!' % file)
    if file[-3:]=='.gz':
        sys.exit('ERROR: %s is compressed, please uncompress it first!!!' % file)
def readGTFfile(fileName):
    if not os.path.exists(fileName):
        sys.exit('ERROR: %s is not exist!!!' % fileName)
    try:
        tabixfile = pysam.TabixFile(fileName)
        tabixfile.close()
    except:
        sys.exit('ERROR: make sure %s is sorted and tabix indexed!!!' % fileName)

def readFaFile(fileName):
    if not os.path.exists(fileName):
        sys.exit('ERROR: %s is not exist!!!' % fileName)
    fff=open(fileName)
    fl=fff.readline()
    fff.close()
    if fl[0]!='>':
        sys.exit('ERROR:  %s need to be Fasta format!!!' % fileName)
    try:
        faFile=pysam.FastaFile(fileName)
        faFile.close()
    except:
        sys.exit('ERROR: make sure %s is fai indexed!!!' % fileName)
def createDir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)
        
def fusionFq(fastqFile,outPrefix):
    oldFq=open(fastqFile)
    newFq=open(outPrefix+'fusion/fusion1.fq','w')

    explainFL_ID2Type=pd.read_csv(outPrefix+'explainFL_ID2Type.txt',sep='\t')
    explainFL_ID2Type_can=explainFL_ID2Type.loc[[i in ['F1','F2','FC1','FC2'] for i in explainFL_ID2Type['type'].values]]
    canID={}
    for i in explainFL_ID2Type_can['ID'].values:
        canID[i]=1
    while True:
        fq1=oldFq.readline()
        if fq1:
            fq2=oldFq.readline()
            fq3=oldFq.readline()
            fq4=oldFq.readline()
            ID=fq1.split(' ')[0][1:]
            if canID.__contains__(ID):
                newFq.write(fq1+fq2+fq3+fq4)
        else:
            break

    oldFq.close()
    newFq.close()

def alignFastq(fastq,genome,sam,thread,strand=False):
    if strand:
        os.system('minimap2 -ax splice -uf -p 0.5 -t %i  %s %s>%s 2>/dev/null' % (thread, genome, fastq,sam) )
    else:
        os.system('minimap2 -ax splice  -p 0.5 -t %i  %s %s>%s 2>/dev/null' % (thread, genome, fastq,sam) )
def sam2bam(sam,bam):
    os.system('samtools sort -m 10G -@ 8 %s -o %s  2>/dev/null' % (sam,bam))
    os.system('samtools index %s  2>/dev/null' % bam)
