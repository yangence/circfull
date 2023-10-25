'''
Usage: circfull DNSC -f fastq  [-t threads] [-c] [-m mem] [-o output]

Options:
    -h --help                   Show help message.
    -v --version                Show version.
    -f fastq                    circFL-seq fastq file.
    -c                          Perform cluster of consensus sequence.
    -m mem                      The maximum memory usage, only worked when -c is designated. [default: 100G]
    -t threads                  Number of threads [default: 20].
    -o output                   Output dir [default: circFL_out].
'''
import os,sys,time
from .genericFun import *
from .DNSC_getNovo import getNovo
from .DNSC_novoCluster import novoCluster

def fastq2fa(fq,fa,outPrefix):
    fin=open(fq)
    fout=open(fa,'w')
    line1=fin.readline()
    while line1:
        line2=fin.readline()
        line3=fin.readline()
        line4=fin.readline()
        fout.write('>'+line1[1:])
        fout.write(line2)
        line1=fin.readline()

def tideHunter(fa,th,thread):
    os.system('TideHunter -f 2 -c 1.5 -p 30 -l -t %i %s>%s' % (thread,fa,th))

def read2read(outPrefix,thread):
    os.system('minimap2 -x  ava-ont -t %i %s %s >%s 2>/dev/null' % (thread,outPrefix+'rawseq.fa',outPrefix+'rawseq.fa',outPrefix+'raw2raw.paf'))

def sortPaf(outPrefix,outPrefix_tmp,thread,mem):
    os.system('sort -k 1,1V -k 6,6V   --parallel=%i   --temporary-directory=%s  -S %s   %s  >%s 2>/dev/null' \
    % (thread,outPrefix_tmp,mem,outPrefix+'raw2raw.paf',outPrefix+'raw2raw.sort.paf'))
def createDir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)

def DNSC(options):
    isCluster=False
    if options['-c']:
        isCluster=True
    fastq=options['-f']
    mem=options['-m']
    thread=int(options['-t'])
    outDir=options['-o']
    plog('Check fastq file')
    checkFastq(fastq)
    if not os.path.exists(outDir):
        os.mkdir(outDir)
    if outDir[-1]!='/':
        outPrefix=outDir+'/'
    else:
        outPrefix=outDir
    DNSC_outPrefix=outPrefix+'DNSC/'
    DNSC_tmp_outPrefix=DNSC_outPrefix+'tmp/'
    createDir(outPrefix);createDir(DNSC_outPrefix);createDir(DNSC_tmp_outPrefix)
    fa=DNSC_outPrefix+'test_1.fa'
    th=DNSC_outPrefix+'TideHunter.tab'
    
    plog('Transform fastq to fasta: fastq2fa')
    fastq2fa(fastq,fa,DNSC_outPrefix)
    plog('Detect consensue sequence: tideHunter')
    tideHunter(fa,th,thread)
    
    plog('Filter consensue sequence: getNovo')
    getNovo(DNSC_outPrefix,fastq,thread,isCluster)
    if isCluster:
        plog('Reads to reads: read2read')
        read2read(DNSC_outPrefix,thread)
        plog('Sort Paf: sortPaf')
        sortPaf(DNSC_outPrefix,DNSC_tmp_outPrefix,thread,mem)
        plog('Cluster TCS: novoCluster')
        novoCluster(DNSC_outPrefix)
    plog('All done!!!')