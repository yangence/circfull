'''
Usage: circfull cRG -f DNSC -g genome -a anno [-t threads] [-u] [-o output]

Options:
    -h --help                   Show help message.
    -v --version                Show version.
    -f DNSC                     DNSC directory.
    -g genome                   Fasta file of genome.
    -a anno                     Tabix indexed gtf file of gene annotation.
    -t threads                  Number of threads [default: 20].
    -u                          Stranded reads.
    -o output                   Output dir [default: circFL_out].
'''

from .genericFun import *
from .RG import RG
def createFastq(th,outPrefix):
    fin=open(th)
    fout=open(outPrefix+'pseudo.fq','w')
    while True:
        line=fin.readline()
        if line:
            each=line.strip().split('\t')
            if each[0]=='ID':
                continue
            id=each[0]
            seq=each[8]
            seqLen=len(seq)
            qual=''.join(['F' for i in range(seqLen)])
            fout.write('>'+id+'\n'+seq+seq+seq+'\n+\n'+qual+qual+qual+'\n')
            line=fin.readline()
        else:
            fin.close()
            break
def cRG(options):
    outDir=options['-o']
    DNSC=options['-f']
    if not os.path.exists(outDir):
        os.mkdir(outDir)
    if outDir[-1]!='/':
        outPrefix=outDir+'/'
    else:
        outPrefix=outDir
    cRG_outPrefix=outPrefix+'cRG/'
    createDir(cRG_outPrefix)
    plog('Check DNSC directory')
    fileCheck(DNSC)
    if DNSC[-1]!='/':
        DNSC=DNSC+'/'
    thpFile=False
    if os.path.exists(DNSC+'TideHunter_Pass.tab'):
        thpFile=DNSC+'TideHunter_Pass.tab'
    if os.path.exists(DNSC+'DNSC/TideHunter_Pass.tab'):
        thpFile=DNSC+'DNSC/TideHunter_Pass.tab'
    if not thpFile:
        sys.exit('ERROR: make sure file of TideHunter_Pass.tab in %s' % DNSC)
    
    plog('Make query sequence: createFastq')
    createFastq(thpFile,cRG_outPrefix)
    
    options['-f']=cRG_outPrefix+'pseudo.fq'
    options['-o']=cRG_outPrefix
    options['-r']=False
    options['-m']=False
    RG(options)
    