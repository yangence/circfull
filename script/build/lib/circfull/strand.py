'''
Usage: circfull strand -f fastq -F fastq2 -a anno [-r RG_dir] [-l hang_len] [-t threads] [-o output]

Options:
    -h --help                   Show help message.
    -v --version                Show version.
    -f fastq                    Raw circFL-seq fastq file with adapter remain.
    -F fastq2                   Trimmed circFL-seq fastq file.
    -a anno                     Tabix indexed gtf file of gene annotation.
    -r RG_dir                   RG output dir. [default: circFL_out]
    -l hang_len                 Length of sequence to find primers in both end. [default: 100]
    -t threads                  Number of threads [default: 20].
    -o output                   Output dir [default: circFL_out].
'''
from .genericFun import *
from .strand_identifyStrand import identifyStrand
from .strand_judgeStrand import judgeStrand
from .strand_predictStrand import predictStrand
from .strand_createStrandFq import createStrandFq

def strand(options):
    fastq=options['-f']
    fastq2=options['-F']
    anno=options['-a']
    RG_dir=options['-r']
    thresholdLen=int(options['-l'])
    thread=int(options['-t'])
    outDir=options['-o']
    if not os.path.exists(outDir):
        os.mkdir(outDir)
    if outDir[-1]!='/':
        outPrefix=outDir+'/'
    else:
        outPrefix=outDir
    strand_outPrefix=outPrefix+'strand/'
    createDir(strand_outPrefix)
    plog('Check fastq')
    checkFastq(fastq)
    plog('Check fastq2')
    checkFastq(fastq2)
    plog('Check anno file')
    readGTFfile(anno)
    
    plog('Check RG directory')
    fileCheck(RG_dir)
    
    if RG_dir[-1]!='/':
        RG_dir=RG_dir+'/'
    RG_file=False
    if os.path.exists(RG_dir+'explainFL_Normal_adj.txt') and os.path.exists(RG_dir+'constructFL_Normal_adj.txt'):
        RG_file=True
    if os.path.exists(RG_dir+'RG/explainFL_Normal_adj.txt') and os.path.exists(RG_dir+'RG/constructFL_Normal_adj.txt'):
        RG_file=True
        RG_dir=RG_dir+'RG/'
    if not RG_file:
        sys.exit('ERROR: make sure file of explainFL_Normal_adj.txt and constructFL_Normal_adj.txt in %s or %sRG/' % (RG_dir,RG_dir))
    
    if not os.path.exists(strand_outPrefix+'primer_score.txt'):
        plog('Get strand match score: identifyStrand')
        identifyStrand(fastq,strand_outPrefix,thread,thresholdLen)
    
    plog('Prepare for training set: judgeStrand')
    judgeStrand(RG_dir,strand_outPrefix)
    plog('Predict strand: predictStrand')
    predictStrand(strand_outPrefix,fastq)
    plog('Create stranded fastq: createStrandFq')
    createStrandFq(RG_dir,strand_outPrefix,fastq2)
    
    plog('All done!!!')