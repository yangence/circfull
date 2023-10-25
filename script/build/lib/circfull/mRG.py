'''
Usage: circfull mRG -f fastq -g genome [-m rmsk] [-r RG] [-c cRG] [-s sRG] [-t threads] [-o output]

Options:
    -h --help                   Show help message.
    -v --version                Show version.
    -f fastq                    circFL-seq fastq file.
    -g genome                   Fasta file of genome.
    -m rmsk                     Filter out low quality circRNA with a huge overlap of repeat region.
    -r RG                       RG directory [default: circFL_out].
    -c cRG                      cRG directory [default: circFL_out].
    -s sRG                      sRG directory [default: circFL_out].
    -t threads                  Number of threads [default: 20].
    -o output                   Output dir [default: circFL_out].
'''

from .genericFun import *
from .mRG_merge import merge
from .RG_filterOut import filterOut
def mRG(options):
    rmskFile=False
    if options['-m']:
        rmskFile=options['-m']
        plog('Check rmsk file')
        fileCheck(rmskFile)
    outDir=options['-o']
    RG=options['-r']
    cRG=options['-c']
    sRG=options['-s']
    fastq=options['-f']
    genome=options['-g']
    thread=int(options['-t'])
    plog('Check fastq file')
    checkFastq(fastq)
    plog('Check genome file')
    readFaFile(genome)
    if not os.path.exists(outDir):
        os.mkdir(outDir)
    if outDir[-1]!='/':
        outPrefix=outDir+'/'
    else:
        outPrefix=outDir
    createDir(outPrefix+'mRG/')
    createDir(outPrefix+'mRG/tmp/')
    if RG[-1]!='/':
        RG=RG+'/'
    if cRG[-1]!='/':
        cRG=cRG+'/'
    if sRG[-1]!='/':
        sRG=sRG+'/'

    plog('Check RG directory')
    fileCheck(RG)
    RG_N=False
    if os.path.exists(RG+'circFL_Normal.txt'):
        RG_N=True
    if os.path.exists(RG+'RG/circFL_Normal.txt'):
        RG_N=True
        RG=RG+'RG/'
    if not RG_N:
        sys.exit('ERROR: make sure file of circFL_Normal.txt in %s or %sRG/' % (RG,RG))
    
    plog('Check cRG directory')
    #fileCheck(cRG)
    cRG_N=False
    if os.path.exists(cRG+'circFL_Normal.txt'):
        cRG_N=True
    if os.path.exists(cRG+'cRG/circFL_Normal.txt'):
        cRG_N=True
        cRG=cRG+'cRG/'
    if os.path.exists(cRG+'cRG/RG/circFL_Normal.txt'):
        cRG_N=True
        cRG=cRG+'cRG/RG/'
    if not cRG_N:
        plog('Warning: make sure file of circFL_Normal.txt in %s or %scRG/' % (cRG,cRG))

    plog('Check sRG directory')
    #fileCheck(sRG)
    sRG_N=False
    if os.path.exists(sRG+'circFL_Normal.txt'):
        sRG_N=True
    if os.path.exists(sRG+'sRG/circFL_Normal.txt'):
        sRG_N=True
        sRG=sRG+'sRG/'

    if not sRG_N:
        plog('Warning: make sure file of circFL_Normal.txt in %s or %ssRG/' % (sRG,sRG))
    if (not cRG_N) and (not sRG_N):
        #sys.exit('ERROR: at least one of cRG and sRG directory contain circFL_Normal.txt')
        plog('No cRG and sRG found: merge')
        merge(RG=RG,cRG=RG,outDir=outPrefix+'mRG/')
    else:
        if cRG_N and (not sRG_N):
            plog('Merge results from RG and cRG: merge')
            merge(RG=RG,cRG=cRG,outDir=outPrefix+'mRG/')
        elif sRG_N and (not cRG_N):
            plog('Merge results from RG and sRG: merge')
            merge(RG=RG,sRG=sRG,outDir=outPrefix+'mRG/')
        else:
            plog('Merge results from RG ,cRG and sRG: merge')
            merge(RG=RG,cRG=cRG,sRG=sRG,outDir=outPrefix+'mRG/')
    
    if not os.path.exists(RG+'test.minimap2.bam'):
        plog('Transform SAM to BAM: sam2bam')
        sam2bam(RG+'test.minimap2.sam',RG+'test.minimap2.bam')
    
    plog('Filter circFL results: filterOut')
    filterOut(RG,outPrefix+'mRG/',fastq,genome,thread,rmskFile)
    
    plog('All done!!!')