'''
Usage: circfull filter -r dir -f fastq -g genome  -c circ [-m rmsk] [-t threads] [-o output] [-k]

Options:
    -h --help                   Show help message.
    -v --version                Show version.
    -r dir                      RG directory.
    -f fastq                    circFL-seq fastq file.
    -g genome                   Fasta file of genome.
    -m rmsk                     Filter out low quality circRNA with a huge overlap of repeat region.
    -c circ                     circFL full-length file.
    -t threads                  Number of threads [default: 20].
    -o output                   Output dir [default: circFL_filter].
    -k                          Skip calcuation of splicing ratio
'''
from .genericFun import *
from .file_filterOut import file_filterOut
def filter(options):
    rmskFile=False
    skipSplice=False
    circ_file=options['-c']
    fastq=options['-f']
    genome=options['-g']
    thread=int(options['-t'])
    RG=options['-r']
    plog('Check fastq file')
    checkFastq(fastq)
    plog('Check genome file')
    readFaFile(genome)
    plog('Check circFL file')
    fileCheck(circ_file)
    
    if options['-m']:
        rmskFile=options['-m']
        plog('Check rmsk file')
        fileCheck(rmskFile)
    if options['-k']:
        skipSplice=True
    outDir=options['-o']
    if not os.path.exists(outDir):
        os.mkdir(outDir)
    if outDir[-1]!='/':
        outPrefix=outDir+'/'
    else:
        outPrefix=outDir
    if RG[-1]!='/':
        RG=RG+'/'
    createDir(outPrefix+'/tmp')
    plog('Filter circFL results: filterOut')
    file_filterOut(RG,circ_file,outPrefix,fastq,genome,thread,rmskFile,skipSplice=skipSplice)
    
    
    
