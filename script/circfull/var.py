'''
Usage: circfull var -f fastq -g genome [-r RG_dir] [-t threads] [-o output]

Options:
    -h --help                   Show help message.
    -v --version                Show version.
    -f fastq                    circFL-seq fastq file.
    -g genome                   Fasta file of genome.
    -r RG_dir                   RG output dir [default: circFL_out].
    -t threads                  Number of threads [default: 20].
    -o output                   Output dir [default: circFL_out].
'''

from .genericFun import *
from .var_getVar import getVar

def var(options):
    fastq=options['-f']
    RG_dir=options['-r']
    thread=int(options['-t'])
    genome=options['-g']
    outDir=options['-o']
    if not os.path.exists(outDir):
        os.mkdir(outDir)
    if outDir[-1]!='/':
        outPrefix=outDir+'/'
    else:
        outPrefix=outDir
    var_outPrefix=outPrefix+'var/'
    createDir(var_outPrefix);createDir(var_outPrefix+'tmp')
    plog('Check fastq')
    checkFastq(fastq)
    plog('Check genome file')
    readFaFile(genome)
    plog('Check RG directory')
    fileCheck(RG_dir)
    if RG_dir[-1]!='/':
        RG_dir=RG_dir+'/'
    RG_file=False
    if os.path.exists(RG_dir+'result_Normal.txt'):
        RG_file=True
    if os.path.exists(RG_dir+'RG/result_Normal.txt'):
        RG_file=True
        RG_dir=RG_dir+'RG/'
    if not RG_file:
        sys.exit('ERROR: make sure file of result_Normal.txt in %s or %sRG/' % (RG_dir,RG_dir))
    plog('Identify variants: getVar')
    getVar(genome,RG_dir,var_outPrefix,fastq,thread)
