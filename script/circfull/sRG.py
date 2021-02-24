'''
Usage: circfull sRG -g genome -a anno [-s strandDir] [-t threads] [-r] [-o output]

Options:
    -h --help                   Show help message.
    -v --version                Show version.
    -g genome                   Fasta file of genome.
    -a anno                     Tabix indexed gtf file of gene annotation.
    -t threads                  Number of threads [default: 20].
    -s dir                      Strand dir [default: circFL_out].
    -r                          Filter out false discovery.
    -o output                   Output dir [default: circFL_out].
'''

import sys,time,pandas as pd,numpy as np,docopt,os,time,pysam
from .genericFun import *
from .RG_explainFL import explainFL
from .RG_filterFL import filterFL
from .RG_adjExplainNormal import adjExplainNormal
from .RG_detectBS import detectBS
from .RG_filterBS import filterBS
from .RG_constructFL import constructFL
from .RG_adjFL import adjFL
from .RG_circSeq import circSeq
from .RG_annotateNormal import annotateNormal
from .RG_detectBS_fusion1 import detectBS_fusion1
from .RG_detectBS_fusion2 import detectBS_fusion2
from .RG_filterBS_fusion1 import filterBS_fusion1
from .RG_filterBS_fusion2 import filterBS_fusion2
from .RG_fusion1ConstructFL import fusion1ConstructFL
from .RG_fusion2ConstructFL import fusion2ConstructFL
from .RG_fusion1AdjFL import fusion1AdjFL
from .RG_fusion2AdjFL import fusion2AdjFL
from .RG_fusion1CombinFL import fusion1CombinFL
from .RG_fusion2CombinFL import fusion2CombinFL
from .RG_circFL_output import circFL_output
from .RG_filterOut import filterOut


def sRG(options,cRG=False):
    isFilter=False
    if options['-r'] or cRG:
        isFilter=True
    
    genome=options['-g']
    anno=options['-a']
    thread=int(options['-t'])
    outDir=options['-o']

    plog('Check anno file')
    readGTFfile(anno)
    plog('Check genome file')
    readFaFile(genome)
    if not os.path.exists(outDir):
        os.mkdir(outDir)
    if outDir[-1]!='/':
        outPrefix=outDir+'/'
    else:
        outPrefix=outDir
    RG_outPrefix=outPrefix+'sRG/'
    RG_tmp_outPrefix=RG_outPrefix+'tmp/'
    fus_outPrefix=RG_outPrefix+'fusion/'
    fus_tmp_outPrefix=fus_outPrefix+'tmp/'
    createDir(outPrefix);createDir(RG_outPrefix);createDir(RG_tmp_outPrefix);
    strandDir=options['-s']+'/'
    strand_outPrefix=strandDir
    strandIn=False
    plog('Check strandDir file')
    if os.path.exists(strandDir+'strandProbability.txt'):
        strandIn=True
        strand_outPrefix=strandDir
    if os.path.exists(strandDir+'strand/strandProbability.txt'):
        strandIn=True
        strand_outPrefix=strandDir+'strand/'
    if not strandIn:
        sys.exit('ERROR: strandProbability.txt not in %s or %s' % (strandDir,strandDir+'strand/'))
    
    
    fastq=strand_outPrefix+'strandedFq.fastq'
    plog('Check fastq file')
    checkFastq(fastq)
    
    sam=RG_outPrefix+'test.minimap2.sam'
    bam=RG_outPrefix+'test.minimap2.bam'
    
    plog('Align fastq to reference genome: alignFastq')
    alignFastq(fastq,genome,sam,thread,True)

    plog('Transform SAM to BAM: sam2bam')
    sam2bam(sam,bam)

    plog('Analyze SAM file: explainFL')
    explainFL(genome,RG_outPrefix,sam)
    plog('Filter and classify candidates: filterFL')
    filterFL(genome,RG_outPrefix,fastq,thread)
    
    plog('Adjust normal: adjExplainNormal')
    strandFile=adjExplainNormal(genome,RG_outPrefix,thread,True)
    #strandFile.to_csv(RG_outPrefix+'strandFile.txt',sep="\t",header=True,index=False)
    
    plog('Re-alignment to pseudo reference: detectBS')
    detectBS([genome,RG_outPrefix,fastq,thread,True])
    
    #strandFile=pd.read_csv(RG_outPrefix+'strandFile.txt',sep="\t")
    plog('Filter BS: filterBS')
    filterBS([genome,anno,RG_outPrefix,thread,strandFile])
    plog('Construct FL-circRNA: constructFL')
    constructFL([genome,RG_outPrefix,thread,strandFile])
    
    plog('Adjust FL-circRNA: adjFL')
    adjFL([genome,anno,RG_outPrefix,thread,strandFile])
    plog('Adjust mistaken circRNA: circSeq')
    circSeq(genome,RG_outPrefix,thread)
    plog('Annotate circRNA: annotateNormal')
    annotateNormal([RG_outPrefix,anno,strandFile])
    
    plog('Format circFL results: circFL_output')
    circFL_output(RG_outPrefix)

    if isFilter:
        plog('Filter circFL results: filterOut')
        filterOut(RG_outPrefix,RG_outPrefix,thread)
    
    plog('All done!!!')