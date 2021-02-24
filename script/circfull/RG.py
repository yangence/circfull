'''
Usage: circfull RG -f fastq -g genome -a anno [-t threads] [-r] [-m rmsk] [-o output]

Options:
    -h --help                   Show help message.
    -v --version                Show version.
    -f fastq                    circFL-seq fastq file.
    -g genome                   Fasta file of genome.
    -a anno                     Tabix indexed gtf file of gene annotation.
    -t threads                  Number of threads [default: 20].
    -r                          Filter out low quality circRNA.
    -m rmsk                     Filter out low quality circRNA with a huge overlap of repeat region. Only works with '-r'.
    -o output                   Output dir [default: circFL_out].
'''
from .genericFun import *
import sys,time,pandas as pd,numpy as np,docopt,os,time,pysam
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


def RG(options):
    isFilter=False
    if options['-r']:
        isFilter=True
    rmskFile=False
    if options['-m']:
        rmskFile=options['-m']
        plog('Check rmsk file')
        fileCheck(rmskFile)
    fastq=options['-f']
    genome=options['-g']
    anno=options['-a']
    thread=int(options['-t'])
    outDir=options['-o']
    plog('Check fastq file')
    checkFastq(fastq)
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
    RG_outPrefix=outPrefix+'RG/'
    RG_tmp_outPrefix=RG_outPrefix+'tmp/'
    fus_outPrefix=RG_outPrefix+'fusion/'
    fus_tmp_outPrefix=fus_outPrefix+'tmp/'
    createDir(outPrefix);createDir(RG_outPrefix);createDir(RG_tmp_outPrefix);createDir(fus_outPrefix);createDir(fus_tmp_outPrefix)
    sam=RG_outPrefix+'test.minimap2.sam'
    bam=RG_outPrefix+'test.minimap2.bam'
    
    plog('Align fastq to reference genome: alignFastq')
    alignFastq(fastq,genome,sam,thread)

    plog('Transform SAM to BAM: sam2bam')
    sam2bam(sam,bam)
    
    plog('Analyze SAM file: explainFL')
    explainFL(genome,RG_outPrefix,sam)
    
    plog('Filter and classify candidates: filterFL')
    filterFL(genome,RG_outPrefix,fastq,thread)
    
    plog('Adjust normal: adjExplainNormal')
    adjExplainNormal(genome,RG_outPrefix,thread)
    
    plog('Re-alignment to pseudo reference: detectBS')
    detectBS([genome,RG_outPrefix,fastq,thread])
    
    plog('Filter BS: filterBS')
    filterBS([genome,anno,RG_outPrefix,thread])
    
    plog('Construct FL-circRNA: constructFL')
    constructFL([genome,RG_outPrefix,thread])
    
    plog('Adjust FL-circRNA: adjFL')
    adjFL([genome,anno,RG_outPrefix,thread])
    plog('Adjust mistaken circRNA: circSeq')
    circSeq(genome,RG_outPrefix,thread)
    
    plog('Annotate circRNA: annotateNormal')
    annotateNormal([RG_outPrefix,anno])
    
    if os.path.exists(RG_outPrefix+'explainFL_Fusion1.txt')|os.path.exists(RG_outPrefix+'explainFL_Fusion2.txt'):
        plog('Analyze fusion circRNA: fusionFq')
        fusionFq(fastq,RG_outPrefix)
    
    if os.path.exists(RG_outPrefix+'explainFL_Fusion1.txt'):
        os.system('cp %s %s' %(RG_outPrefix+'explainFL_Fusion1.txt',fus_outPrefix+'explainFL_Fusion1.txt'))
        plog('Re-alignment to pseudo reference: detectBS_fusion1')
        BS_Fusion1=detectBS_fusion1([genome,fus_outPrefix,thread])
        if BS_Fusion1:
            plog('Filter fusion BS: filterBS_fusion1')
            filterBS_fusion1(genome,anno,fus_outPrefix)
            plog('Construct fusion FL-circRNA: fusion1ConstructFL')
            f1_cf_1=fusion1ConstructFL(genome,fus_outPrefix,'1',thread)
            f1_cf_2=fusion1ConstructFL(genome,fus_outPrefix,'2',thread)
            if f1_cf_1 and f1_cf_2:
                 plog('Adjust fusion FL-circRNA: fusion1AdjFL')
                 fusion1AdjFL([genome,anno,fus_outPrefix,'1',thread])
                 fusion1AdjFL([genome,anno,fus_outPrefix,'2',thread])
                 plog('Annotate fusion FL-circRNA: fusion1CombinFL')
                 fusion1CombinFL([genome,fus_outPrefix,anno])
            else:
                 plog('No fusion circRNAs from different chromosomes were found!')
        else:
            plog('No fusion circRNAs from different chromosomes were found!')
    if os.path.exists(RG_outPrefix+'explainFL_Fusion2.txt'):
        os.system('cp %s %s' %(RG_outPrefix+'explainFL_Fusion2.txt',fus_outPrefix+'explainFL_Fusion2.txt'))
        plog('Re-alignment to pseudo reference: detectBS_fusion2')
        BS_Fusion2=detectBS_fusion2([genome,fus_outPrefix,thread])
        if BS_Fusion2:
            plog('Filter fusion BS: filterBS_fusion2')
            filterBS_fusion2(genome,anno,fus_outPrefix)
            plog('Construct fusion FL-circRNA: fusion2ConstructFL')
            f2_cf_1=fusion2ConstructFL(genome,fus_outPrefix,'1',thread)
            f2_cf_2=fusion2ConstructFL(genome,fus_outPrefix,'2',thread)
            if f2_cf_1 and f2_cf_2:
                 plog('Adjust fusion FL-circRNA: fusion2AdjFL')
                 fusion2AdjFL([genome,anno,fus_outPrefix,'1',thread])
                 fusion2AdjFL([genome,anno,fus_outPrefix,'2',thread])
                 plog('Annotate fusion FL-circRNA: fusion2CombinFL')
                 fusion2CombinFL([genome,fus_outPrefix,anno])
            else:
                 plog('No fusion circRNAs in a same chromosome were found!')
        else:
            plog('No fusion circRNAs in a same chromosome were found!')
    
    plog('Format circFL results: circFL_output')
    circFL_output(RG_outPrefix)

    if isFilter:
        plog('Filter circFL results: filterOut')
        filterOut(RG_outPrefix,RG_outPrefix,thread,rmskFile)
    
    plog('All done!!!')
    