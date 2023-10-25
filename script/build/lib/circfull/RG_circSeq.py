import pandas as pd, numpy as np, pyfasta,sys,os
from multiprocessing import Pool


def getSeq(i):
    tmp=FL.iloc[i].copy()
    exonStart=[int(i) for i in tmp['exon_start'].split(',')]
    exonEnd=[int(i) for i in tmp['exon_end'].split(',')]
    chr=tmp['chr']
    seq=''
    for i in range(len(exonStart)):
        seq+=genome.sequence({'chr': chr, 'start':exonStart[i], 'stop':exonEnd[i]})
    return([tmp['ID'],seq])

def tidehunter(faFile,thFile,thread):
    cmd="TideHunter -f 2 -c 1.2 -l -t %i %s >%s 2>/dev/null" % (thread,faFile,thFile)
    os.system(cmd)


def circSeq(genomeFile,outPrefix,thread):
    global genome,FL
    genome=pyfasta.Fasta(genomeFile)
    FL=pd.read_csv(outPrefix+'constructFL_Normal_adj.txt',sep='\t',dtype={'exon_start':str,'exon_end':str})
    fout=open(outPrefix+'circSeq.fa','w')
    pool=Pool(processes=thread)
    seq=pool.map(getSeq,range(FL.shape[0]))
    pool.close()
    pool.join()
    for i in range(len(seq)):
        fout.write('>'+seq[i][0]+'\n'+seq[i][1]+'\n')
    fout.close()
    tidehunter(outPrefix+'circSeq.fa',outPrefix+'circSeq.th',thread)
    