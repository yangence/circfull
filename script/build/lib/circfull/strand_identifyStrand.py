import numpy as np, pandas as pd, os,sys
from multiprocessing import Pool
from progressbar import *
from ssw import Aligner
matchS = 3; mismatchS = -6; gapopenS = -5; gapextendS = -2
    
def getALL(seq1,seq2,isSeq=False):
    read_align=Aligner(reference=seq1,molecule='dna',gap_open=5,gap_extend=2)
    read_align.matrix.match=3
    read_align.matrix.mimatch=-6
    read_align_result=read_align.align(seq2,revcomp=False)
    return(read_align_result.score)

def getScoremat(seq):
    [ID,leftSeq,rightSeq] = seq
    score=[getALL(AnchorX,leftSeq),getALL(AnchorY,leftSeq),getALL(polyT,leftSeq),getALL(AnchorX_rev,rightSeq),getALL(AnchorY_rev,rightSeq),getALL(polyA,rightSeq)]
    return(ID+'\t'+'\t'.join([str(i) for i in score])+'\n')
    #print((ID+'\t'+'\t'.join([str(i) for i in score])+'\n'))
    #fout_score.write(ID+'\t'+'\t'.join([str(i) for i in score])+'\n')

def identifyStrand(rawFastq,outPrefix,thread,ThresholdLen=100):
    global AnchorX,AnchorY,AnchorX_rev,AnchorY_rev,polyT,polyA,fout_score
    scoreDict={}
    scoreFileName=outPrefix+"primer_score.txt"
    fout_score=open(scoreFileName,'w')
    #if os.path.exists(scoreFileName):
    #
    #    scoreFile=pd.read_csv(scoreFileName,sep='\t')
    #    for i in range(scoreFile.shape[0]):
    #        scoreDict[scoreFile.iloc[i,0]]=1
    AnchorX='GTCGACGGCGCGCCGGATCCATA'
    AnchorY='ATATCTCGAGGGCGCGCCGGATCC'
    AnchorX_rev=AnchorX[::-1].replace('A','t').replace('T','a').replace('C','g').replace('G','c').upper()
    AnchorY_rev=AnchorY[::-1].replace('A','t').replace('T','a').replace('C','g').replace('G','c').upper()
    polyT='TTTTTTTTTTTTTTTTTTTTTTTT'
    Barcode={'BC1':'CACAAAGACACCGACAACTTTCTT',
    'BC2':'ACAGACGACTACAAACGGAATCGA',
    'BC3':'CCTGGTAACTGGGACACAAGACTC',
    'BC4':'TAGGGAAACACGATAGAATCCGAA',
    'BC5':'AAGGTTACACAAACCCTGGACAAG',
    'BC6':'GACTACTTTCTGCCTTTGCGAGAA',
    'BC7':'AAGGATTCATTCCCACGGTAACAC',
    'BC8':'ACGTAACTTGGTTTGTTCCCTGAA',
    'BC9':'AACCAAGACTCGCTGTGCCTAGTT',
    'BC10':'GAGAGGACAAAGGTTTCAACGCTT',
    'BC11':'TCCATTCCCTCCGATAGATGAAAC',
    'BC12':'TCCGATTCTGCTTCTTTCTACCTG'
    }
    adapter='ACGTATTGCT'
    #AnchorX=adapter+Barcode[BC]+AnchorX
    #AnchorY=adapter+Barcode[BC]+AnchorY

    polyA=polyT.replace('T','A')
    AnchorX_len=len(AnchorX)
    AnchorY_len=len(AnchorY)
    poly_len=len(polyT)
    fastqFile=open(rawFastq)
    seqList=[]
    while True:
        fq1=fastqFile.readline()
        if fq1:
            fq2=fastqFile.readline().strip()
            fq3=fastqFile.readline()
            fq4=fastqFile.readline()
            ID=fq1.split()[0].split('_')[0][1:]
            if scoreDict.__contains__(ID):
                continue
            leftSeq=fq2[0:ThresholdLen]
            rightSeq=fq2[::-1][0:ThresholdLen][::-1]
            seqList.append([ID,leftSeq,rightSeq])
        else:
            break
    fastqFile.close()

    outLine=100*thread
    rline=[i  for i in range(0,len(seqList),outLine)]
    number_of_entry=len(rline)
    n=0
    widgets = [ 'Test: ' , Percentage () , ' ' , Bar ( marker = RotatingMarker ( ) ) ,' ' , ETA ( ) , ' '  ]
    bar=progressbar.ProgressBar(widgets=widgets,maxval=number_of_entry).start()
    for i in rline:
        pool=Pool(processes=thread)
        scoreList=pool.map(getScoremat,seqList[i:min(i+outLine,len(seqList))])
        pool.close()
        pool.join()
        for j in scoreList:
            fout_score.write(j)
        bar.update(n)
        n+=1
    '''
    pool=Pool(processes=thread)
    pool.map(getScoremat,seqList)
    pool.close()
    pool.join()
    '''
    fout_score.close()
