import numpy as np, pandas as pd, os,sys
from multiprocessing import Pool
from progressbar import *
def SWalign(seq1,seq2):
    matchS = 3; mismatchS = -6; gapopenS = -5; gapextendS = -2
    len1=len(seq1)
    len2=len(seq2)
    scoreMat=np.zeros([len1+1,len2+1])
    gapMat_right=np.zeros([len1+1,len2+1])
    gapMat_down=np.zeros([len1+1,len2+1])
    indexMat=[[[] for i in range(len2+1)] for i in range(len1+1)]
    for i in range(len1):
        for j in range(len2):
            if seq1[i] == seq2[j]:
                baseS=scoreMat[i,j]+matchS
            else:
                baseS=scoreMat[i,j]+mismatchS
            if gapMat_right[i+1,j] == 0:
                rightS=scoreMat[i+1,j]+gapopenS
            else:
                rightS=scoreMat[i+1,j]+gapextendS
            if gapMat_down[i,j+1] == 0:
                downS=scoreMat[i,j+1]+gapopenS
            else:
                downS=scoreMat[i,j+1]+gapextendS
            arrS=[0,baseS,rightS,downS]
            arrSMax= np.max(arrS)
            scoreMat[i+1,j+1] =arrSMax
            if arrS[2] == arrSMax:
                gapMat_right[i+1,j+1]=1
                indexMat[i+1][j+1].append([i+1,j])
            if arrS[3] == arrSMax:
                gapMat_down[i+1,j+1]=1
                indexMat[i+1][j+1].append([i,j+1])
            if arrS[1] == arrSMax:
                indexMat[i+1][j+1].append([i,j])
            if arrSMax == 0:
                indexMat[i+1][j+1]=[]
            if scoreMat[i,j] == 0:
                indexMat[i][j]=[]
    indexMat=np.array(indexMat,dtype=object)
    scoreMat_max=scoreMat.max()
    max_idx=[]
    for i in range(len1):
        for j in range(len2):
            if scoreMat[i+1,j+1]==scoreMat_max:
                max_idx.append([i+1,j+1])
    return(indexMat,max_idx)
def traceBack(indexMat,max_idx):
    if len(max_idx) == 0:
        return([])
    else:
        midx=[]
        for i in max_idx:           
            new_idx=traceBack(indexMat,indexMat[i[0],i[1]])
            if len(new_idx)==0:
                all_idx=[i]
                midx.append(all_idx)
            for j in new_idx:
                all_idx=[i]
                all_idx.extend(j)
                midx.append(all_idx)
            
           
        return(midx)

def calIdentify(seq1,seq2,midx,isSeq=False):
    seq1_match=[]
    seq2_match=[]
    for i in midx:
        seq1_each=''
        seq2_each=''
        i=i[::-1]
        if len(i)<=2:
            continue
        lastSeq1=i[1][0]
        lastSeq2=i[1][1]
        seq1_each+=seq1[lastSeq1-1]
        seq2_each+=seq2[lastSeq2-1]
        for j in i[2:]:
            if j[0] == lastSeq1:
                seq1_each+='-'
            else:
                lastSeq1=j[0]
                seq1_each+=seq1[lastSeq1-1]
            if j[1] == lastSeq2:
                seq2_each+='-'
            else:
                lastSeq2=j[1]
                seq2_each+=seq2[lastSeq2-1]
        seq1_match.append(seq1_each)
        seq2_match.append(seq2_each)
    identify_num=[0]
    for i in range(len(seq1_match)):
        seq1_each=seq1_match[i]
        seq2_each=seq2_match[i]
        each_identify=0
        for j in range(len(seq1_each)):
            if seq1_each[j] == seq2_each[j]:
                each_identify+=1
        identify_num.append(each_identify)
    if isSeq:
        return([seq1_match,seq2_match,identify_num])
    return(max(identify_num))
    
def getALL(seq1,seq2,isSeq=False):
    [indexMat,max_idx]=SWalign(seq1,seq2)
    midx=traceBack(indexMat,max_idx)
    return(calIdentify(seq1,seq2,midx,isSeq))

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
