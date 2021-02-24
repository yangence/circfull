import sys, os, pysam, pyfasta

def getReadInfo(read):
    rS=read.reference_start
    rE=read.reference_end
    qS=read.query_alignment_start
    qE=read.query_alignment_end
    readcigar=read.cigar
    # 5 HARD_CLIP for quary_alignment_start and quary_alignment_end(This the index just past the last base in seq that is not soft-clipped.)
    if readcigar[0][0] == 5:
        qS+=readcigar[0][1]
        qE+=readcigar[0][1]
    exonS=[rS]
    exonE=[]
    currentPos=rS-1
    i=0
    for each in readcigar:
        i=i+1
        if i == len(readcigar):
            if each[0] in (0,2):
                currentPos+=each[1]
            exonE.append(currentPos+1)
        else:
            if each[0] ==0:
                currentPos+=each[1]
            elif each[0] == 2:
                currentPos+=each[1]
            elif each[0] == 3:
                exonE.append(currentPos+1)
                currentPos+=each[1]
                exonS.append(currentPos+1)
    exonLen=sum([exonE[i]-exonS[i] for i in range(len(exonS))])
    strand="+"            
    if read.flag & 16 == 16:
        strand="-"
    return([read.query_name,read.reference_name,rS,rE,qS,qE,','.join([str(i) for i in exonS]),','.join([str(i) for i in exonE]),exonLen,strand])
def explainFL(genomeFile,outPrefix,sam):
    genome = pyfasta.Fasta(genomeFile)
    samfile=pysam.AlignmentFile(sam,"r")
    fout=open(outPrefix+'explainFL.txt','w')
    fout_nop=open(outPrefix+'explainFL_noprimary.txt','w')
    for read in samfile.fetch():
        #if read.mapping_quality < 20: 
        #   continue
        if len(read.cigar) == 0:
            continue
        readInfo=getReadInfo(read)
        leftSeq=genome.sequence({'chr': readInfo[1], 'start':readInfo[2]-1, 'stop':readInfo[2]}).upper()
        rightSeq=genome.sequence({'chr': readInfo[1], 'start':readInfo[3]+1 , 'stop':readInfo[3]+2}).upper()
        readInfo.append(leftSeq)
        readInfo.append(rightSeq)
        if read.flag & 256 >0:
            fout_nop.write('\t'.join([str(i) for i in readInfo])+"\n")
        else:
            fout.write('\t'.join([str(i) for i in readInfo])+"\n")
    samfile.close()
    fout.close()
    fout_nop.close()